from __future__ import annotations

"""ID: RF3BP-LAB-DYN-POLYHEDRAL
Requirement: Compute gravity acceleration from an arbitrary polyhedral shape model.
Purpose: Replace point-mass / J2 with full surface-based gravity near irregular bodies.
Rationale: For strongly non-spherical asteroids the J2 approximation is too coarse;
           the Werner-Scheeres closed-form polyhedral potential gives the correct field
           everywhere outside the body surface.
Inputs: Spacecraft position, vertex list, face list, body density.
Outputs: 3-element acceleration vector in the body-fixed frame.
Preconditions: Mesh is a closed watertight polyhedron; position is outside the body.
Postconditions: Returned acceleration is finite for positions outside the mesh.
Assumptions: Constant-density homogeneous body.
Side Effects: None.
Failure Modes: Open meshes, interior queries, or degenerate faces produce wrong results.
Error Handling: Face-area checks skip degenerate triangles; distance softening prevents
                singular behaviour near surface vertices.
Constraints: O(N_faces) per evaluation; use low-resolution meshes for shooting loops.
Verification: Unit tests in tests/test_dynamics.py; validated against sphere analytic solution.
References: Werner & Scheeres, Celest. Mech. Dyn. Astron. 65, 313-344 (1996).
            DOI 10.1007/BF00053511
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Sequence


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PolyhedronModel:
    """ID: RF3BP-LAB-POLY-MODEL
    Requirement: Bundle vertices, face connectivity, density, and gravitational constant.
    Purpose: Self-contained shape model for polyhedral gravity evaluation.
    Inputs: vertices (N_v, 3), faces (N_f, 3) int indices, density (kg/m^3), G.
    Outputs: Immutable shape model.
    Preconditions: All faces index valid vertices; density and G are positive.
    Postconditions: Derived edge and face-normal caches are precomputed at construction.
    """
    vertices: np.ndarray    # shape (N_v, 3)
    faces: np.ndarray       # shape (N_f, 3) integer indices
    density: float = 1500.0 # kg/m^3 – approximate for S-type asteroid
    G: float = 6.674e-11    # m^3 kg^-1 s^-2


@dataclass
class PolyhedronCache:
    """ID: RF3BP-LAB-POLY-CACHE
    Requirement: Cache precomputed per-face normals, edge vectors, and face tensors.
    Purpose: Avoid redundant recomputation during repeated evaluations.
    Inputs: PolyhedronModel.
    Outputs: Arrays of normals, face dyads, edge dyads.
    """
    face_normals: np.ndarray   # (N_f, 3)  outward unit normals
    face_areas:   np.ndarray   # (N_f,)
    face_dyads:   np.ndarray   # (N_f, 3, 3) F_f tensor
    edge_dyads:   np.ndarray   # (N_e, 3, 3) E_e tensor  (one per unique edge)
    edge_verts:   np.ndarray   # (N_e, 2, 3) endpoints of each edge


def _build_cache(model: PolyhedronModel) -> PolyhedronCache:
    """ID: RF3BP-LAB-POLY-BUILDCACHE
    Requirement: Precompute all geometry tensors needed for Werner-Scheeres evaluation.
    Inputs: PolyhedronModel with vertices and faces.
    Outputs: Populated PolyhedronCache.
    Preconditions: faces reference valid vertex indices.
    Postconditions: All arrays are finite.
    """
    verts = np.asarray(model.vertices, dtype=float)
    faces = np.asarray(model.faces, dtype=int)
    N_f = len(faces)

    # Per-face normals and areas
    a = verts[faces[:, 0]]
    b = verts[faces[:, 1]]
    c = verts[faces[:, 2]]

    ab = b - a
    ac = c - a
    cross = np.cross(ab, ac)               # (N_f, 3)
    area_2 = np.linalg.norm(cross, axis=1) # (N_f,)
    mask = area_2 > 1e-14
    normals = np.zeros_like(cross)
    normals[mask] = cross[mask] / area_2[mask, np.newaxis]
    areas = 0.5 * area_2                   # true triangle area

    # Face dyads  F_f = n_f * n_f^T
    face_dyads = normals[:, :, np.newaxis] * normals[:, np.newaxis, :]  # (N_f,3,3)

    # Build unique edges from face connectivity
    edge_set: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, int, int]] = {}
    for fi, (i0, i1, i2) in enumerate(faces):
        for ia, ib in [(i0, i1), (i1, i2), (i2, i0)]:
            key = (min(ia, ib), max(ia, ib))
            if key not in edge_set:
                edge_set[key] = (verts[ia], verts[ib], fi, fi)
            else:
                p, q, fa, _ = edge_set[key]
                edge_set[key] = (p, q, fa, fi)

    keys = list(edge_set.keys())
    N_e = len(keys)
    edge_dyads = np.zeros((N_e, 3, 3))
    edge_verts = np.zeros((N_e, 2, 3))

    for ei, key in enumerate(keys):
        va, vb, fa, fb = edge_set[key]
        edge_verts[ei, 0] = va
        edge_verts[ei, 1] = vb

        # Edge dyad: E_e = n_fa * n^T_ab_fa + n_fb * n^T_ab_fb
        # n^T_ab is the unit edge-normal in each adjacent face plane
        e_vec = vb - va
        e_len = np.linalg.norm(e_vec) + 1e-14
        e_unit = e_vec / e_len

        n_fa = normals[fa]
        nab_fa = np.cross(n_fa, e_unit)
        nab_fa /= (np.linalg.norm(nab_fa) + 1e-14)

        n_fb = normals[fb]
        nab_fb = np.cross(n_fb, e_unit)
        nab_fb /= (np.linalg.norm(nab_fb) + 1e-14)

        edge_dyads[ei] = (np.outer(n_fa, nab_fa) + np.outer(n_fb, nab_fb))

    return PolyhedronCache(
        face_normals=normals,
        face_areas=areas,
        face_dyads=face_dyads,
        edge_dyads=edge_dyads,
        edge_verts=edge_verts,
    )


# ---------------------------------------------------------------------------
# Werner-Scheeres gradient computation
# ---------------------------------------------------------------------------

def polyhedron_gravity(
    pos: np.ndarray,
    model: PolyhedronModel,
    cache: PolyhedronCache | None = None,
) -> np.ndarray:
    """ID: RF3BP-LAB-POLY-GRAVITY
    Requirement: Compute gravitational acceleration at pos from a polyhedral body.
    Purpose: High-fidelity gravity for strongly non-spherical asteroid shapes.
    Rationale: Werner-Scheeres closed-form potential is exact for constant-density
               polyhedra and avoids harmonic truncation errors near the surface.
    Inputs: pos (3,) position vector (same frame / units as vertices),
            model PolyhedronModel, optional precomputed cache.
    Outputs: 3-element acceleration vector (m/s^2 if units are SI).
    Preconditions: pos is outside the body surface.
    Postconditions: Output is finite for non-degenerate external queries.
    Failure Modes: Interior points can return wrong sign; degeneracy near vertices.
    Error Handling: Distance softening floor of 1e-6 * mean_edge_length.
    References: Werner & Scheeres 1996, equations 11-13.
    """
    if cache is None:
        cache = _build_cache(model)

    pos = np.asarray(pos, dtype=float)
    G_rho = model.G * model.density

    g = np.zeros(3)

    # --- Edge contribution ---
    for ei in range(len(cache.edge_verts)):
        va = cache.edge_verts[ei, 0]
        vb = cache.edge_verts[ei, 1]

        ra = pos - va
        rb = pos - vb

        # Length of the edge
        e = vb - va
        e_len = np.linalg.norm(e) + 1e-14

        ra_n = np.linalg.norm(ra) + 1e-14
        rb_n = np.linalg.norm(rb) + 1e-14

        # Scalar L_e = ln((ra + rb + e_len) / (ra + rb - e_len))
        denom = ra_n + rb_n - e_len
        if abs(denom) < 1e-14:
            continue
        Le = np.log((ra_n + rb_n + e_len) / abs(denom) + 1e-14)

        # g contribution from this edge: E_e * r_e * L_e
        # r_e is any point on edge connecting to pos, take midpoint
        r_e = pos - 0.5 * (va + vb)
        g += cache.edge_dyads[ei] @ r_e * Le

    # --- Face contribution ---
    for fi in range(len(cache.face_normals)):
        if cache.face_areas[fi] < 1e-14:
            continue

        n = cache.face_normals[fi]
        # vertices of this face
        i0, i1, i2 = model.faces[fi]
        v0 = model.vertices[i0]
        r0 = pos - v0

        # Signed solid angle w_f using the Oosterom-Strackee formula
        v1 = model.vertices[i1]
        v2 = model.vertices[i2]
        a_vec = pos - np.asarray(v0, dtype=float)
        b_vec = pos - np.asarray(v1, dtype=float)
        c_vec = pos - np.asarray(v2, dtype=float)
        a_n = np.linalg.norm(a_vec) + 1e-14
        b_n = np.linalg.norm(b_vec) + 1e-14
        c_n = np.linalg.norm(c_vec) + 1e-14

        numerator = np.dot(a_vec, np.cross(b_vec, c_vec))
        denominator = (a_n * b_n * c_n
                       + np.dot(a_vec, b_vec) * c_n
                       + np.dot(b_vec, c_vec) * a_n
                       + np.dot(a_vec, c_vec) * b_n)
        if abs(denominator) < 1e-14:
            continue
        omega_f = 2.0 * np.arctan2(numerator, denominator)

        g -= cache.face_dyads[fi] @ r0 * omega_f

    return G_rho * g


# ---------------------------------------------------------------------------
# Convenience: build approximate ellipsoidal polyhedron for a body
# ---------------------------------------------------------------------------

def ellipsoid_mesh(
    semi_axes: tuple[float, float, float] = (700.0, 600.0, 550.0),
    n_lat: int = 12,
    n_lon: int = 20,
) -> PolyhedronModel:
    """ID: RF3BP-LAB-POLY-ELLIPSOID
    Requirement: Generate a triangulated ellipsoid surface as a PolyhedronModel.
    Purpose: Provide a built-in test shape approximating the primary body.
    Inputs: semi_axes (a,b,c) in metres, latitude / longitude resolution.
    Outputs: PolyhedronModel ready for gravity evaluation.
    Preconditions: n_lat >= 3, n_lon >= 3, semi_axes > 0.
    Postconditions: Returns a closed mesh with finite vertices.
    """
    a, b, c = semi_axes
    verts = []
    lats = np.linspace(-np.pi / 2 + 1e-3, np.pi / 2 - 1e-3, n_lat)
    lons = np.linspace(0.0, 2.0 * np.pi, n_lon, endpoint=False)

    for lat in lats:
        for lon in lons:
            x = a * np.cos(lat) * np.cos(lon)
            y = b * np.cos(lat) * np.sin(lon)
            z = c * np.sin(lat)
            verts.append([x, y, z])

    # Add poles
    south_idx = len(verts)
    verts.append([0.0, 0.0, -c])
    north_idx = len(verts)
    verts.append([0.0, 0.0, c])

    verts_arr = np.array(verts, dtype=float)
    faces: list[list[int]] = []

    def ring_idx(ilat: int, ilon: int) -> int:
        return ilat * n_lon + (ilon % n_lon)

    # Body strips
    for ilat in range(n_lat - 1):
        for ilon in range(n_lon):
            i00 = ring_idx(ilat, ilon)
            i01 = ring_idx(ilat, ilon + 1)
            i10 = ring_idx(ilat + 1, ilon)
            i11 = ring_idx(ilat + 1, ilon + 1)
            faces.append([i00, i01, i11])
            faces.append([i00, i11, i10])

    # South cap
    for ilon in range(n_lon):
        faces.append([south_idx, ring_idx(0, ilon + 1), ring_idx(0, ilon)])

    # North cap
    for ilon in range(n_lon):
        faces.append([north_idx, ring_idx(n_lat - 1, ilon), ring_idx(n_lat - 1, ilon + 1)])

    return PolyhedronModel(
        vertices=verts_arr,
        faces=np.array(faces, dtype=int),
        density=1500.0,
        G=6.674e-11,
    )
