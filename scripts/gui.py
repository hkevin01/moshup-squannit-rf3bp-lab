from __future__ import annotations

# ---------------------------------------------------------------------------
# Ensure src/ is on sys.path when the GUI is launched without pip install.
# ---------------------------------------------------------------------------
import pathlib, sys
_repo_root = pathlib.Path(__file__).resolve().parent.parent
_src = _repo_root / "src"
if str(_src) not in sys.path:
    sys.path.insert(0, str(_src))
# ---------------------------------------------------------------------------

"""RF3BP Lab - PyQt5 graphical interface.

Allows the user to:
  - Edit all SystemParams fields with labelled spin-boxes
  - Choose where figures are saved
  - Run the full demo pipeline in a background thread with live log output
  - View each generated PNG in a tabbed image viewer
  - Open the output folder in the system file manager
"""

import os
import subprocess
import threading

import matplotlib
matplotlib.use("Agg")  # must come before any pyplot import
import matplotlib.pyplot as plt
import numpy as np

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QPixmap
from PyQt5.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------

class DemoWorker(QThread):
    """Run the demo pipeline on a background thread and emit progress signals."""

    log_line = pyqtSignal(str)       # one log line
    figure_ready = pyqtSignal(str)   # absolute path to a saved PNG
    finished_ok = pyqtSignal(dict)   # summary_metrics dict
    finished_err = pyqtSignal(str)   # error message

    def __init__(self, params_dict: dict, out_dir: str) -> None:
        super().__init__()
        self._params = params_dict
        self._out_dir = out_dir

    def run(self) -> None:
        try:
            self._do_run()
        except Exception as exc:
            import traceback
            self.finished_err.emit(traceback.format_exc())

    def _emit(self, msg: str) -> None:
        self.log_line.emit(msg)

    def _do_run(self) -> None:
        from rf3bp_lab.dynamics.models import (
            FidelityWeights, compare_cr3bp_rf3bp, cr3bp_rhs,
            propagate, rf3bp_breakdown, rf3bp_pulsating_rhs,
        )
        from rf3bp_lab.dynamics.params import SystemParams
        from rf3bp_lab.shooting.hierarchical import HierarchicalShooter
        from rf3bp_lab.utils.plotting import (
            plot_model_gap, plot_perturbation_norms,
            plot_result_dashboard, plot_stage_convergence, plot_trajectory,
        )
        import json

        os.makedirs(self._out_dir, exist_ok=True)

        p = SystemParams(**self._params)
        seed = np.array([0.55, 0.0, 0.0, 0.0, 0.42, 0.02], dtype=float)

        self._emit("Starting hierarchical shooter ...")
        shooter = HierarchicalShooter(p)
        result = shooter.solve(seed)

        s0 = result["state0"]
        period = float(result["period"])
        stage_data = result["stages"]

        self._emit(f"  period          = {period:.6f}")
        self._emit(f"  residual_norm   = {result['residual_norm']:.3e}")
        self._emit(f"  stage1_cost     = {result['stage1_cost']:.3e}")
        self._emit(f"  final_cost      = {result['stage2_cost']:.3e}")
        self._emit(f"  stages          = {len(stage_data)}")

        t_eval = np.linspace(0.0, period * 6.0, 4000)
        self._emit("Propagating CR3BP orbit ...")
        sol_cr3bp = propagate(cr3bp_rhs, s0, (0.0, t_eval[-1]), p, t_eval=t_eval)
        self._emit("Propagating RF3BP orbit ...")
        sol_rf3bp = propagate(rf3bp_pulsating_rhs, s0, (0.0, t_eval[-1]), p, t_eval=t_eval)

        def _save(fig_name: str) -> str:
            path = os.path.join(self._out_dir, fig_name)
            plt.savefig(path, dpi=150, bbox_inches="tight")
            plt.close()
            self.figure_ready.emit(path)
            self._emit(f"  saved {fig_name}")
            return path

        self._emit("Generating trajectory plots ...")
        plot_trajectory(sol_cr3bp.y[:3], "CR3BP Reference Orbit")
        _save("trajectory_cr3bp.png")

        plot_trajectory(sol_rf3bp.y[:3], "RF3BP Orbit - Pulsation + Perturbations")
        _save("trajectory_rf3bp.png")

        self._emit("Computing perturbation norms ...")
        keys = ["nonspherical", "pulsation", "solar_gravity", "srp"]
        history = {k: np.zeros_like(sol_rf3bp.t) for k in keys}
        for i in range(sol_rf3bp.t.size):
            b_full = rf3bp_breakdown(float(sol_rf3bp.t[i]), sol_rf3bp.y[:, i], p)
            b_nj2 = rf3bp_breakdown(float(sol_rf3bp.t[i]), sol_rf3bp.y[:, i], p,
                                    FidelityWeights(1.0, 0.0, 1.0, 1.0))
            history["nonspherical"][i] = np.linalg.norm(
                (b_full.grav_primary + b_full.grav_secondary)
                - (b_nj2.grav_primary + b_nj2.grav_secondary))
            history["pulsation"][i] = np.linalg.norm(b_full.pulsation)
            history["solar_gravity"][i] = np.linalg.norm(b_full.solar_gravity)
            history["srp"][i] = np.linalg.norm(b_full.srp)

        plot_perturbation_norms(sol_rf3bp.t, history)
        _save("perturbation_norms.png")

        stage_labels = [s.label for s in stage_data]
        stage_residuals = np.array([s.residual_norm for s in stage_data])
        plot_stage_convergence(stage_labels, stage_residuals)
        _save("continuation_convergence.png")

        self._emit("Computing RF3BP vs CR3BP model gap ...")
        gap_abs = np.zeros_like(sol_rf3bp.t)
        gap_rel = np.zeros_like(sol_rf3bp.t)
        for i in range(sol_rf3bp.t.size):
            g = compare_cr3bp_rf3bp(float(sol_rf3bp.t[i]), sol_rf3bp.y[:, i], p)
            gap_abs[i] = g.delta_norm
            gap_rel[i] = g.relative_gap

        plot_model_gap(sol_rf3bp.t, gap_abs, gap_rel)
        _save("model_gap_cr3bp_vs_rf3bp.png")

        perturbation_peaks = {k: float(np.max(v)) for k, v in history.items()}
        summary = {
            "period": period,
            "residual_norm": float(result["residual_norm"]),
            "stage1_cost": float(result["stage1_cost"]),
            "final_cost": float(result["stage2_cost"]),
            "max_abs_gap": float(np.max(gap_abs)),
            "max_rel_gap": float(np.max(gap_rel)),
            "mean_abs_gap": float(np.mean(gap_abs)),
            "mean_rel_gap": float(np.mean(gap_rel)),
            "n_steps": int(sol_rf3bp.t.size),
        }

        plot_result_dashboard(summary, perturbation_peaks,
                              stage_labels=stage_labels, stage_residuals=stage_residuals)
        _save("result_snapshot_dashboard.png")

        results_dir = str(_repo_root / "docs" / "results")
        os.makedirs(results_dir, exist_ok=True)
        json_path = os.path.join(results_dir, "latest_demo_metrics.json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump({"summary_metrics": summary,
                       "perturbation_peaks": perturbation_peaks,
                       "stage_residuals": [float(s.residual_norm) for s in stage_data],
                       "stage_labels": stage_labels}, f, indent=2, sort_keys=True)
        self._emit(f"Metrics saved -> {json_path}")
        self._emit("Done.")
        self.finished_ok.emit(summary)


# ---------------------------------------------------------------------------
# Parameter form
# ---------------------------------------------------------------------------

_PARAM_META = [
    ("mu",                 "Mass ratio mu",              0.001,  0.1,   5, 0.001),
    ("r12_mean_m",         "Mean separation (m)",        100.0,  1e5,   1, 10.0),
    ("pulsation_e",        "Pulsation eccentricity e",   0.0,    0.5,   4, 0.01),
    ("pulsation_nu",       "Pulsation frequency nu",     0.01,   5.0,   4, 0.01),
    ("j2_primary",         "J2 primary",                 0.0,    0.5,   4, 0.005),
    ("j2_secondary",       "J2 secondary",               0.0,    0.5,   4, 0.005),
    ("r_primary_m",        "Primary radius (m)",         10.0,   5000.0,1, 10.0),
    ("r_secondary_m",      "Secondary radius (m)",       5.0,    2000.0,1, 5.0),
    ("sun_mu_scaled",      "Solar mu (scaled)",          1e-6,   1e-2,  8, 1e-5),
    ("sun_distance_scaled","Solar distance (scaled)",    100.0,  1e5,   1, 10.0),
    ("srp_accel_scaled",   "SRP accel (scaled)",         1e-9,   1e-3,  9, 1e-7),
    ("omega0",             "Frame rotation omega0",      0.1,    10.0,  4, 0.1),
]

def _default_params() -> dict:
    from rf3bp_lab.dynamics.params import SystemParams
    p = SystemParams()
    return {f: getattr(p, f) for f, *_ in _PARAM_META}


class ParamPanel(QGroupBox):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__("System Parameters", parent)
        form = QFormLayout(self)
        self._spinboxes: dict[str, QDoubleSpinBox] = {}
        defaults = _default_params()
        for name, label, lo, hi, decimals, step in _PARAM_META:
            sb = QDoubleSpinBox()
            sb.setRange(lo, hi)
            sb.setDecimals(decimals)
            sb.setSingleStep(step)
            sb.setValue(defaults[name])
            sb.setMinimumWidth(130)
            self._spinboxes[name] = sb
            form.addRow(label + ":", sb)

        btn_reset = QPushButton("Reset to defaults")
        btn_reset.clicked.connect(self._reset)
        form.addRow(btn_reset)

    def _reset(self) -> None:
        defaults = _default_params()
        for name, sb in self._spinboxes.items():
            sb.setValue(defaults[name])

    def get_params(self) -> dict:
        return {name: sb.value() for name, sb in self._spinboxes.items()}


# ---------------------------------------------------------------------------
# Image viewer tab widget
# ---------------------------------------------------------------------------

class ImageViewer(QScrollArea):
    def __init__(self) -> None:
        super().__init__()
        self._label = QLabel()
        self._label.setAlignment(Qt.AlignCenter)
        self.setWidget(self._label)
        self.setWidgetResizable(True)
        self._path: str | None = None

    def load(self, path: str) -> None:
        self._path = path
        self._reload()

    def _reload(self) -> None:
        if not self._path or not os.path.exists(self._path):
            return
        pix = QPixmap(self._path)
        available = self.size()
        scaled = pix.scaled(available.width() - 10, available.height() - 10,
                            Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self._label.setPixmap(scaled)

    def resizeEvent(self, event):  # type: ignore[override]
        super().resizeEvent(event)
        self._reload()


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

_FIG_NAMES = [
    ("trajectory_cr3bp.png",          "CR3BP Orbit"),
    ("trajectory_rf3bp.png",          "RF3BP Orbit"),
    ("perturbation_norms.png",        "Perturbations"),
    ("continuation_convergence.png",  "Convergence"),
    ("model_gap_cr3bp_vs_rf3bp.png",  "Model Gap"),
    ("result_snapshot_dashboard.png", "Dashboard"),
]


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("RF3BP Lab - Moshup-Squannit")
        self.resize(1280, 820)
        self._worker: DemoWorker | None = None
        self._out_dir = str(_repo_root / "docs" / "figures")

        # ---- central splitter ----
        splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(splitter)

        # ---- LEFT panel ----
        left = QWidget()
        left.setMaximumWidth(320)
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(6, 6, 6, 6)

        self._param_panel = ParamPanel()
        scroll = QScrollArea()
        scroll.setWidget(self._param_panel)
        scroll.setWidgetResizable(True)
        left_lay.addWidget(scroll, stretch=1)

        # output directory
        dir_grp = QGroupBox("Output directory")
        dir_lay = QHBoxLayout(dir_grp)
        self._dir_label = QLabel(self._out_dir)
        self._dir_label.setWordWrap(True)
        btn_browse = QPushButton("Browse...")
        btn_browse.setFixedWidth(80)
        btn_browse.clicked.connect(self._browse_dir)
        dir_lay.addWidget(self._dir_label, stretch=1)
        dir_lay.addWidget(btn_browse)
        left_lay.addWidget(dir_grp)

        # run / open buttons
        btn_row = QHBoxLayout()
        self._btn_run = QPushButton("Run Demo")
        self._btn_run.setFixedHeight(36)
        self._btn_run.setFont(QFont("", 11, QFont.Bold))
        self._btn_run.clicked.connect(self._run)
        self._btn_open = QPushButton("Open Figures Folder")
        self._btn_open.setFixedHeight(36)
        self._btn_open.clicked.connect(self._open_dir)
        btn_row.addWidget(self._btn_run)
        btn_row.addWidget(self._btn_open)
        left_lay.addLayout(btn_row)

        self._progress = QProgressBar()
        self._progress.setRange(0, 0)  # indeterminate
        self._progress.setVisible(False)
        left_lay.addWidget(self._progress)

        splitter.addWidget(left)

        # ---- RIGHT panel ----
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(4, 4, 4, 4)
        right_lay.setSpacing(0)

        # Vertical splitter so the user can drag the divider between
        # the figure tabs and the log panel - they can never overlap.
        right_splitter = QSplitter(Qt.Vertical)

        self._tabs = QTabWidget()
        self._viewers: dict[str, ImageViewer] = {}
        for fname, tab_label in _FIG_NAMES:
            viewer = ImageViewer()
            self._viewers[fname] = viewer
            self._tabs.addTab(viewer, tab_label)
        right_splitter.addWidget(self._tabs)

        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        self._log.setMaximumBlockCount(2000)
        self._log.setFont(QFont("Monospace", 9))
        self._log.setMinimumHeight(80)
        right_splitter.addWidget(self._log)

        # Give tabs ~75 % of the height, log ~25 %
        right_splitter.setStretchFactor(0, 3)
        right_splitter.setStretchFactor(1, 1)
        right_splitter.setSizes([580, 180])

        right_lay.addWidget(right_splitter)

        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        # pre-load any existing figures
        self._reload_existing_figures()

    # ------------------------------------------------------------------
    def _browse_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Choose output directory", self._out_dir)
        if d:
            self._out_dir = d
            self._dir_label.setText(d)

    def _open_dir(self) -> None:
        path = self._out_dir
        if not os.path.isdir(path):
            QMessageBox.warning(self, "Not found", f"Directory does not exist:\n{path}")
            return
        if sys.platform == "win32":
            os.startfile(path)  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])

    def _reload_existing_figures(self) -> None:
        for fname, _ in _FIG_NAMES:
            full = os.path.join(self._out_dir, fname)
            if os.path.exists(full):
                self._viewers[fname].load(full)

    def _run(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        self._log.clear()
        self._btn_run.setEnabled(False)
        self._progress.setVisible(True)
        self._log.appendPlainText("=" * 60)
        self._log.appendPlainText(f"Output dir: {self._out_dir}")

        params = self._param_panel.get_params()
        self._worker = DemoWorker(params, self._out_dir)
        self._worker.log_line.connect(self._on_log)
        self._worker.figure_ready.connect(self._on_figure)
        self._worker.finished_ok.connect(self._on_done)
        self._worker.finished_err.connect(self._on_error)
        self._worker.start()

    def _on_log(self, msg: str) -> None:
        self._log.appendPlainText(msg)
        self._log.verticalScrollBar().setValue(self._log.verticalScrollBar().maximum())

    def _on_figure(self, path: str) -> None:
        fname = os.path.basename(path)
        if fname in self._viewers:
            self._viewers[fname].load(path)
            # Switch to that tab so the user sees results arriving live
            for i, (f, _) in enumerate(_FIG_NAMES):
                if f == fname:
                    self._tabs.setCurrentIndex(i)
                    break

    def _on_done(self, summary: dict) -> None:
        self._progress.setVisible(False)
        self._btn_run.setEnabled(True)
        self._log.appendPlainText("")
        self._log.appendPlainText("Results summary:")
        for k, v in summary.items():
            self._log.appendPlainText(f"  {k:<22} = {v}")
        self._log.appendPlainText("=" * 60)

    def _on_error(self, tb: str) -> None:
        self._progress.setVisible(False)
        self._btn_run.setEnabled(True)
        self._log.appendPlainText("\n[ERROR]\n" + tb)
        QMessageBox.critical(self, "Run failed",
                             "The demo raised an exception.\nSee the log panel for details.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
