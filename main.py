"""
AI-Guided Molecular Docking Pipeline
=====================================
main.py — Tkinter GUI + orchestration entry point

Coordinates:
  - File selection (receptor PDB + multi-ligand SDF)
  - Settings dialog (gear button)
  - Background pipeline execution with live log output
  - Progress tracking per stage

Dependencies (besides stdlib):
  pip install openai numpy scipy rdkit
  + PyMOL, Meeko, AutoDock Vina installed separately (paths set in Settings)
"""

from __future__ import annotations

import os
import platform
import queue
import subprocess
import threading
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk
from typing import Any, Callable, Optional

# ── Local pipeline modules ──────────────────────────────────────────────────
try:
    from pipeline import run_pipeline
    from settings import SettingsDialog, load_config, validate_config
except ImportError as _e:
    # Allow the GUI to open even if sibling modules aren't ready yet
    _MISSING_MODULES: Optional[str] = str(_e)
else:
    _MISSING_MODULES = None


# ════════════════════════════════════════════════════════════════════════════
#  Main Application Window
# ════════════════════════════════════════════════════════════════════════════

class DockingApp:
    """Root Tkinter window for the docking pipeline."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("AI GUIDED DOCKING PROTOCOL (AGDP)")
        self.root.geometry("980x740")
        self.root.minsize(820, 620)

        # State Variables
        self.config: dict[str, Any] = load_config() if not _MISSING_MODULES else {}
        self.log_queue: queue.Queue[tuple[str, Any, str]] = queue.Queue()
        self.is_running: bool = False
        self._stop_requested: bool = False

        self._setup_styles()
        self._build_ui()
        self._poll()  # Start GUI ↔ thread message pump

        # Warn if sibling modules couldn't be imported
        if _MISSING_MODULES:
            self._write_log(
                f"⚠  Could not import pipeline modules: {_MISSING_MODULES}\n"
                "   Make sure settings.py, pipeline.py, ai_identify.py, "
                "receptor_prep.py,\n"
                "   ligand_prep.py, and docking.py are in the same folder.",
                "warning",
            )

    # ── Styles ──────────────────────────────────────────────────────────────
    def _setup_styles(self) -> None:
        s = ttk.Style()
        try:
            s.theme_use("clam")
        except tk.TclError:
            pass  # Fall back to default theme on platforms lacking 'clam'

        bg = "#ecf0f1"
        s.configure("TFrame", background=bg)
        s.configure("TLabelframe", background=bg)
        s.configure(
            "TLabelframe.Label",
            font=("Segoe UI", 9, "bold"),
            foreground="#2c3e50",
            background=bg,
        )
        s.configure(
            "Header.TLabel",
            font=("Segoe UI", 13, "bold"),
            foreground="#2c3e50",
            background=bg,
        )
        s.configure(
            "Sub.TLabel",
            font=("Segoe UI", 8),
            foreground="#7f8c8d",
            background=bg,
        )
        s.configure("Run.TButton", font=("Segoe UI", 10, "bold"), padding=7)
        s.configure("Small.TButton", font=("Segoe UI", 9), padding=4)
        s.configure("TEntry", fieldbackground="white")
        s.configure("TSpinbox", fieldbackground="white")

    # ── UI Construction ─────────────────────────────────────────────────────
    def _build_ui(self) -> None:
        self.root.configure(bg="#ecf0f1")
        outer = ttk.Frame(self.root, padding="12")
        outer.pack(fill=tk.BOTH, expand=True)

        self._build_header(outer)
        self._build_inputs(outer)
        self._build_status_bar(outer)
        self._build_log(outer)
        self._build_buttons(outer)

    # ── Header row ──────────────────────────────────────────────────────────
    def _build_header(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=(0, 8))

        ttk.Label(
            row,
            text="  AI GUIDED DOCKING PROTOCOL (AGDP)",
            style="Header.TLabel",
        ).pack(side=tk.LEFT)

        ttk.Button(
            row,
            text="⚙  Settings",
            style="Small.TButton",
            command=self._open_settings,
        ).pack(side=tk.RIGHT)

        ttk.Label(
            row,
            text="AutoDock Vina  ·  Meeko  ·  PyMOL  ·  AI ",
            style="Sub.TLabel",
        ).pack(side=tk.RIGHT, padx=12)

    # ── Input-file section ───────────────────────────────────────────────────
    def _build_inputs(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Input Files", padding="10")
        frame.pack(fill=tk.X, pady=(0, 8))
        frame.columnconfigure(1, weight=1)

        # Receptor row
        ttk.Label(frame, text="Receptor PDB:").grid(
            row=0, column=0, sticky=tk.W, padx=(0, 8), pady=(0, 6)
        )

        self.receptor_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.receptor_var, state="readonly").grid(
            row=0, column=1, sticky=tk.EW, padx=(0, 6)
        )

        ttk.Button(
            frame,
            text="Browse…",
            style="Small.TButton",
            command=self._browse_receptor,
        ).grid(row=0, column=2)

        # Ligand rows
        ttk.Label(frame, text="Ligand SDF(s):").grid(
            row=1, column=0, sticky=tk.NW, padx=(0, 8), pady=(4, 0)
        )

        lb_wrap = ttk.Frame(frame)
        lb_wrap.grid(row=1, column=1, sticky=tk.NSEW, padx=(0, 6), pady=(4, 0))
        frame.rowconfigure(1, weight=1)

        self.lig_listbox = tk.Listbox(
            lb_wrap,
            height=5,
            selectmode=tk.EXTENDED,
            bg="white",
            relief="solid",
            bd=1,
            font=("Segoe UI", 9),
            activestyle="dotbox",
        )
        self.lig_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        ysb = ttk.Scrollbar(lb_wrap, orient=tk.VERTICAL, command=self.lig_listbox.yview)
        ysb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lig_listbox.configure(yscrollcommand=ysb.set)

        # Ligand action buttons
        lig_btns = ttk.Frame(frame)
        lig_btns.grid(row=1, column=2, sticky=tk.N, pady=(4, 0))

        buttons: list[tuple[str, Callable[[], None]]] = [
            ("Add…", self._add_ligands),
            ("Remove", self._remove_ligands),
            ("Clear", self._clear_ligands),
        ]
        for label, cmd in buttons:
            ttk.Button(
                lig_btns, text=label, style="Small.TButton", command=cmd, width=8
            ).pack(pady=(0, 4))

        # Ligand count badge
        self.lig_count_var = tk.StringVar(value="0 file(s)")
        ttk.Label(
            lig_btns,
            textvariable=self.lig_count_var,
            foreground="gray",
            font=("Segoe UI", 8),
        ).pack()

    # ── Status bar / progress ────────────────────────────────────────────────
    def _build_status_bar(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=(0, 4))

        self.status_var = tk.StringVar(
            value="Ready — configure settings then select files."
        )
        ttk.Label(row, textvariable=self.status_var, font=("Segoe UI", 9)).pack(
            side=tk.LEFT
        )

        self.pct_var = tk.StringVar(value="0%")
        ttk.Label(row, textvariable=self.pct_var, font=("Segoe UI", 9), width=5).pack(
            side=tk.RIGHT
        )

        self.progress_var = tk.DoubleVar(value=0)
        ttk.Progressbar(
            row,
            variable=self.progress_var,
            maximum=100,
            length=240,
            mode="determinate",
        ).pack(side=tk.RIGHT, padx=(0, 4))

    # ── Log panel ────────────────────────────────────────────────────────────
    def _build_log(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Pipeline Log", padding="4")
        frame.pack(fill=tk.BOTH, expand=True, pady=(0, 8))

        self.log_widget = scrolledtext.ScrolledText(
            frame,
            wrap=tk.WORD,
            state="disabled",
            bg="#0d1117",
            fg="#c9d1d9",
            font=("Consolas", 9),
            height=22,
            insertbackground="white",
            selectbackground="#264f78",
            relief="flat",
            padx=6,
            pady=6,
        )
        self.log_widget.pack(fill=tk.BOTH, expand=True)

        # Colour tags
        _tags: dict[str, tuple[str, str]] = {
            "info": ("#79c0ff", "normal"),
            "success": ("#56d364", "normal"),
            "warning": ("#e3b341", "normal"),
            "error": ("#ff7b72", "normal"),
            "header": ("#d2a8ff", "bold"),
            "dim": ("#484f58", "normal"),
            "bold": ("#f0f6fc", "bold"),
        }
        for tag, (color, weight) in _tags.items():
            self.log_widget.tag_config(
                tag, foreground=color, font=("Consolas", 9, weight)
            )

    # ── Bottom button row ────────────────────────────────────────────────────
    def _build_buttons(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X)

        self.run_btn = ttk.Button(
            row,
            text="▶   Run Docking Pipeline",
            style="Run.TButton",
            command=self._start_pipeline,
        )
        self.run_btn.pack(side=tk.LEFT, padx=(0, 8))

        self.stop_btn = ttk.Button(
            row,
            text="■  Stop",
            style="Small.TButton",
            command=self._stop_pipeline,
            state="disabled",
        )
        self.stop_btn.pack(side=tk.LEFT)

        ttk.Button(
            row, text="Clear Log", style="Small.TButton", command=self._clear_log
        ).pack(side=tk.RIGHT)

        ttk.Button(
            row,
            text="Open Output Folder",
            style="Small.TButton",
            command=self._open_output_folder,
        ).pack(side=tk.RIGHT, padx=(0, 8))

    # ════════════════════════════════════════════════════════════════════════
    #  File picker callbacks
    # ════════════════════════════════════════════════════════════════════════

    def _browse_receptor(self) -> None:
        path = filedialog.askopenfilename(
            title="Select Receptor PDB File",
            parent=self.root,
            filetypes=[("PDB files", "*.pdb *.PDB"), ("All files", "*")],
        )
        if path:
            self.receptor_var.set(path)
            self._write_log(f"Receptor selected: {Path(path).name}", "dim")

    def _add_ligands(self) -> None:
        paths = filedialog.askopenfilenames(
            title="Select Ligand SDF File(s)  [hold Ctrl/Shift for multi-select]",
            parent=self.root,
            filetypes=[
                ("SDF files", "*.sdf *.SDF"),
                ("MOL files", "*.mol"),
                ("All files", "*"),
            ],
        )
        existing = set(self.lig_listbox.get(0, tk.END))
        added = 0
        for p in paths:
            if p not in existing:
                self.lig_listbox.insert(tk.END, p)
                existing.add(p)
                added += 1
        if added:
            self._refresh_lig_count()
            self._write_log(
                f"Added {added} ligand(s). Total: {self.lig_listbox.size()}", "dim"
            )

    def _remove_ligands(self) -> None:
        sel = self.lig_listbox.curselection()
        for idx in reversed(sel):
            self.lig_listbox.delete(idx)
        self._refresh_lig_count()

    def _clear_ligands(self) -> None:
        self.lig_listbox.delete(0, tk.END)
        self._refresh_lig_count()

    def _refresh_lig_count(self) -> None:
        n = self.lig_listbox.size()
        self.lig_count_var.set(f"{n} file(s)")

    # ════════════════════════════════════════════════════════════════════════
    #  Settings
    # ════════════════════════════════════════════════════════════════════════

    def _open_settings(self) -> None:
        if _MISSING_MODULES:
            messagebox.showerror(
                "Import Error",
                f"Cannot open Settings — module import failed:\n{_MISSING_MODULES}",
            )
            return
        SettingsDialog(self.root, self.config)
        self.config = load_config()  # reload after dialog closes
        self._write_log("Settings updated.", "dim")

    # ════════════════════════════════════════════════════════════════════════
    #  Pipeline execution
    # ════════════════════════════════════════════════════════════════════════

    def _start_pipeline(self) -> None:
        if _MISSING_MODULES:
            messagebox.showerror(
                "Import Error",
                f"Pipeline modules not loaded:\n{_MISSING_MODULES}",
            )
            return

        receptor = self.receptor_var.get().strip()
        ligands = list(self.lig_listbox.get(0, tk.END))

        # Input validation
        if not receptor:
            messagebox.showerror("Missing Input", "Please select a Receptor PDB file.")
            return
        if not Path(receptor).exists():
            messagebox.showerror(
                "File Not Found", f"Receptor file does not exist:\n{receptor}"
            )
            return
        if not ligands:
            messagebox.showerror(
                "Missing Input", "Please add at least one Ligand SDF file."
            )
            return
        for lf in ligands:
            if not Path(lf).exists():
                messagebox.showerror("File Not Found", f"Ligand file not found:\n{lf}")
                return

        errors = validate_config(self.config)
        if errors:
            messagebox.showerror(
                "Configuration Incomplete",
                "Please fix the following in ⚙ Settings:\n\n"
                + "\n".join(f"  •  {e}" for e in errors),
            )
            return

        # Reset UI
        self._clear_log()
        self.progress_var.set(0)
        self.pct_var.set("0%")
        self.status_var.set("Starting pipeline…")
        self.is_running = True
        self._stop_requested = False
        self.run_btn.configure(state="disabled")
        self.stop_btn.configure(state="normal")

        # Banner
        self._write_log("═" * 60, "header")
        self._write_log("  AI GUIDED DOCKING PROTOCOL (AGDP)", "header")
        self._write_log(
            f"  Started: {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}", "header"
        )
        self._write_log("═" * 60, "header")
        self._write_log(f"Receptor : {Path(receptor).name}", "info")
        self._write_log(
            f"Ligands  : {len(ligands)} file(s)  "
            + "  |  ".join(Path(l).name for l in ligands[:3])
            + (" …" if len(ligands) > 3 else ""),
            "info",
        )
        self._write_log("", "dim")

        # Launch background thread
        t = threading.Thread(
            target=self._pipeline_thread,
            args=(receptor, ligands),
            daemon=True,
        )
        t.start()

    def _stop_pipeline(self) -> None:
        self._stop_requested = True
        self._queue(
            "log", "⏹  Stop requested — will halt after the current step.", "warning"
        )

    # ── Background worker ────────────────────────────────────────────────────
    def _pipeline_thread(self, receptor: str, ligands: list[str]) -> None:
        try:
            run_pipeline(
                receptor_file=receptor,
                ligand_files=ligands,
                config=self.config,
                log_cb=lambda msg, lvl="info": self._queue("log", msg, lvl),
                progress_cb=self._queue_progress,
                status_cb=lambda s: self._queue("status", s),
                stop_flag=lambda: self._stop_requested,
            )
        except Exception as exc:
            import traceback

            self._queue("log", f"❌  Fatal error: {exc}", "error")
            self._queue("log", traceback.format_exc(), "dim")
        finally:
            self._queue("done", None)

    # ════════════════════════════════════════════════════════════════════════
    #  Thread → GUI message pump (queue)
    # ════════════════════════════════════════════════════════════════════════

    def _queue(self, kind: str, value: Any = None, level: str = "info") -> None:
        self.log_queue.put((kind, value, level))

    def _queue_progress(self, pct: float) -> None:
        self.log_queue.put(("progress", pct, "info"))

    def _poll(self) -> None:
        """Drain the queue and update GUI — runs periodically in the main thread."""
        try:
            while True:
                kind, value, level = self.log_queue.get_nowait()

                if kind == "log":
                    self._write_log(str(value), level)

                elif kind == "progress":
                    val = float(value)
                    self.progress_var.set(val)
                    self.pct_var.set(f"{int(val)}%")

                elif kind == "status":
                    self.status_var.set(str(value))

                elif kind == "done":
                    self.is_running = False
                    self.run_btn.configure(state="normal")
                    self.stop_btn.configure(state="disabled")
                    if not self._stop_requested:
                        self.status_var.set("Pipeline complete ✓")
                        self.progress_var.set(100)
                        self.pct_var.set("100%")
                    else:
                        self.status_var.set("Stopped by user.")

        except queue.Empty:
            pass
        finally:
            self.root.after(80, self._poll)

    # ════════════════════════════════════════════════════════════════════════
    #  Log helpers
    # ════════════════════════════════════════════════════════════════════════

    def _write_log(self, message: str, level: str = "info") -> None:
        self.log_widget.configure(state="normal")
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_widget.insert(tk.END, f"[{ts}]  {message}\n", level)
        self.log_widget.see(tk.END)
        self.log_widget.configure(state="disabled")

    def _clear_log(self) -> None:
        self.log_widget.configure(state="normal")
        self.log_widget.delete("1.0", tk.END)
        self.log_widget.configure(state="disabled")

    # ════════════════════════════════════════════════════════════════════════
    #  Utility actions
    # ════════════════════════════════════════════════════════════════════════

    def _open_output_folder(self) -> None:
        out_dir = self.config.get("output_dir", "")
        if not out_dir or not Path(out_dir).exists():
            messagebox.showinfo(
                "Output Folder",
                "Output directory not set or doesn't exist yet.\nSet it via ⚙ Settings.",
            )
            return
        _open_folder(out_dir)


# ════════════════════════════════════════════════════════════════════════════
#  Cross-platform folder opener
# ════════════════════════════════════════════════════════════════════════════

def _open_folder(path: str) -> None:
    p = platform.system()
    try:
        if p == "Windows":
            os.startfile(path)  # type: ignore[attr-defined]
        elif p == "Darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception as e:
        messagebox.showwarning("Open Folder", f"Could not open folder:\n{e}")


# ════════════════════════════════════════════════════════════════════════════
#  Entry point
# ════════════════════════════════════════════════════════════════════════════

def main() -> None:
    root = tk.Tk()

    # App icon (optional — silently ignored if not available)
    try:
        icon = tk.PhotoImage(
            data=(
                "R0lGODlhEAAQAIAAAAAAAP///yH5BAEAAAAALAAAAAAQABAAAAIghI+py+0Po5y02ouz"
                "3rz7D4biSJbmiabqyrbuC8fyTAUAOw=="
            )
        )
        root.iconphoto(True, icon)
    except Exception:
        pass

    _app = DockingApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()