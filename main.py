"""
Molecular Docking Pipeline
=====================================
main.py — Tkinter GUI + orchestration entry point

Coordinates:
  - File selection (receptor structure PDB/mmCIF + multi-ligand SDF)
  - Settings dialog (gear button)
  - Background pipeline execution with live log output
  - Progress tracking per stage

Dependencies (besides stdlib):
  pip install numpy scipy rdkit ttkbootstrap
  + PyMOL, Meeko, AutoDock Vina installed separately (paths set in Settings)
"""

from __future__ import annotations

import json
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

try:
    import ttkbootstrap as ttkb

    _HAS_TTKB = True
except ImportError:  # graceful fallback when ttkbootstrap isn't installed
    ttkb = None  # type: ignore[assignment]
    _HAS_TTKB = False


class LogCollector:
    """Captures every log message with a timestamp for later file export."""

    def __init__(self) -> None:
        self._entries: list[tuple[str, str, str]] = []  # (timestamp, message, level)

    def __call__(self, message: str, level: str = "info") -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        self._entries.append((ts, message, level))

    def save(self, path: Path | str) -> None:
        """Write all collected entries to a plain-text log file."""
        lines: list[str] = []
        for ts, msg, _level in self.entries:
            lines.append(f"[{ts}]  {msg}")
        Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")

    @property
    def entries(self) -> list[tuple[str, str, str]]:
        return list(self._entries)

    def clear(self) -> None:
        self._entries.clear()


# ── Local pipeline modules ──────────────────────────────────────────────────
try:
    from hetatm_parser import build_chain_ligand_map, parse_hetatm_records
    from llm_export import ExportOptions, export_experiment
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
        self.root.title("MOLECULAR DOCKING PIPELINE")
        self.root.geometry("980x760")
        self.root.minsize(820, 620)

        # State Variables
        self.config: dict[str, Any] = load_config() if not _MISSING_MODULES else {}
        self.log_queue: queue.Queue[tuple[str, Any, str]] = queue.Queue()
        self.is_running: bool = False
        self._stop_requested: bool = False
        self._log_collector = LogCollector()
        self._log_buffer: list[tuple[str, str, str]] = []  # (timestamp, message, level)
        self.profile_name: str = "default"

        # Active-site selection state (populated at receptor load time)
        self.selected_chain: str = ""
        self.selected_ligand_code: str = ""
        self.selected_ligand_instance: dict = {}
        self._chain_map: dict[str, list[dict]] = {}
        self._ligand_agg: dict[str, dict] = {}
        self._lig_index_resn: dict[int, str] = {}
        self._lig_index_instance: dict[int, dict] = {}
        self._lig_greyed: dict[int, bool] = {}
        self._resn_chains: dict[str, list[str]] = {}
        self.active_site_chain_var = tk.StringVar()
        self.active_site_note_var = tk.StringVar(
            value="Load a receptor structure (PDB or mmCIF) to populate the "
            "active-site selectors."
        )

        self._setup_styles()
        self._build_ui()
        self._poll()  # Start GUI ↔ thread message pump

        # Warn if sibling modules couldn't be imported
        if _MISSING_MODULES:
            self._write_log(
                f"⚠  Could not import pipeline modules: {_MISSING_MODULES}\n"
                "   Make sure settings.py, pipeline.py, hetatm_parser.py, "
                "receptor_prep.py,\n"
                "   ligand_prep.py, and docking.py are in the same folder.",
                "warning",
            )

    # ── Styles ──────────────────────────────────────────────────────────────
    def _setup_styles(self) -> None:
        if _HAS_TTKB:
            self._style = ttkb.Style(theme="darkly")
        else:
            # Graceful fallback when ttkbootstrap isn't installed
            self._style = ttk.Style()
            try:
                self._style.theme_use("clam")
            except tk.TclError:
                pass  # Fall back to default theme on platforms lacking 'clam'

        s = self._style
        colors = getattr(s, "colors", None)
        bg = getattr(colors, "bg", "#ecf0f1")
        fg = getattr(colors, "fg", "#2c3e50")
        sub = getattr(colors, "secondary", "#7f8c8d")

        s.configure("TFrame", background=bg)
        s.configure("TLabelframe", background=bg)
        s.configure(
            "TLabelframe.Label",
            font=("Segoe UI", 9, "bold"),
            foreground=fg,
            background=bg,
        )
        s.configure("Header.TLabel", font=("Segoe UI", 13, "bold"), foreground=fg)
        s.configure("Sub.TLabel", font=("Segoe UI", 8), foreground=sub)
        s.configure("CardHead.TLabel", font=("Segoe UI", 9, "bold"), foreground=fg)
        s.configure("Card.TFrame", background=getattr(colors, "inputbg", bg))
        s.configure("Run.TButton", font=("Segoe UI", 10, "bold"), padding=7)
        s.configure("Small.TButton", font=("Segoe UI", 9), padding=4)
        s.configure("TEntry", fieldbackground=getattr(colors, "inputbg", "white"))
        s.configure("TSpinbox", fieldbackground=getattr(colors, "inputbg", "white"))
        if _HAS_TTKB:
            s.configure("success.TButton", font=("Segoe UI", 10, "bold"), padding=7)

    # ── UI Construction ─────────────────────────────────────────────────────
    def _build_ui(self) -> None:
        colors = getattr(self._style, "colors", None)
        self.root.configure(bg=getattr(colors, "bg", "#ecf0f1"))
        outer = ttk.Frame(self.root, padding="12")
        outer.pack(fill=tk.BOTH, expand=True)

        self._build_header(outer)
        self._build_inputs(outer)
        self._build_stepper(outer)

        # Pack bottom controls first using side=tk.BOTTOM to prevent clipping
        self._build_buttons(outer)
        self._build_status_bar(outer)
        self._build_log(outer)

    # ── Header row ──────────────────────────────────────────────────────────
    def _build_header(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=(0, 8))

        ttk.Label(
            row,
            text="  MOLECULAR DOCKING PIPELINE",
            style="Header.TLabel",
        ).pack(side=tk.LEFT)

        ttk.Button(
            row,
            text="⚙  Settings",
            style="Small.TButton",
            command=self._open_settings,
        ).pack(side=tk.RIGHT)

        ttk.Button(
            row,
            text="Export LLM JSON",
            style="Small.TButton",
            command=self._export_llm_json,
        ).pack(side=tk.RIGHT, padx=(0, 6))

        ttk.Label(
            row,
            text="AutoDock Vina  ·  Meeko  ·  PyMOL",
            style="Sub.TLabel",
        ).pack(side=tk.RIGHT, padx=12)

        self.profile_lbl = ttk.Label(
            row, text=f"Profile: {self.profile_name}", style="Sub.TLabel"
        )
        self.profile_lbl.pack(side=tk.RIGHT, padx=12)

    # ── Input-file section ───────────────────────────────────────────────────
    def _build_inputs(self, parent: ttk.Frame) -> None:
        card = self._card(parent, "INPUT FILES")
        card.pack(fill=tk.X, pady=(0, 8))
        frame = card.body
        frame.columnconfigure(1, weight=1)

        # Receptor row
        ttk.Label(frame, text="Receptor Structure:").grid(
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

        # Active Site Selection panel (row 1)
        active_site_panel = ttk.Frame(frame)
        active_site_panel.grid(
            row=1, column=0, columnspan=3, sticky=tk.EW, pady=(0, 8)
        )
        self._build_active_site_panel(active_site_panel)

        # Ligand rows (row 2)
        ttk.Label(frame, text="Ligand SDF(s):").grid(
            row=2, column=0, sticky=tk.NW, padx=(0, 8), pady=(4, 0)
        )

        lb_wrap = ttk.Frame(frame)
        lb_wrap.grid(row=2, column=1, sticky=tk.NSEW, padx=(0, 6), pady=(4, 0))
        frame.rowconfigure(2, weight=1)

        self.lig_tree = ttk.Treeview(
            lb_wrap,
            columns=("name", "size", "status"),
            show="headings",
            selectmode=tk.EXTENDED,
            height=4,
        )
        self.lig_tree.heading("name", text="Name")
        self.lig_tree.heading("size", text="Size")
        self.lig_tree.heading("status", text="Status")
        self.lig_tree.column("name", width=360, anchor=tk.W)
        self.lig_tree.column("size", width=80, anchor=tk.E)
        self.lig_tree.column("status", width=150, anchor=tk.W)
        self.lig_tree.tag_configure("pending", foreground="#8b949e")
        self.lig_tree.tag_configure("running", foreground="#79c0ff")
        self.lig_tree.tag_configure("ok", foreground="#56d364")
        self.lig_tree.tag_configure("err", foreground="#ff7b72")
        self.lig_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        ysb = ttk.Scrollbar(lb_wrap, orient=tk.VERTICAL, command=self.lig_tree.yview)
        ysb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lig_tree.configure(yscrollcommand=ysb.set)

        # Ligand action buttons (row 2, column 2)
        lig_btns = ttk.Frame(frame)
        lig_btns.grid(row=2, column=2, sticky=tk.N, pady=(4, 0))

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
            foreground=getattr(getattr(self._style, "colors", None), "secondary", "gray"),
            font=("Segoe UI", 8),
        ).pack()

    def _card(self, parent: tk.Widget, title: str) -> tk.Frame:
        """A dark card container: tk.Frame with a 1px highlight border + title.

        Returns a plain tk.Frame so the 1px border renders reliably across
        ttkbootstrap themes; ttk widgets are placed in ``card.body``.
        """
        colors = getattr(self._style, "colors", None)
        bg = getattr(colors, "inputbg", "#ecf0f1")
        border = getattr(colors, "border", "#454545")
        card = tk.Frame(
            parent,
            bg=bg,
            highlightbackground=border,
            highlightcolor=border,
            highlightthickness=1,
        )
        ttk.Label(card, text=title, style="CardHead.TLabel").pack(
            anchor=tk.W, padx=10, pady=(8, 0)
        )
        body = ttk.Frame(card, style="Card.TFrame")
        body.pack(fill=tk.BOTH, expand=True, padx=10, pady=(4, 10))
        card.body = body
        return card

    # ── Active Site Selection panel ──────────────────────────────────────────
    def _build_active_site_panel(self, panel: ttk.Frame) -> None:
        panel.columnconfigure(1, weight=1)

        ttk.Label(
            panel, text="Active Site Selection", style="CardHead.TLabel"
        ).grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 4))

        # Chain selector
        ttk.Label(panel, text="Chain:").grid(
            row=1, column=0, sticky=tk.W, padx=(0, 8), pady=(0, 4)
        )
        self.chain_combo = ttk.Combobox(
            panel,
            textvariable=self.active_site_chain_var,
            state="readonly",
            width=6,
        )
        self.chain_combo.grid(row=1, column=1, sticky=tk.W, pady=(0, 4))
        self.chain_combo.bind("<<ComboboxSelected>>", self._on_chain_selected)

        # Ligand selector
        ttk.Label(panel, text="Ligand:").grid(
            row=2, column=0, sticky=tk.NW, padx=(0, 8), pady=(0, 4)
        )
        lig_wrap = ttk.Frame(panel)
        lig_wrap.grid(row=2, column=1, sticky=tk.EW, pady=(0, 4))

        self.lig_sel = tk.Listbox(
            lig_wrap,
            height=3,
            selectmode=tk.SINGLE,
            exportselection=False,
            bg="#0d1117",
            fg="#c9d1d9",
            selectbackground="#264f78",
            relief="solid",
            bd=1,
            font=("Segoe UI", 9),
            activestyle="dotbox",
        )
        self.lig_sel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lig_ysb = ttk.Scrollbar(
            lig_wrap, orient=tk.VERTICAL, command=self.lig_sel.yview
        )
        lig_ysb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lig_sel.configure(yscrollcommand=lig_ysb.set)
        self.lig_sel.bind("<<ListboxSelect>>", self._on_ligand_select)

        # Multi-chain guidance
        ttk.Label(
            panel,
            textvariable=self.active_site_note_var,
            style="Sub.TLabel",
            wraplength=560,
            justify=tk.LEFT,
        ).grid(row=3, column=0, columnspan=3, sticky=tk.W, pady=(2, 0))

    def _on_chain_selected(self, *_args) -> None:
        chain = self.active_site_chain_var.get()
        self.selected_chain = chain
        self.selected_ligand_code = ""
        self.selected_ligand_instance = {}
        self._render_ligand_list(chain)
        self.active_site_note_var.set("")

    def _on_ligand_select(self, _event=None) -> None:
        sel = self.lig_sel.curselection()
        if not sel:
            self.selected_ligand_code = ""
            self.selected_ligand_instance = {}
            self.active_site_note_var.set("")
            return
        idx = sel[0]
        if self._lig_greyed.get(idx, False):
            self.lig_sel.selection_clear(0, tk.END)
            self.selected_ligand_code = ""
            self.selected_ligand_instance = {}
            return
        resn = self._lig_index_resn.get(idx, "")
        self.selected_ligand_code = resn
        self.selected_ligand_instance = dict(self._lig_index_instance.get(idx, {}))
        chains = self._resn_chains.get(resn, [])
        if len(chains) > 1:
            self.active_site_note_var.set(
                f"{resn} is present in {len(chains)} chains "
                f"({', '.join(chains)}). Select the chain and ligand you "
                f"want to centre the docking grid on."
            )
        else:
            self.active_site_note_var.set("")

    def _render_ligand_list(self, chain: str) -> None:
        """Rebuild the ligand selector for *chain*; ligands absent there are greyed."""
        self.lig_sel.delete(0, tk.END)
        self._lig_index_resn = {}
        self._lig_index_instance = {}
        self._lig_greyed = {}
        self._resn_chains = {}

        for resn in sorted(self._ligand_agg.keys()):
            entry = self._ligand_agg[resn]
            all_chains = entry["chains"]
            self._resn_chains[resn] = all_chains

            instances = entry["per_chain"].get(chain, [])
            greyed = not instances
            if not instances:
                instances = [entry["per_chain"][all_chains[0]][0]]

            for inst in instances:
                key = inst["instance_key"]
                altloc = key.get("altloc") or "-"
                label = (
                    f"{resn}  │  resi {inst['resi']}  │  alt {altloc}  │  "
                    f"{inst['heavy_atom_count']} heavy atoms  │  chains: {', '.join(all_chains)}"
                )
                idx = self.lig_sel.size()
                self.lig_sel.insert(tk.END, label)
                self._lig_index_resn[idx] = resn
                self._lig_index_instance[idx] = dict(key)
                self._lig_greyed[idx] = greyed
                if greyed:
                    self.lig_sel.itemconfig(idx, foreground="#aab2b8")

    def _aggregate_ligands(self, chain_map: dict[str, list[dict]]) -> dict[str, dict]:
        agg: dict[str, dict] = {}
        for chain, hetatms in chain_map.items():
            for h in hetatms:
                resn = h["resn"]
                entry = agg.setdefault(resn, {"chains": [], "per_chain": {}})
                entry["per_chain"].setdefault(chain, []).append(h)
                if chain not in entry["chains"]:
                    entry["chains"].append(chain)
        for entry in agg.values():
            entry["chains"].sort()
            for instances in entry["per_chain"].values():
                instances.sort(key=lambda item: (item["resi"], item["instance_id"]))
        return agg

    # ── Status bar / progress ────────────────────────────────────────────────
    def _build_status_bar(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(side=tk.BOTTOM, fill=tk.X, pady=(4, 8))

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
        card = self._card(parent, "PIPELINE LOG")
        card.pack(fill=tk.BOTH, expand=True, pady=(0, 8))
        frame = card.body

        # Filter toggle row
        filter_row = ttk.Frame(frame)
        filter_row.pack(fill=tk.X, pady=(0, 4))
        self._errors_only_var = tk.BooleanVar(value=False)
        if _HAS_TTKB:
            filter_cb = ttkb.Checkbutton(
                filter_row,
                text="Warnings & errors only",
                variable=self._errors_only_var,
                bootstyle="round-toggle",
                command=self._re_render_log,
            )
        else:
            filter_cb = ttk.Checkbutton(
                filter_row,
                text="Warnings & errors only",
                variable=self._errors_only_var,
                command=self._re_render_log,
            )
        filter_cb.pack(side=tk.LEFT)

        self.log_widget = scrolledtext.ScrolledText(
            frame,
            wrap=tk.WORD,
            bg="#0d1117",
            fg="#c9d1d9",
            font=("Consolas", 9),
            height=10,
            insertbackground="white",
            selectbackground="#264f78",
            relief="flat",
            padx=6,
            pady=6,
        )
        self.log_widget.pack(fill=tk.BOTH, expand=True)

        # Keep the log read-only while preserving standard copy shortcuts.
        self.log_widget.bind("<Control-c>", self._copy_log_selection)
        if platform.system() == "Darwin":
            self.log_widget.bind("<Command-c>", self._copy_log_selection)
        self.log_widget.bind("<Key>", lambda e: "break")

        _tags: dict[str, tuple[str, str]] = {
            "info": ("#79c0ff", "normal"),
            "success": ("#56d364", "normal"),
            "warning": ("#e3b341", "normal"),
            "error": ("#ff7b72", "normal"),
            "header": ("#d2a8ff", "bold"),
            "dim": ("#8b949e", "normal"),
            "bold": ("#f0f6fc", "bold"),
        }
        for tag, (color, weight) in _tags.items():
            self.log_widget.tag_config(
                tag, foreground=color, font=("Consolas", 9, weight)
            )

    # ── Stage stepper ───────────────────────────────────────────────────────
    # (name, progress % at which the stage is considered complete)
    _STAGES: list[tuple[str, float]] = [
        ("Setup", 3.0),
        ("Active Site", 15.0),
        ("Clean Receptor", 28.0),
        ("Grid Box", 33.0),
        ("Receptor PDBQT", 40.0),
        ("Native Ligand", 42.0),
        ("Ligand Prep", 55.0),
        ("Redock", 65.0),
        ("Docking", 95.0),
        ("Interactions", 99.0),
        ("Summary", 100.0),
    ]

    def _build_stepper(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=(0, 8))

        colors = getattr(self._style, "colors", None)
        dot_bg = getattr(colors, "bg", "#222222")
        self._steps: list[tuple[tk.Label, ttk.Label]] = []
        for name, _done_at in self._STAGES:
            step = ttk.Frame(row)
            step.pack(side=tk.LEFT, expand=True, fill=tk.X)
            dot = tk.Label(
                step,
                text="○",
                fg="#8b949e",
                bg=dot_bg,
                font=("Segoe UI", 11, "bold"),
            )
            dot.pack(side=tk.LEFT, padx=(0, 2))
            nm = ttk.Label(step, text=name, style="Sub.TLabel")
            nm.pack(side=tk.LEFT)
            self._steps.append((dot, nm))
        self._update_stage(0)

    def _update_stage(self, pct: float) -> None:
        """Recolor the stepper dots: done=success ●, current=primary ●, pending=○."""
        if not hasattr(self, "_steps"):
            return
        current = len(self._STAGES)
        for i, (_name, done_at) in enumerate(self._STAGES):
            if pct < done_at:
                current = i
                break
        colors = getattr(self._style, "colors", None)
        done = getattr(colors, "success", "#00bc8c")
        active = getattr(colors, "primary", "#375a7f")
        pending = getattr(colors, "secondary", "#444444")
        for i, (dot, _nm) in enumerate(self._steps):
            if i < current:
                dot.configure(text="●", fg=done)
            elif i == current:
                dot.configure(text="●", fg=active)
            else:
                dot.configure(text="○", fg=pending)

    # ── Bottom button row ────────────────────────────────────────────────────
    def _build_buttons(self, parent: ttk.Frame) -> None:
        row = ttk.Frame(parent)
        row.pack(side=tk.BOTTOM, fill=tk.X, pady=(4, 0))

        if _HAS_TTKB:
            self.run_btn = ttkb.Button(
                row,
                text="▶   Run Docking Pipeline",
                bootstyle="success",
                command=self._start_pipeline,
            )
        else:
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
            title="Select Receptor Structure File",
            parent=self.root,
            filetypes=[
                ("Structure files", "*.pdb *.PDB *.cif *.CIF *.mmcif"),
                ("All files", "*"),
            ],
        )
        if path:
            self.receptor_var.set(path)
            self._write_log(f"Receptor selected: {Path(path).name}", "dim")
            self._load_receptor_ligands(path)

    def _load_receptor_ligands(self, path: str) -> None:
        self.selected_chain = ""
        self.selected_ligand_code = ""
        self.selected_ligand_instance = {}
        self.active_site_chain_var.set("")
        self.active_site_note_var.set("Parsing receptor HETATM records…")
        self.lig_sel.delete(0, tk.END)
        self._lig_index_resn = {}
        self._lig_index_instance = {}
        self._lig_greyed = {}
        self._resn_chains = {}

        try:
            hetatms = parse_hetatm_records(path)
            chain_map = build_chain_ligand_map(hetatms)
        except Exception as exc:
            self.active_site_note_var.set(
                "Failed to parse receptor ligands — see the log for details."
            )
            self._write_log(f"  ✘  Receptor ligand parsing failed: {exc}", "error")
            self.chain_combo["values"] = []
            return

        self._chain_map = chain_map
        self._ligand_agg = self._aggregate_ligands(chain_map)

        chains = sorted(chain_map.keys())
        self.chain_combo["values"] = chains

        if not chains:
            self.active_site_note_var.set(
                "No candidate ligands found in this receptor "
                "(only water/ions/buffers were present)."
            )
            self._write_log(
                "  ⚠  No candidate ligands after filtering non-ligands — "
                "select a different receptor or re-check the structure file.",
                "warning",
            )
            return

        self._write_log(
            f"  Active-site panel: {len(chains)} chain(s), "
            f"{len(self._ligand_agg)} candidate ligand type(s).",
            "info",
        )

        if len(chains) == 1:
            self.chain_combo.set(chains[0])
            self.selected_chain = chains[0]
            self._render_ligand_list(chains[0])
        else:
            self.active_site_note_var.set(
                f"Candidate ligands found in chains: {', '.join(chains)}. "
                "Select a chain to enable its ligands."
            )

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
        added = 0
        for p in paths:
            if self.lig_tree.exists(p):
                continue
            size = os.path.getsize(p)
            self.lig_tree.insert(
                "",
                tk.END,
                iid=p,
                values=(Path(p).name, f"{size / 1024:.1f} KB", "Pending"),
                tags=("pending",),
            )
            added += 1
        if added:
            self._refresh_lig_count()
            self._write_log(
                f"Added {added} ligand(s). Total: {len(self.lig_tree.get_children())}",
                "dim",
            )

    def _remove_ligands(self) -> None:
        for iid in self.lig_tree.selection():
            self.lig_tree.delete(iid)
        self._refresh_lig_count()

    def _clear_ligands(self) -> None:
        self.lig_tree.delete(*self.lig_tree.get_children())
        self._refresh_lig_count()

    def _refresh_lig_count(self) -> None:
        n = len(self.lig_tree.get_children())
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
        dlg = SettingsDialog(
            self.root,
            self.config,
            style=self._style,
            profile=getattr(self, "profile_name", "default"),
        )
        self.profile_name = getattr(dlg, "saved_profile", "default")
        self.config = load_config(self.profile_name)
        self.profile_lbl.configure(text=f"Profile: {self.profile_name}")
        self._write_log("Settings updated.", "dim")

    def _export_llm_json(self) -> None:
        """Export a compact self-contained docking summary from a run manifest."""
        manifest = filedialog.askopenfilename(
            title="Select DockLLM Source Manifest",
            parent=self.root,
            filetypes=[("DockLLM source", "dockllm-source.json"), ("JSON files", "*.json")],
        )
        if not manifest:
            return
        output = filedialog.asksaveasfilename(
            title="Export LLM Analysis File",
            parent=self.root,
            defaultextension=".json",
            initialfile="docking_experiment.llm.json",
            filetypes=[("LLM analysis JSON", "*.json")],
        )
        if not output:
            return
        try:
            context = json.loads(Path(manifest).read_text(encoding="utf-8"))
            export_experiment(context, output, ExportOptions())
        except Exception as exc:
            messagebox.showerror("Export Failed", str(exc), parent=self.root)
            return
        messagebox.showinfo(
            "Export Complete",
            "Compact self-contained docking summary exported for language-model analysis.",
            parent=self.root,
        )

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
        ligands = list(self.lig_tree.get_children())

        if not receptor:
            messagebox.showerror("Missing Input", "Please select a Receptor structure file.")
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

        if not self.selected_chain or not self.selected_ligand_code or not self.selected_ligand_instance:
            messagebox.showerror(
                "No Active Site Selected",
                "No active site selected. Load a receptor structure (PDB or "
                "mmCIF) and select a chain and ligand before running.",
            )
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
        self._log_collector.clear()
        self.progress_var.set(0)
        self.pct_var.set("0%")
        self.status_var.set("Starting pipeline…")
        self.is_running = True
        self._stop_requested = False
        self.run_btn.configure(state="disabled")
        self.stop_btn.configure(state="normal")

        # Banner
        self._write_log("═" * 60, "header")
        self._write_log("  MOLECULAR DOCKING PIPELINE", "header")
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
        self._write_log(
            f"Active site : [{self.selected_ligand_code}] in chain "
            f"{self.selected_chain}",
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
            def _log_through_collector(msg: str, lvl: str = "info") -> None:
                self._log_collector(msg, lvl)
                self._queue("log", msg, lvl)

            run_pipeline(
                receptor_file=receptor,
                ligand_files=ligands,
                config=self.config,
                log_cb=_log_through_collector,
                progress_cb=self._queue_progress,
                status_cb=lambda s: self._queue("status", s),
                stop_flag=lambda: self._stop_requested,
                write_log_fn=lambda path: self._log_collector.save(path),
                selected_chain=self.selected_chain,
                selected_ligand_code=self.selected_ligand_code,
                selected_ligand_instance=self.selected_ligand_instance,
                ligand_status_cb=lambda key, status, detail: self._queue(
                    "lig_status", (key, status, detail)
                ),
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
                    self._update_stage(val)

                elif kind == "status":
                    self.status_var.set(str(value))

                elif kind == "lig_status":
                    key, status, detail = value
                    self._set_lig_status(key, status, detail)

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

    _LOG_ICONS = {"success": "✓ ", "warning": "⚠ ", "error": "✗ "}
    _LOG_ICON_CHARS = ("✓", "✗", "⚠", "✔", "✘", "❌", "⏹")

    def _copy_log_selection(self, _event=None) -> str:
        """Copy selected log text while keeping the widget read-only."""
        self.log_widget.event_generate("<<Copy>>")
        return "break"

    def _write_log(self, message: str, level: str = "info") -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        self._log_buffer.append((ts, message, level))
        errors_only = bool(
            getattr(self, "_errors_only_var", None) and self._errors_only_var.get()
        )
        if errors_only and level not in ("warning", "error"):
            return
        self._insert_log(ts, message, level)
        self.log_widget.see(tk.END)

    def _insert_log(self, ts: str, message: str, level: str) -> None:
        msg = message.lstrip()
        icon = self._LOG_ICONS.get(level, "")
        if icon and not msg.startswith(self._LOG_ICON_CHARS):
            message = icon + message
        self.log_widget.insert(tk.END, f"[{ts}]  {message}\n", level)

    def _re_render_log(self) -> None:
        """Re-render the log widget from the buffer (used by the filter toggle)."""
        self.log_widget.delete("1.0", tk.END)
        errors_only = self._errors_only_var.get()
        for ts, msg, level in self._log_buffer:
            if errors_only and level not in ("warning", "error"):
                continue
            self._insert_log(ts, msg, level)
        self.log_widget.see(tk.END)

    def _clear_log(self) -> None:
        self.log_widget.delete("1.0", tk.END)
        self._log_buffer.clear()

    _LIG_STATUS_TEXT = {
        "preparing": "Preparing…",
        "prep_failed": "Prep failed",
        "docking": "Docking…",
        "ok": "Done",
        "failed": "Failed",
    }

    def _set_lig_status(self, key: str, status: str, detail: str = "") -> None:
        if not hasattr(self, "lig_tree") or not self.lig_tree.exists(key):
            return
        text = self._LIG_STATUS_TEXT.get(status, status)
        if status == "ok" and detail:
            text = f"{text} · {detail} kcal/mol"
        tag = {
            "preparing": "running",
            "docking": "running",
            "ok": "ok",
            "prep_failed": "err",
            "failed": "err",
        }.get(status, "pending")
        row = self.lig_tree.item(key)
        values = list(row["values"])
        if len(values) >= 3:
            values[2] = text
        self.lig_tree.item(key, values=values, tags=(tag,))

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
