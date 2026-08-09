"""
Molecular Docking Pipeline
=====================================
settings.py — Settings dialog + config persistence

Handles:
  - Loading / saving config.json
  - Validating all required paths before a run
  - Tabbed Tkinter settings dialog with Browse buttons for executables & parameters
"""

import json
import os
import shutil
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, simpledialog, ttk

try:
    import ttkbootstrap as ttkb

    _HAS_TTKB = True
except ImportError:  # graceful fallback when ttkbootstrap isn't installed
    ttkb = None  # type: ignore[assignment]
    _HAS_TTKB = False

# ════════════════════════════════════════════════════════════════════════════
#  Config file location (same folder as this script)
# ════════════════════════════════════════════════════════════════════════════

CONFIG_FILE = Path(__file__).resolve().parent / "config.json"

DEFAULT_CONFIG: dict = {
    # ── Executables ──────────────────────────────────────────────────────
    "vina_exe": "",  # AutoDock Vina binary
    "pymol_exe": "",  # PyMOL binary (used with -c flag)
    "mk_prepare_ligand_cmd": "",  # mk_prepare_ligand script/executable
    "mk_prepare_receptor_cmd": "",  # mk_prepare_receptor script/executable
    "plip_cmd": shutil.which("plip") or "",  # PLIP executable or plipcmd.py
    # ── Output ───────────────────────────────────────────────────────────
    "output_dir": str(Path.home() / "docking_results"),
    # ── Docking parameters ───────────────────────────────────────────────
    "grid_padding": 4.0,  # Docking box padding in Angstroms
    "exhaustiveness": 8,
    "num_modes": 9,
    "cpu": 0,  # 0 = Vina auto-detects
}


# ════════════════════════════════════════════════════════════════════════════
#  Public helpers (used by main.py and pipeline.py)
# ════════════════════════════════════════════════════════════════════════════


def _profile_path(name: str = "default") -> Path:
    """Path of a named config profile. "default" → config.json, else config_<name>.json."""
    if not name or name == "default":
        return CONFIG_FILE
    return CONFIG_FILE.with_name(f"config_{name}.json")


def list_profiles() -> list[str]:
    """Names of all config profiles present next to config.json ("default" first)."""
    profiles = ["default"]
    if CONFIG_FILE.parent.exists():
        for p in sorted(CONFIG_FILE.parent.glob("config_*.json")):
            profiles.append(p.name[len("config_") : -len(".json")])
    return profiles


def load_config(name: str = "default") -> dict:
    """Load a named profile and merge with DEFAULT_CONFIG so new keys are always present.

    Creates the profile file with defaults if it does not yet exist.
    """
    merged = DEFAULT_CONFIG.copy()
    path = _profile_path(name)
    if path.exists():
        try:
            with open(path, "r", encoding="utf-8") as f:
                stored = json.load(f)
            merged.update(stored)
        except (json.JSONDecodeError, OSError):
            pass
    else:
        _write_json(merged, path)
    return merged


def save_config(config: dict, name: str = "default") -> None:
    """Persist config dict to the given profile file."""
    _write_json(config, _profile_path(name))


def validate_config(config: dict) -> list[str]:
    """Return a list of error strings for missing / invalid settings.

    An empty list means the config is ready to run.
    """
    errors: list[str] = []

    # Required executables
    exe_fields = {
        "vina_exe": "AutoDock Vina executable",
        "pymol_exe": "PyMOL executable",
        "mk_prepare_ligand_cmd": "mk_prepare_ligand command/script",
        "mk_prepare_receptor_cmd": "mk_prepare_receptor command/script",
        "plip_cmd": "PLIP command/script",
    }
    for key, label in exe_fields.items():
        val = str(config.get(key, "")).strip()
        if not val:
            errors.append(f"{label} path is not set.")
        elif not Path(val).expanduser().exists() and shutil.which(val) is None:
            errors.append(f"{label} not found at: {val}")

    # API key is optional — active-site selection is now done manually in the
    # GUI and requires no network calls.

    # Output directory
    out = str(config.get("output_dir", "")).strip()
    if not out:
        errors.append("Output directory is not set.")
    else:
        out_path = Path(out)
        if out_path.exists() and not out_path.is_dir():
            errors.append(f"Output path exists but is not a directory: {out}")

    return errors


# ════════════════════════════════════════════════════════════════════════════
#  Settings Dialog
# ════════════════════════════════════════════════════════════════════════════


class SettingsDialog:
    """Modal Toplevel settings GUI dialog."""

    def __init__(
        self,
        parent: tk.Tk,
        config: dict,
        style=None,
        profile: str = "default",
    ):
        self._config = config.copy()
        self._vars: dict = {}
        self._style = style
        self._profile = profile
        self.saved_profile = profile

        dlg = tk.Toplevel(parent)
        dlg.title("Settings — Molecular Docking Pipeline")
        dlg.geometry("740x660")
        dlg.resizable(True, False)
        dlg.transient(parent)
        dlg.grab_set()
        self._dlg = dlg

        self._setup_styles()
        self._build(dlg)
        dlg.wait_window()

    def _setup_styles(self):
        if self._style is None:
            self._style = ttkb.Style(theme="darkly") if _HAS_TTKB else ttk.Style()
        s = self._style
        colors = getattr(s, "colors", None)
        bg = getattr(colors, "bg", "#f5f6fa")
        fg = getattr(colors, "fg", "#000000")
        hint = getattr(colors, "secondary", "#7f8c8d")
        s.configure("Dialog.TFrame", background=bg)
        s.configure("Dialog.TLabelframe", background=bg)
        s.configure(
            "Dialog.TLabelframe.Label",
            font=("Segoe UI", 9, "bold"),
            background=bg,
            foreground=fg,
        )
        s.configure(
            "Hint.TLabel",
            foreground=hint,
            font=("Segoe UI", 8),
            background=bg,
        )
        s.configure(
            "SectionHead.TLabel",
            font=("Segoe UI", 9, "bold"),
            background=bg,
            foreground=fg,
        )

    def _build(self, dlg: tk.Toplevel):
        colors = getattr(self._style, "colors", None)
        dlg.configure(bg=getattr(colors, "bg", "#f5f6fa"))
        outer = ttk.Frame(dlg, padding="14", style="Dialog.TFrame")
        outer.pack(fill=tk.BOTH, expand=True)

        self._build_profile_row(outer)

        nb = ttk.Notebook(outer)
        nb.pack(fill=tk.BOTH, expand=True, pady=(0, 12))

        self._tab_executables(nb)
        self._tab_docking(nb)

        sep = ttk.Separator(outer, orient=tk.HORIZONTAL)
        sep.pack(fill=tk.X, pady=(0, 8))

        btn_row = ttk.Frame(outer, style="Dialog.TFrame")
        btn_row.pack(fill=tk.X)

        ttk.Button(
            btn_row, text="Validate Paths", command=self._validate_and_report
        ).pack(side=tk.LEFT)
        ttk.Button(
            btn_row, text="Cancel", command=self._dlg.destroy
        ).pack(side=tk.RIGHT, padx=(6, 0))
        ttk.Button(btn_row, text="Save", command=self._save).pack(side=tk.RIGHT)

    # ── Tab 1 — Executables & Paths ─────────────────────────────────────────

    def _tab_executables(self, nb: ttk.Notebook):
        frame = ttk.Frame(nb, padding="12", style="Dialog.TFrame")
        nb.add(frame, text=" Executables & Paths ")
        frame.columnconfigure(1, weight=1)

        rows = [
            (
                "AutoDock Vina:",
                "vina_exe",
                "The vina binary. e.g. /usr/local/bin/vina or C:\\Program Files\\AutoDock Vina\\vina.exe",
                "file",
            ),
            (
                "PyMOL:",
                "pymol_exe",
                "PyMOL binary (headless mode). e.g. /usr/bin/pymol or C:\\Program Files\\PyMOL\\PyMOLWin.exe",
                "file",
            ),
            (
                "mk_prepare_ligand:",
                "mk_prepare_ligand_cmd",
                "Meeko ligand prep script or executable. e.g. mk_prepare_ligand (on PATH) or mk_prepare_ligand.py",
                "file",
            ),
            (
                "mk_prepare_receptor:",
                "mk_prepare_receptor_cmd",
                "Meeko receptor prep script or executable. e.g. mk_prepare_receptor (on PATH)",
                "file",
            ),
            (
                "PLIP:",
                "plip_cmd",
                "PLIP executable or plipcmd.py. Used for production-pose interaction reports.",
                "file",
            ),
            (
                "Output directory:",
                "output_dir",
                "Root folder for results. A sub-folder is created per run.",
                "dir",
            ),
        ]

        for i, (label, key, hint, kind) in enumerate(rows):
            base_row = i * 3

            ttk.Label(frame, text=label, style="SectionHead.TLabel").grid(
                row=base_row, column=0, sticky=tk.W, padx=(0, 10), pady=(10, 0)
            )

            var = tk.StringVar(value=str(self._config.get(key, "")))
            self._vars[key] = var

            entry = ttk.Entry(frame, textvariable=var)
            entry.grid(
                row=base_row,
                column=1,
                sticky=tk.EW,
                padx=(0, 6),
                pady=(10, 0),
            )

            browse_cmd = (
                (lambda v=var: self._browse_file(v))
                if kind == "file"
                else (lambda v=var: self._browse_dir(v))
            )
            ttk.Button(frame, text="Browse…", width=9, command=browse_cmd).grid(
                row=base_row, column=2, pady=(10, 0)
            )

            ttk.Label(frame, text=hint, style="Hint.TLabel").grid(
                row=base_row + 1,
                column=0,
                columnspan=3,
                sticky=tk.W,
                padx=(0, 10),
                pady=(1, 0),
            )

            ttk.Frame(frame, height=2, style="Dialog.TFrame").grid(
                row=base_row + 2, column=0
            )

    # ── Tab 3 — Docking Parameters ──────────────────────────────────────────

    def _tab_docking(self, nb: ttk.Notebook):
        frame = ttk.Frame(nb, padding="12", style="Dialog.TFrame")
        nb.add(frame, text=" Docking Parameters ")
        frame.columnconfigure(1, weight=0)
        frame.columnconfigure(2, weight=1)

        # Updated parameters list including grid_padding (float)
        params = [
            (
                "Grid Padding (Å):",
                "grid_padding",
                float(self._config.get("grid_padding", 4.0)),
                3.0,
                12.0,
                0.5,
                "float",
                "Distance added to each side of the native ligand bounding box (default 4.0 Å).",
            ),
            (
                "Exhaustiveness:",
                "exhaustiveness",
                int(self._config.get("exhaustiveness", 8)),
                1,
                64,
                1,
                "int",
                "Search thoroughness (default 8). Higher = better pose search but slower.",
            ),
            (
                "Number of poses:",
                "num_modes",
                int(self._config.get("num_modes", 9)),
                1,
                20,
                1,
                "int",
                "Maximum binding modes (poses) Vina will generate per ligand.",
            ),
            (
                "CPU cores (0 = auto):",
                "cpu",
                int(self._config.get("cpu", 0)),
                0,
                128,
                1,
                "int",
                "Cores assigned to each Vina run. 0 lets Vina auto-detect available cores.",
            ),
        ]

        for i, (label, key, val, lo, hi, step, dtype, hint) in enumerate(params):
            base = i * 3

            ttk.Label(frame, text=label, style="SectionHead.TLabel").grid(
                row=base, column=0, sticky=tk.W, padx=(0, 10), pady=(10, 0)
            )

            if dtype == "float":
                var = tk.DoubleVar(value=float(val))
            else:
                var = tk.IntVar(value=int(val))

            self._vars[key] = var

            ttk.Spinbox(
                frame,
                textvariable=var,
                from_=lo,
                to=hi,
                increment=step,
                width=8,
            ).grid(
                row=base, column=1, sticky=tk.W, padx=(0, 10), pady=(10, 0)
            )

            ttk.Label(
                frame, text=hint, style="Hint.TLabel", wraplength=420
            ).grid(
                row=base + 1,
                column=0,
                columnspan=3,
                sticky=tk.W,
                padx=(0, 10),
                pady=(2, 0),
            )

            ttk.Frame(frame, height=2, style="Dialog.TFrame").grid(
                row=base + 2, column=0
            )

        row_offset = len(params) * 3
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(
            row=row_offset,
            column=0,
            columnspan=3,
            sticky=tk.EW,
            pady=(12, 8),
        )

        notes = (
            "ℹ  Redocking RMSD validation threshold: 2.0 Å (< 2.0 Å = PASS).\n"
            "ℹ  Dynamic grid padding allows customisation based on ligand size & active site flexibility."
        )
        ttk.Label(
            frame, text=notes, style="Hint.TLabel", justify=tk.LEFT
        ).grid(row=row_offset + 1, column=0, columnspan=3, sticky=tk.W)

    # ── Profile helpers ─────────────────────────────────────────────────────

    def _build_profile_row(self, outer: ttk.Frame) -> None:
        row = ttk.Frame(outer, style="Dialog.TFrame")
        row.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(row, text="Profile:", style="SectionHead.TLabel").pack(
            side=tk.LEFT, padx=(0, 6)
        )
        self._profile_var = tk.StringVar(value=self._profile)
        self._profile_combo = ttk.Combobox(
            row,
            textvariable=self._profile_var,
            values=list_profiles(),
            state="readonly",
            width=26,
        )
        self._profile_combo.pack(side=tk.LEFT)
        self._profile_combo.bind("<<ComboboxSelected>>", self._on_profile_selected)

        ttk.Button(
            row, text="New…", width=8, command=self._new_profile
        ).pack(side=tk.LEFT, padx=(6, 0))

        ttk.Label(
            row,
            text="Save writes to the selected profile file.",
            style="Hint.TLabel",
        ).pack(side=tk.LEFT, padx=(10, 0))

    def _on_profile_selected(self, _event=None) -> None:
        name = self._profile_var.get()
        if not name or name == self._profile:
            return
        self._profile = name
        self._config = load_config(name)
        self._populate()

    def _new_profile(self) -> None:
        name = simpledialog.askstring(
            "New Profile", "Profile name:", parent=self._dlg
        )
        if not name:
            return
        name = "".join(c for c in name.strip() if c.isalnum() or c in "-_ ").strip()
        if not name:
            return
        if name in list_profiles():
            messagebox.showwarning(
                "Profile Exists",
                f"Profile '{name}' already exists.",
                parent=self._dlg,
            )
            return
        self._profile = name
        self._config = load_config(name)  # creates config_<name>.json with defaults
        self._profile_var.set(name)
        self._profile_combo["values"] = list_profiles()
        self._populate()

    def _populate(self) -> None:
        """Refresh all bound vars from self._config (used when switching profiles)."""
        for key, var in self._vars.items():
            val = self._config.get(key, DEFAULT_CONFIG.get(key, ""))
            try:
                var.set(val)
            except tk.TclError:
                pass

    # ── Internal Helpers ────────────────────────────────────────────────────

    def _browse_file(self, var: tk.StringVar):
        initial = _initial_dir(var.get())
        path = filedialog.askopenfilename(
            title="Select File",
            parent=self._dlg,
            initialdir=initial,
            filetypes=[
                ("All files", "*"),
                ("Executable / Script", "*.exe *.py *.bat *.cmd *.sh"),
            ],
        )
        if path:
            var.set(path)

    def _browse_dir(self, var: tk.StringVar):
        initial = _initial_dir(var.get())
        path = filedialog.askdirectory(
            title="Select Directory", parent=self._dlg, initialdir=initial
        )
        if path:
            var.set(path)

    def _validate_and_report(self):
        tmp = self._collect_values()
        errors = validate_config(tmp)
        if errors:
            messagebox.showwarning(
                "Validation Issues",
                "The following problems were found:\n\n"
                + "\n".join(f"  •  {e}" for e in errors),
                parent=self._dlg,
            )
        else:
            messagebox.showinfo(
                "Validation Passed",
                "All paths exist and required fields are filled in.\n"
                "You can now run the pipeline.",
                parent=self._dlg,
            )

    def _collect_values(self) -> dict:
        result = self._config.copy()
        for key, var in self._vars.items():
            try:
                result[key] = var.get()
            except tk.TclError:
                result[key] = DEFAULT_CONFIG.get(key, "")

        return result

    def _save(self):
        new_cfg = self._collect_values()

        # Coerce type fields safely
        for key in ("exhaustiveness", "num_modes", "cpu"):
            try:
                new_cfg[key] = int(new_cfg[key])
            except (ValueError, TypeError):
                new_cfg[key] = DEFAULT_CONFIG[key]

        try:
            new_cfg["grid_padding"] = float(new_cfg.get("grid_padding", 4.0))
        except (ValueError, TypeError):
            new_cfg["grid_padding"] = DEFAULT_CONFIG["grid_padding"]

        save_config(new_cfg, self._profile)
        self._config.update(new_cfg)
        self.saved_profile = self._profile
        self._dlg.destroy()


# ════════════════════════════════════════════════════════════════════════════
#  Utility functions
# ════════════════════════════════════════════════════════════════════════════


def _write_json(data: dict, path: Path | None = None) -> None:
    """Atomically write a config file."""
    path = path or CONFIG_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        tmp.replace(path)
    except OSError:
        tmp.unlink(missing_ok=True)
        raise


def _initial_dir(current_val: str) -> str:
    if current_val:
        p = Path(current_val)
        candidate = p.parent if p.is_file() else p
        if candidate.exists():
            return str(candidate)
    return str(Path.home())

