"""
AI-Guided Molecular Docking Pipeline
=====================================
settings.py — Settings dialog + config persistence

Handles:
  - Loading / saving config.json
  - Validating all required paths before a run
  - Tabbed Tkinter settings dialog with Browse buttons for executables & parameters
"""

import json
import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

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
    # ── AI / API ─────────────────────────────────────────────────────────
    "api_key": "",
    "api_provider": "openrouter",  # openrouter | openai | nvidia | custom
    "api_base_url": "https://openrouter.ai/api/v1",
    "model": "meta-llama/llama-3.1-70b-instruct",
    # ── Output ───────────────────────────────────────────────────────────
    "output_dir": str(Path.home() / "docking_results"),
    # ── Docking parameters ───────────────────────────────────────────────
    "grid_padding": 6.0,  # Docking box padding in Angstroms
    "exhaustiveness": 8,
    "num_modes": 9,
    "cpu": 0,  # 0 = Vina auto-detects
}


# ════════════════════════════════════════════════════════════════════════════
#  Public helpers (used by main.py and pipeline.py)
# ════════════════════════════════════════════════════════════════════════════


def load_config() -> dict:
    """Load config.json and merge with DEFAULT_CONFIG so new keys are always present.

    Creates config.json with defaults if it does not yet exist.
    """
    merged = DEFAULT_CONFIG.copy()
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, "r", encoding="utf-8") as f:
                stored = json.load(f)
            merged.update(stored)
        except (json.JSONDecodeError, OSError):
            pass
    else:
        _write_json(merged)
    return merged


def save_config(config: dict) -> None:
    """Persist config dict to config.json."""
    _write_json(config)


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
    }
    for key, label in exe_fields.items():
        val = str(config.get(key, "")).strip()
        if not val:
            errors.append(f"{label} path is not set.")
        elif not Path(val).exists():
            errors.append(f"{label} not found at: {val}")

    # API key
    if not str(config.get("api_key", "")).strip():
        errors.append("API key is not set.")

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

    def __init__(self, parent: tk.Tk, config: dict):
        self._config = config.copy()
        self._vars: dict = {}
        self._api_entry: ttk.Entry | None = None
        self._model_combo: ttk.Combobox | None = None
        self._or_status_var: tk.StringVar | None = None

        dlg = tk.Toplevel(parent)
        dlg.title("Settings — AI GUIDED DOCKING PROTOCOL (AGDP)")
        dlg.geometry("740x580")
        dlg.resizable(True, False)
        dlg.transient(parent)
        dlg.grab_set()
        self._dlg = dlg

        self._setup_styles()
        self._build(dlg)
        dlg.wait_window()

    def _setup_styles(self):
        s = ttk.Style()
        bg = "#f5f6fa"
        s.configure("Dialog.TFrame", background=bg)
        s.configure("Dialog.TLabelframe", background=bg)
        s.configure(
            "Dialog.TLabelframe.Label",
            font=("Segoe UI", 9, "bold"),
            background=bg,
        )
        s.configure(
            "Hint.TLabel",
            foreground="#7f8c8d",
            font=("Segoe UI", 8),
            background=bg,
        )
        s.configure(
            "Link.TLabel",
            foreground="#2980b9",
            font=("Segoe UI", 8),
            background=bg,
            cursor="hand2",
        )
        s.configure(
            "SectionHead.TLabel", font=("Segoe UI", 9, "bold"), background=bg
        )

    def _build(self, dlg: tk.Toplevel):
        dlg.configure(bg="#f5f6fa")
        outer = ttk.Frame(dlg, padding="14", style="Dialog.TFrame")
        outer.pack(fill=tk.BOTH, expand=True)

        nb = ttk.Notebook(outer)
        nb.pack(fill=tk.BOTH, expand=True, pady=(0, 12))

        self._tab_executables(nb)
        self._tab_api(nb)
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
                "The vina binary. e.g. /usr/local/bin/vina",
                "file",
            ),
            (
                "PyMOL:",
                "pymol_exe",
                "PyMOL binary called with -c flag. e.g. /usr/bin/pymol",
                "file",
            ),
            (
                "mk_prepare_ligand:",
                "mk_prepare_ligand_cmd",
                "Meeko ligand prep script or executable.",
                "file",
            ),
            (
                "mk_prepare_receptor:",
                "mk_prepare_receptor_cmd",
                "Meeko receptor prep script or executable.",
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

    # ── Tab 2 — AI / API ────────────────────────────────────────────────────

    def _tab_api(self, nb: ttk.Notebook):
        frame = ttk.Frame(nb, padding="12", style="Dialog.TFrame")
        nb.add(frame, text=" AI / API ")
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Provider:", style="SectionHead.TLabel").grid(
            row=0, column=0, sticky=tk.W, padx=(0, 10), pady=(8, 0)
        )

        PROVIDERS = {
            "OpenRouter": (
                "https://openrouter.ai/api/v1",
                "meta-llama/llama-3.1-70b-instruct",
            ),
            "OpenAI": ("https://api.openai.com/v1", "gpt-4o-mini"),
            "NVIDIA AI": (
                "https://integrate.api.nvidia.com/v1",
                "meta/llama-3.1-70b-instruct",
            ),
            "Custom": ("", ""),
        }

        stored_provider = (
            str(self._config.get("api_provider", "openrouter")).lower().strip()
        )
        label_map = {
            "openrouter": "OpenRouter",
            "openai": "OpenAI",
            "nvidia": "NVIDIA AI",
            "custom": "Custom",
        }

        provider_var = tk.StringVar(
            value=label_map.get(stored_provider, "OpenRouter")
        )
        self._vars["api_provider"] = provider_var

        provider_combo = ttk.Combobox(
            frame,
            textvariable=provider_var,
            state="readonly",
            values=list(PROVIDERS.keys()),
            width=20,
        )
        provider_combo.grid(row=0, column=1, sticky=tk.W, pady=(8, 0))

        # API Key
        ttk.Label(frame, text="API Key:", style="SectionHead.TLabel").grid(
            row=1, column=0, sticky=tk.W, padx=(0, 10), pady=(8, 0)
        )

        api_var = tk.StringVar(value=str(self._config.get("api_key", "")))
        self._vars["api_key"] = api_var

        self._api_entry = ttk.Entry(frame, textvariable=api_var, show="*")
        self._api_entry.grid(
            row=1, column=1, sticky=tk.EW, padx=(0, 6), pady=(8, 0)
        )

        ttk.Button(
            frame, text="👁", width=3, command=self._toggle_api_visibility
        ).grid(row=1, column=2, pady=(8, 0))

        # Base URL
        ttk.Label(frame, text="Base URL:", style="SectionHead.TLabel").grid(
            row=2, column=0, sticky=tk.W, padx=(0, 10), pady=(6, 0)
        )

        base_url_var = tk.StringVar(
            value=str(
                self._config.get(
                    "api_base_url", "https://openrouter.ai/api/v1"
                )
            )
        )
        self._vars["api_base_url"] = base_url_var

        ttk.Entry(frame, textvariable=base_url_var).grid(
            row=2, column=1, columnspan=2, sticky=tk.EW, pady=(6, 0)
        )

        # Model
        ttk.Label(frame, text="Model:", style="SectionHead.TLabel").grid(
            row=3, column=0, sticky=tk.W, padx=(0, 10), pady=(6, 0)
        )

        model_var = tk.StringVar(value=str(self._config.get("model", "")))
        self._vars["model"] = model_var

        self._model_combo = ttk.Combobox(
            frame, textvariable=model_var, width=48
        )
        self._model_combo.grid(row=3, column=1, sticky=tk.EW, pady=(6, 0))

        ttk.Button(
            frame,
            text="⟳ Fetch",
            width=8,
            command=lambda: self._fetch_or_models(
                api_var, base_url_var, model_var
            ),
        ).grid(row=3, column=2, pady=(6, 0))

        self._or_status_var = tk.StringVar(value="")
        ttk.Label(
            frame, textvariable=self._or_status_var, style="Hint.TLabel"
        ).grid(row=4, column=0, columnspan=3, sticky=tk.W, pady=(2, 0))

        hints = {
            "OpenRouter": (
                "→ openrouter.ai — free and paid model access",
                "https://openrouter.ai/keys",
            ),
            "OpenAI": (
                "→ platform.openai.com/api-keys",
                "https://platform.openai.com/api-keys",
            ),
            "NVIDIA AI": (
                "→ build.nvidia.com/explore/discover",
                "https://build.nvidia.com/explore/discover",
            ),
            "Custom": ("Enter any OpenAI-compatible base URL above.", ""),
        }

        link_text_var = tk.StringVar()
        link_label = ttk.Label(
            frame, textvariable=link_text_var, style="Link.TLabel"
        )
        link_label.grid(
            row=5, column=0, columnspan=3, sticky=tk.W, pady=(6, 0)
        )
        self._link_url = ""

        def _on_provider_change(*_):
            label = provider_var.get()
            info = PROVIDERS.get(label, ("", ""))
            base_url_var.set(info[0])
            if not model_var.get():
                model_var.set(info[1])
            hint, url = hints.get(label, ("", ""))
            link_text_var.set(hint)
            self._link_url = url

        provider_combo.bind("<<ComboboxSelected>>", _on_provider_change)
        link_label.bind("<Button-1>", lambda _: _open_url(self._link_url))

        _on_provider_change()
        model_var.set(str(self._config.get("model", model_var.get())))

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
                float(self._config.get("grid_padding", 6.0)),
                3.0,
                12.0,
                0.5,
                "float",
                "Distance added to each side of the native ligand bounding box (default 6.0 Å).",
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

    # ── Internal Helpers ────────────────────────────────────────────────────

    def _browse_file(self, var: tk.StringVar):
        initial = _initial_dir(var.get())
        path = filedialog.askopenfilename(
            title="Select File",
            parent=self._dlg,
            initialdir=initial,
            filetypes=[
                ("All files", "*"),
                ("Executable / Script", "*.py *.sh *.bat"),
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

    def _toggle_api_visibility(self):
        if self._api_entry is None:
            return
        currently_hidden = self._api_entry.cget("show") == "*"
        self._api_entry.configure(show="" if currently_hidden else "*")

    def _fetch_or_models(self, api_var, base_url_var, model_var):
        import json as _json
        import threading
        import urllib.request

        base = base_url_var.get().rstrip("/")
        key = api_var.get().strip()
        if not base:
            self._or_status_var.set("No base URL set.")
            return
        self._or_status_var.set("Fetching models…")

        def _fetch():
            try:
                req = urllib.request.Request(
                    f"{base}/models",
                    headers={
                        "Authorization": f"Bearer {key}",
                        "Content-Type": "application/json",
                    },
                )
                with urllib.request.urlopen(req, timeout=10) as r:
                    data = _json.loads(r.read())
                ids = sorted(
                    m.get("id", "") for m in data.get("data", []) if m.get("id")
                )

                def _apply():
                    if ids:
                        self._model_combo["values"] = ids
                        if model_var.get() not in ids:
                            model_var.set(ids[0])
                        self._or_status_var.set(f"{len(ids)} models loaded.")
                    else:
                        self._or_status_var.set("No models returned.")

                self._dlg.after(0, _apply)
            except Exception as e:
                self._dlg.after(
                    0, lambda: self._or_status_var.set(f"Fetch failed: {e}")
                )

        threading.Thread(target=_fetch, daemon=True).start()

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

        label_to_key = {
            "OpenRouter": "openrouter",
            "OpenAI": "openai",
            "NVIDIA AI": "nvidia",
            "Custom": "custom",
        }
        if "api_provider" in result:
            result["api_provider"] = label_to_key.get(
                result["api_provider"], str(result["api_provider"]).lower()
            )

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
            new_cfg["grid_padding"] = float(new_cfg.get("grid_padding", 6.0))
        except (ValueError, TypeError):
            new_cfg["grid_padding"] = DEFAULT_CONFIG["grid_padding"]

        save_config(new_cfg)
        self._config.update(new_cfg)
        self._dlg.destroy()


# ════════════════════════════════════════════════════════════════════════════
#  Utility functions
# ════════════════════════════════════════════════════════════════════════════


def _write_json(data: dict) -> None:
    """Atomically write config.json."""
    CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp = CONFIG_FILE.with_suffix(".tmp")
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        tmp.replace(CONFIG_FILE)
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


def _open_url(url: str) -> None:
    import webbrowser

    if url:
        try:
            webbrowser.open(url)
        except Exception:
            pass