# AI-Guided Molecular Docking Protocol (AGDP)

An automated, end-to-end virtual screening pipeline that takes a raw receptor structure (**PDB or mmCIF**) and one or more ligand SDF files, identifies the active site deterministically, and runs a fully automated **AutoDock Vina** docking workflow — orchestrated through a clean Tkinter GUI.

---

## Overview

AGDP runs a 9-stage pipeline (Stages 0–8):

1. **Stage 0 — Output setup** — a timestamped run directory is created inside your output folder
2. **Stage 1 — Active-site selection** — HETATM records are parsed from the structure (pure Python) and the co-crystallised ligand is chosen from the list. Selection is deterministic and user-confirmed — **no AI / API key required**
3. **Stage 2 — Receptor cleaning + grid box** — PyMOL removes solvent, ions and other chains; the native ligand's crystal coordinates are extracted directly from the original structure (pure Python) and used to compute the docking box
4. **Stage 3 — Receptor PDBQT preparation** — Meeko converts the cleaned receptor to PDBQT
5. **Stage 4 — Native ligand 3D retrieval & PDBQT** — a clean 3D structure of the native ligand is fetched from RCSB PDB and prepared for redocking
6. **Stage 5 — Input ligand PDBQT preparation** — Meeko processes each input SDF (2D SDFs are embedded via ETKDGv3 + MMFF94) with explicit hydrogen addition
7. **Stage 6 — Validation redocking** — the native ligand is re-docked into the pocket; RMSD against the crystal pose is reported
8. **Stage 7 — Production docking** — all input ligands are docked and ranked by binding affinity
9. **Stage 8 — Summary report** — a plain-text report is written with all scores, RMSD, and file locations

---

## Features

- **PDB and mmCIF input** — receptors can be `.pdb`, `.cif` or `.mmcif` (mmCIF is parsed by a minimal, dependency-free parser; the cleaned receptor is always written as PDB downstream)
- **Deterministic active-site selection** — HETATM records are read directly from the structure and you pick the ligand; no LLM calls, no API key
- **Pure-Python native-ligand extraction** — crystal coordinates are read from the original structure for grid centring and RMSD reference
- **Automatic ligand fetching** — retrieves clean 3D structures from RCSB PDB and PubChem for validation redocking
- **RMSD validation** — Hungarian-algorithm optimal atom matching (scipy) with a direct fallback; PASS / MARGINAL / FAIL thresholds
- **Intelligent receptor prep** — preserves essential cofactors (NAD, FAD, HEM, PLP, CoA) while removing solvent, ions, and other chains
- **Live GUI log** — colour-coded, real-time pipeline output with progress bar (ttkbootstrap themes when installed, graceful plain-Tk fallback)
- **Configurable docking parameters** — exhaustiveness, number of poses, CPU count, grid padding
- **Config profiles** — named `config_<name>.json` profiles alongside the default `config.json`
- **Stop at any stage** — the pipeline can be interrupted cleanly between steps

---

## Requirements

### Python Packages

```
pip install numpy scipy rdkit ttkbootstrap
```

(`ttkbootstrap` is optional — the GUI falls back to plain Tk if it isn't installed.)

### External Tools (paths configured via Settings)

| Tool | Purpose |
|---|---|
| [AutoDock Vina](https://vina.scripps.edu/) ≥ 1.2 | Molecular docking engine |
| [PyMOL](https://pymol.org/) | Receptor cleaning |
| [Meeko](https://github.com/forlilab/Meeko) ≥ 0.5 | PDBQT preparation (`mk_prepare_ligand`, `mk_prepare_receptor`) |

---

## Installation

```bash
git clone https://github.com/mustafajaved745-source/AI-Guided-Docking-Protocol-AGDP-.git
cd AI-Guided-Docking-Protocol-AGDP-
pip install numpy scipy rdkit ttkbootstrap
```

Then launch:

```bash
python main.py
```

On first launch, open **⚙ Settings** to configure the tool paths.

---

## Configuration

All settings are stored in `config.json` (auto-created with defaults on first run). Configure via the GUI Settings dialog, or edit directly:

```json
{
  "vina_exe":                "/path/to/vina",
  "pymol_exe":               "/path/to/pymol",
  "mk_prepare_ligand_cmd":   "/path/to/mk_prepare_ligand.py",
  "mk_prepare_receptor_cmd": "/path/to/mk_prepare_receptor.py",
  "output_dir":              "/path/to/docking_results",
  "grid_padding":            4.0,
  "exhaustiveness":          8,
  "num_modes":               9,
  "cpu":                     0
}
```

- `grid_padding` — docking-box padding per side in Å (default 4.0; minimum box size enforced at 12 Å per axis)
- `cpu` — `0` lets Vina auto-detect the core count
- **Named profiles** — a `config_<name>.json` placed next to `config.json` becomes selectable from the Settings dialog

---

## Usage

1. Launch `python main.py`
2. Click **⚙ Settings**, set the tool paths, click **Save**
3. Browse for a **receptor structure** — PDB (`*.pdb`) or mmCIF (`*.cif`, `*.mmcif`) — containing a co-crystallised ligand
4. Review the detected HETATM ligands and confirm the active-site ligand code
5. Add one or more **ligand SDF** files
6. Click **▶ Run Docking Pipeline**
7. Monitor live progress in the log panel
8. Click **Open Output Folder** when complete

---

## Output Structure

Each run creates a timestamped folder inside your configured output directory:

```
docking_results/
└── docking_run_20240115_143022/
    ├── cleaned_receptor.pdb       # Protein-only structure (single chain)
    ├── cleaned_receptor.pdbqt     # Meeko-prepared receptor
    ├── native_ligand.sdf          # Native ligand in SDF format
    ├── grid_box.txt               # Vina grid configuration
    ├── redocking/
    │   ├── <LIG>_online.sdf       # Clean 3D structure from RCSB PDB / PubChem
    │   ├── <LIG>_3d.sdf           # 3D-embedded native ligand
    │   ├── <LIG>_native.pdbqt     # Prepared native ligand
    │   └── native_redock_poses.pdbqt
    ├── ligands/
    │   ├── <name>_3d.sdf          # 3D-embedded input ligands
    │   └── <name>.pdbqt           # Prepared input ligands
    ├── docking/
    │   ├── <name>_poses.pdbqt     # Docked poses per ligand
    │   └── <name>_vina.log        # Vina run logs
    ├── docking_summary.txt        # Full ranked results report
    └── MDP-log.txt                # UI / pipeline log
```

---

## Validation — Redocking RMSD Thresholds

| RMSD | Status | Interpretation |
|---|---|---|
| ≤ 2.0 Å | ✅ PASS | Pocket correctly centred; results are reliable |
| 2.0 – 3.0 Å | ⚠️ MARGINAL | Interpret with caution; check ligand code |
| > 3.0 Å | ❌ FAIL | Pocket may be mis-centred; verify and rerun |

---

## Project Structure

```
├── main.py           # Tkinter GUI and pipeline orchestration entry point
├── pipeline.py       # Stage orchestrator — coordinates all modules
├── hetatm_parser.py  # Dependency-free HETATM parsing (PDB / mmCIF) for active-site selection
├── mmcif_parser.py   # Minimal mmCIF atom-site parser and PDB/mmCIF detection
├── receptor_prep.py  # PyMOL receptor cleaning, native-ligand extraction, grid box
├── ligand_prep.py    # Meeko ligand + receptor PDBQT preparation, 3D embedding
├── docking.py        # AutoDock Vina execution, redock RMSD, summary report
├── settings.py       # Settings dialog, config.json + named-profile management
├── requirements.txt  # Python package dependencies
└── config.json       # Auto-generated user configuration (not committed)
```

---

## Receptor Requirements

- Standard **PDB** or **mmCIF** file with at least one co-crystallised **non-trivial HETATM** record
- Waters-only and apo structures are **not supported** — a reference ligand is needed to centre the docking grid
- Multi-chain structures are handled automatically; the chain containing the selected ligand is isolated
- If redocking RMSD is poor, check that the correct HETATM ligand was selected in Stage 1

---

## Troubleshooting

**"No HETATM records found"** — The receptor structure has no ligand. Use a holo (ligand-bound) structure from the PDB.

**PyMOL fails with unknown option** — AGDP automatically tries `-cq`, `-c`, and no flags in sequence. Ensure PyMOL is callable from the configured path.

**Meeko fails with sanitization error** — The input SDF may have invalid valences. Pre-process with OpenBabel: `obabel input.sdf -O fixed.sdf -h`

**RMSD > 3.0 Å** — The active-site ligand may have been selected incorrectly, or the grid padding is too small. Check the HETATM table in the log and re-run with the correct ligand code.

**mmCIF file rejected** — Ensure the file uses a `.cif` / `.mmcif` extension (detection is by extension with a content sniff fallback).

**Plain (unthemed) GUI** — `ttkbootstrap` isn't installed; AGDP falls back to standard Tk. Install it with `pip install ttkbootstrap` for themed widgets.

---

## License

GNU Affero General Public License v3.0 — see `LICENSE` for details.

## Commercial License

The software in this repository is available under the AGPL v3 license for open-source and academic use.

For commercial use — see `COMMERCIAL` for details.
