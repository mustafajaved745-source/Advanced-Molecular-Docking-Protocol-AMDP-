# AI-Guided Molecular Docking Protocol (AGDP)

An automated, end-to-end virtual screening pipeline that combines **AI-powered active-site identification** with **AutoDock Vina** molecular docking — orchestrated through a clean Tkinter GUI.

---

## Overview

AGDP takes a raw receptor PDB and one or more ligand SDF files, then runs a fully automated 8-stage pipeline:

1. **AI Active-Site Identification** — an LLM analyses the receptor's HETATM records and identifies the biologically relevant co-crystallised ligand
2. **Receptor Cleaning** — PyMOL removes solvent, ions, and buffer molecules; isolates the target chain
3. **Grid Box Calculation** — docking box is centred on the native ligand's centroid with configurable padding
4. **Receptor PDBQT Preparation** — Meeko converts the cleaned PDB to PDBQT
5. **Ligand PDBQT Preparation** — Meeko processes each input SDF with explicit hydrogen addition
6. **Validation Redocking** — the native ligand is re-docked back into the pocket; RMSD is reported
7. **Production Docking** — all input ligands are docked and ranked by binding affinity
8. **Summary Report** — a plain-text report is written with all scores, RMSD, and file locations

---

## Features

- **AI provider-agnostic** — works with OpenRouter, OpenAI, NVIDIA AI Endpoints, or any OpenAI-compatible API
- **Automatic ligand fetching** — retrieves clean 3D structures from RCSB PDB and PubChem for validation redocking
- **Intelligent receptor prep** — preserves essential cofactors (NAD, FAD, HEM, PLP, CoA) while removing noise
- **RMSD validation** — uses Hungarian-algorithm optimal atom matching (scipy) with a direct fallback
- **Live GUI log** — colour-coded, real-time pipeline output with progress bar
- **Configurable docking parameters** — exhaustiveness, number of poses, CPU count, grid padding
- **Stop at any stage** — pipeline can be interrupted cleanly between steps

---

## Requirements

### Python Packages
```
pip install openai numpy scipy rdkit
```

### External Tools (paths configured via Settings)

| Tool | Purpose |
|---|---|
| [AutoDock Vina](https://vina.scripps.edu/) ≥ 1.2 | Molecular docking engine |
| [PyMOL](https://pymol.org/) | Receptor cleaning and ligand extraction |
| [Meeko](https://github.com/forlilab/Meeko) ≥ 0.5 | PDBQT preparation (`mk_prepare_ligand`, `mk_prepare_receptor`) |

---

## Installation

```bash
git clone <git clone https://github.com/mustafajaved745-source/AI-Guided-Docking-Protocol-AGDP-.git>
cd agdp
pip install openai numpy scipy rdkit
```

Then launch:
```bash
python main.py
```

On first launch, open **⚙ Settings** to configure tool paths and your API key.

---

## Configuration

All settings are stored in `config.json` (auto-created on first run). Configure via the GUI Settings dialog, or edit directly:

```json
{
  "vina_exe":               "/path/to/vina",
  "pymol_exe":              "/path/to/pymol",
  "mk_prepare_ligand_cmd":  "/path/to/mk_prepare_ligand.py",
  "mk_prepare_receptor_cmd":"/path/to/mk_prepare_receptor.py",
  "api_key":                "your-api-key-here",
  "api_provider":           "openrouter",
  "api_base_url":           "https://openrouter.ai/api/v1",
  "model":                  "meta-llama/llama-3.1-70b-instruct",
  "output_dir":             "/path/to/docking_results",
  "grid_padding":           6.0,
  "exhaustiveness":         8,
  "num_modes":              9,
  "cpu":                    0
}
```

### API Providers

| Provider | Base URL | Notes |
|---|---|---|
| OpenRouter | `https://openrouter.ai/api/v1` | Free & paid models; recommended |
| OpenAI | `https://api.openai.com/v1` | GPT-4o-mini works well |
| NVIDIA AI | `https://integrate.api.nvidia.com/v1` | Llama models |
| Custom | Any OpenAI-compatible URL | Self-hosted LLMs supported |

---

## Usage

1. Launch `python main.py`
2. Click **⚙ Settings**, set all tool paths and API key, click **Save**
3. Browse for a **Receptor PDB** file (must contain a co-crystallised ligand)
4. Add one or more **Ligand SDF** files
5. Click **▶ Run Docking Pipeline**
6. Monitor live progress in the log panel
7. Click **Open Output Folder** when complete

---

## Output Structure

Each run creates a timestamped folder inside your configured output directory:

```
docking_results/
└── docking_run_20240115_143022/
    ├── cleaned_receptor.pdb       # Protein-only structure (single chain)
    ├── cleaned_receptor.pdbqt     # Meeko-prepared receptor
    ├── native_ligand.pdb          # Extracted crystal ligand
    ├── native_ligand.sdf          # Native ligand in SDF format
    ├── grid_box.txt               # Vina grid configuration
    ├── redocking/
    │   ├── <LIG>_online.sdf       # Clean 3D structure from RCSB/PubChem
    │   ├── <LIG>_native.pdbqt     # Prepared native ligand
    │   └── native_redock_poses.pdbqt
    ├── ligands/
    │   └── <name>.pdbqt           # Prepared input ligands
    ├── docking/
    │   └── <name>_poses.pdbqt     # Docked poses per ligand
    └── docking_summary.txt        # Full ranked results report
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
├── ai_identify.py    # HETATM parsing and LLM-based ligand identification
├── receptor_prep.py  # PyMOL receptor cleaning and grid box calculation
├── ligand_prep.py    # Meeko ligand and receptor PDBQT preparation
├── docking.py        # AutoDock Vina execution, RMSD, and summary report
├── settings.py       # Settings dialog and config.json management
└── config.json       # Auto-generated user configuration
```

---

## Receptor Requirements

- Must be a standard PDB file with at least one co-crystallised **non-trivial HETATM** record
- Waters-only and apo structures are **not supported** — the AI needs a reference ligand to centre the docking grid
- Multi-chain structures are handled automatically; the chain containing the identified ligand is isolated

---

## Troubleshooting

**"No HETATM records found"** — The receptor PDB has no ligand. Use a holo (ligand-bound) structure from the PDB.

**PyMOL fails with unknown option** — AGDP automatically tries `-cq`, `-c`, and no flags in sequence. Ensure PyMOL is callable from the configured path.

**Meeko fails with sanitization error** — The input SDF may have invalid valences. Pre-process with OpenBabel: `obabel input.sdf -O fixed.sdf -h`

**RMSD > 3.0 Å** — The AI may have identified the wrong ligand code, or the grid padding is too small. Check the HETATM table in the log and re-run with a corrected ligand code if needed.

**API key errors** — Verify your key is valid for the selected provider and has sufficient credits. Use the **⟳ Fetch** button to test model connectivity.

---

## License

GNU Affero General Public License v3.0 — see `LICENSE` for details.

## Commercial License

The software in this repository is available under the AGPL v3 license for open-source and academic use.

For commercial use — see `COMMERCIAL` for details.
