# Advanced Molecular Docking Protocol (AMDP)

Desktop Tkinter application for reproducible AutoDock Vina screening. It accepts a ligand-bound receptor structure (`.pdb`, `.cif`, or `.mmcif`) and one or more SDF ligands, prepares inputs, validates a reference-ligand redock, ranks production dockings, and writes PLIP interaction and DockLLM analysis files.

## What it does

AMDP runs these stages in one timestamped output folder:

0. Creates run folders and records selected settings.
1. Parses receptor HETATM records and requires you to choose exact chain and native-ligand instance. No LLM or API key is used for site selection.
2. Cleans selected receptor chain with PyMOL, optionally preserves essential cofactor, extracts crystal-ligand coordinates, and builds validated Vina grid.
3. Prepares receptor PDBQT with Meeko.
4. Fetches native ligand 3D structure from PubChem/RCSB/NIH NCI fallbacks, then prepares it for redocking.
5. Ensures input SDFs have 3D coordinates (ETKDGv3 + MMFF94 when needed), adds hydrogens, and prepares ligand PDBQT files.
6. Redocks native ligand and reports RMSD against crystal coordinates when validation preparation succeeds.
7. Docks supplied ligands with AutoDock Vina and ranks successful results by affinity.
8. Runs PLIP for validation and production poses, producing per-ligand XML, text, and CSV interaction reports plus combined CSV.
9. Writes docking summary, complete UI log, DockLLM source manifest, and compact LLM-analysis JSON.

## Requirements

### Python

- Python 3.10+ with Tk support
- Packages in [`requirements.txt`](requirements.txt): NumPy, SciPy, RDKit, and `ttkbootstrap`

```bash
python -m pip install -r requirements.txt
```

`ttkbootstrap` provides themed widgets; GUI falls back to standard Tk if unavailable.

### External tools

All tools below are required by settings validation before a run.

| Tool | Used for |
| --- | --- |
| [AutoDock Vina](https://vina.scripps.edu/) | Redocking and production docking |
| [PyMOL](https://pymol.org/) | Receptor cleanup and PLIP complex creation |
| [Meeko](https://github.com/forlilab/Meeko) | `mk_prepare_ligand` and `mk_prepare_receptor` PDBQT preparation |
| [PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) | Protein-ligand interaction reports |

Set each executable/script path in app Settings. Paths can point to a binary, a Python script, or a command available on `PATH` where supported.

## Install and launch

```bash
git clone https://github.com/mustafajaved745-source/Advanced-Molecular-Docking-Protocol-AMDP-.git
cd Advanced-Molecular-Docking-Protocol-AMDP-
python -m pip install -r requirements.txt
python main.py
```

On first launch, open **Settings** and configure Vina, PyMOL, both Meeko commands, PLIP, output directory, and docking parameters. `config.json` is created beside `settings.py`; named profiles are saved as `config_<name>.json`.

## Configuration

```json
{
  "vina_exe": "/path/to/vina",
  "pymol_exe": "/path/to/pymol",
  "mk_prepare_ligand_cmd": "/path/to/mk_prepare_ligand.py",
  "mk_prepare_receptor_cmd": "/path/to/mk_prepare_receptor.py",
  "plip_cmd": "/path/to/plip",
  "output_dir": "/path/to/docking_results",
  "grid_padding": 4.0,
  "exhaustiveness": 8,
  "num_modes": 9,
  "cpu": 0
}
```

- `grid_padding`: clearance around each native-ligand axis in Å; grid enforces at least 12 Å per axis and validates reference-atom clearance.
- `exhaustiveness`: base Vina search effort. AMDP scales it by grid volume, capped at 64.
- `num_modes`: maximum Vina poses per ligand.
- `cpu`: `0` lets Vina choose CPU count.
- Profiles: create/select `config_<name>.json` from Settings to keep tool paths and docking parameters separate.

Do not commit machine-specific `config.json` files with private paths.

## Run docking experiment

1. Start `python main.py`.
2. Open **Settings**, complete all tool paths, choose output folder, then save.
3. Choose receptor `.pdb`, `.cif`, or `.mmcif` file.
4. Select target chain, then select exact co-crystallized ligand instance shown by app.
5. Add one or more ligand `.sdf` files.
6. Click **Run Docking Pipeline**. Live log, stage progress, and per-ligand status update while background thread runs.
7. Review `docking_summary.txt`, redocking RMSD, ranked affinities, PLIP reports, and DockLLM JSON in run folder.

You can request clean stop; pipeline stops after current step.

## Input rules and behavior

- Receptor must contain non-trivial co-crystallized ligand HETATM records. Apo structures and water-only structures cannot define reference grid.
- AMDP supports PDB and mmCIF receptors. Downstream cleaned receptor is PDB.
- Chain and ligand instance are user-selected and revalidated before docking. Repeated ligand codes are disambiguated by structural instance.
- Essential cofactors such as NAD, FAD, HEM, PLP, and CoA can be preserved during receptor cleanup when detected.
- If native-ligand retrieval/preparation fails, native redocking is skipped; production docking can still proceed.
- PLIP failures for individual poses are recorded in log/summary without discarding successful docking results.

## Output

Each run is saved as `docking_run_YYYYMMDD_HHMMSS` under configured output directory:

```text
docking_run_YYYYMMDD_HHMMSS/
├── cleaned_receptor.pdb
├── cleaned_receptor.pdbqt
├── grid_box.txt
├── redocking/
│   ├── <LIG>_online.sdf
│   ├── <LIG>_3d.sdf
│   ├── <LIG>_native.pdbqt
│   └── native_redock_poses.pdbqt
├── ligands/
│   ├── <name>_3d.sdf
│   └── <name>.pdbqt
├── docking/
│   ├── <name>_poses.pdbqt
│   └── <name>_vina.log
├── interactions/
│   ├── <ligand>/                 # PLIP complex, XML, text report, CSV
│   └── interaction_summary.csv
├── docking_summary.txt
├── MDP-log.txt
├── dockllm-source.json
└── docking_experiment.llm.json
```

### Redocking RMSD

| RMSD | Status | Meaning |
| --- | --- | --- |
| ≤ 2.0 Å | PASS | Reference pose reproduced well |
| 2.0–3.0 Å | MARGINAL | Check ligand selection and grid before interpreting results |
| > 3.0 Å | FAIL | Pocket may be mis-centred |

RMSD compares heavy atoms using optimal matching; atom-count, protonation, or tautomer differences can make it unreliable.

## DockLLM export

At pipeline completion AMDP writes `dockllm-source.json` and compact `docking_experiment.llm.json`. From GUI, use **Export LLM JSON** to regenerate analysis file from source manifest.

CLI export:

```bash
python -m llm_export \
  --experiment /path/to/dockllm-source.json \
  --output docking_experiment.llm.json \
  --detail standard
```

Use `--detail full|standard|compact`, plus `--pretty`, `--include-coordinates`, or `--include-all-pose-interactions` when needed.

## Project layout

```text
main.py                  Tkinter UI, file selection, background-run orchestration
pipeline.py              Stages 0–9, online ligand lookup, output writing
hetatm_parser.py         PDB/mmCIF HETATM parsing and exact target context
mmcif_parser.py          Minimal mmCIF atom-site parsing
receptor_prep.py         PyMOL cleanup, native coordinates, grid calculation
ligand_prep.py           RDKit 3D preparation and Meeko PDBQT conversion
docking.py               Vina execution, RMSD, text summary
interaction_analysis.py  PLIP complex/report/CSV generation
settings.py              Configuration and profile dialog
llm_export/              DockLLM JSON builder and CLI
requirements.txt         Python dependencies
```

## Troubleshooting

- **No HETATM records / no active site:** use ligand-bound receptor containing real co-crystallized ligand.
- **Tool path validation fails:** set Vina, PyMOL, both Meeko commands, and PLIP in Settings; verify file exists or command is on `PATH`.
- **PyMOL unknown option:** app retries compatible headless flag forms. Confirm configured PyMOL executable runs outside app.
- **Meeko ligand error:** SDF may have invalid valence/bonds. Repair molecule and retry; inspect UI log.
- **Redocking skipped or RMSD poor:** confirm exact ligand instance, chain, ligand protonation/tautomer, and grid padding.
- **PLIP failed:** docking output remains usable; inspect `MDP-log.txt` and interaction error in summary, then confirm PLIP/PyMOL paths.
- **Plain GUI:** install `ttkbootstrap` in same Python environment used for `python main.py`.

## License

GNU Affero General Public License v3.0. See [LICENSE](LICENSE). Commercial use that cannot comply with AGPL requires separate terms; see [COMMERCIAL.md](COMMERCIAL.md).
