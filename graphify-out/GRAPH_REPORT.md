# Graph Report - .  (2026-08-09)

## Corpus Check
- Corpus is ~33,568 words - fits in a single context window. You may not need a graph.

## Summary
- 308 nodes · 553 edges · 10 communities (9 shown, 1 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 3 edges (avg confidence: 0.65)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Docking Execution and RMSD
- Desktop GUI and Logging
- Application Configuration
- Target Context Parsing
- Ligand Identity and Receptor Prep
- Interaction Analysis
- PDBQT Preparation
- Graphify Knowledge Pipeline
- Docking Protocol Documentation
- Licensing and Distribution

## God Nodes (most connected - your core abstractions)
1. `DockingApp` - 41 edges
2. `run_pipeline()` - 28 edges
3. `SettingsDialog` - 20 edges
4. `parse_hetatm_records()` - 14 edges
5. `build_docking_target_context()` - 13 edges
6. `LigandInstanceKey` - 12 edges
7. `run_plip_interaction_analysis()` - 11 edges
8. `prepare_ligand_pdbqt()` - 10 edges
9. `run_vina_docking()` - 9 edges
10. `iter_hetatm_atoms()` - 9 edges

## Surprising Connections (you probably didn't know these)
- `LogCollector` --uses--> `SettingsDialog`  [INFERRED]
  main.py → settings.py
- `DockingApp` --uses--> `SettingsDialog`  [INFERRED]
  main.py → settings.py
- `Python Runtime Dependencies` --conceptually_related_to--> `Advanced Molecular Docking Protocol`  [INFERRED]
  requirements.txt → README.md
- `run_pipeline()` --calls--> `parse_hetatm_records()`  [EXTRACTED]
  pipeline.py → hetatm_parser.py
- `run_pipeline()` --calls--> `build_docking_target_context()`  [EXTRACTED]
  pipeline.py → hetatm_parser.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **AMDP Nine-Stage Workflow** — readme_active_site_selection, readme_receptor_cleaning_grid_box, readme_receptor_pdbqt_preparation, readme_native_ligand_preparation, readme_input_ligand_preparation, readme_validation_redocking, readme_production_docking, readme_summary_report [EXTRACTED 1.00]
- **Graphify Build Pipeline** — _agents_skills_graphify_skill_structural_extraction, _agents_skills_graphify_skill_semantic_extraction, _agents_skills_graphify_skill_community_detection, _agents_skills_graphify_skill_graph_health_check, _agents_skills_graphify_references_exports_graph_exports [EXTRACTED 1.00]

## Communities (10 total, 1 thin omitted)

### Community 0 - "Docking Execution and RMSD"
Cohesion: 0.06
Nodes (59): calculate_redock_rmsd(), generate_summary_report(), _optimal_rmsd(), _parse_pdbqt_model1_heavy_coords(), _parse_vina_scores(), Path, Molecular Docking Pipeline ===================================== docking.py —…, Calculate the RMSD between the crystal-pose native ligand and the top-ranked… (+51 more)

### Community 1 - "Desktop GUI and Logging"
Cohesion: 0.06
Nodes (18): DockingApp, LogCollector, main(), Any, Frame, Path, Tk, Copy selected log text while keeping the widget read-only. (+10 more)

### Community 2 - "Application Configuration"
Cohesion: 0.09
Nodes (28): build_chain_ligand_map(), Group candidate ligands by chain. Returns {chain_id: [list of candidate hetatm…, Group candidate ligands by chain. Returns {chain_id: [list of candidate hetatm…, _open_folder(), Molecular Docking Pipeline ===================================== main.py —…, Notebook, _initial_dir(), list_profiles() (+20 more)

### Community 3 - "Target Context Parsing"
Cohesion: 0.08
Nodes (31): build_docking_target_context(), DockingTargetContext, _float_field(), iter_hetatm_atoms(), LigandAtomReference, Path, Molecular Docking Pipeline ========================== hetatm_parser.py — HETATM…, Yield one normalized dict per HETATM atom in *structure_file*. Format is chosen… (+23 more)

### Community 4 - "Ligand Identity and Receptor Prep"
Cohesion: 0.11
Nodes (21): LigandInstanceKey, parse_hetatm_records(), Read a receptor structure file (PDB or mmCIF) and return one dict per unique…, Full structural identity of one selected ligand conformer., Full structural identity of one selected ligand conformer., Read a receptor structure file (PDB or mmCIF) and return one dict per unique…, clean_receptor_with_pymol(), extract_native_ligand_coords() (+13 more)

### Community 5 - "Interaction Analysis"
Cohesion: 0.16
Nodes (24): Element, _append_value(), _build_complex_pdb(), combine_interaction_csvs(), _find_plip_report(), _flatten_xml(), _format_score(), _is_docked_ligand_site() (+16 more)

### Community 6 - "PDBQT Preparation"
Cohesion: 0.14
Nodes (24): _build_cmd(), ensure_3d_sdf(), _ensure_explicit_hs(), _meeko_hint(), prepare_ligand_pdbqt(), prepare_receptor_pdbqt(), CompletedProcess, Path (+16 more)

### Community 7 - "Graphify Knowledge Pipeline"
Cohesion: 0.15
Nodes (14): Folder Watcher, URL Ingestion, Graph Export Formats, Provenance Confidence Model, Cross-Repository Graph Merge, Post-Commit Graph Rebuild, Graph-Aware Query Traversal, Media Transcription (+6 more)

### Community 8 - "Docking Protocol Documentation"
Cohesion: 0.15
Nodes (14): Deterministic Active-Site Selection, Advanced Molecular Docking Protocol, AutoDock Vina, Input Ligand Preparation, Meeko, Native Ligand Preparation, Production Docking, PyMOL (+6 more)

## Knowledge Gaps
- **15 isolated node(s):** `Native Ligand Preparation`, `Docking Summary Report`, `AutoDock Vina`, `PyMOL`, `Commercial License` (+10 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `DockingApp` connect `Desktop GUI and Logging` to `Application Configuration`?**
  _High betweenness centrality (0.216) - this node is a cross-community bridge._
- **Why does `run_pipeline()` connect `Docking Execution and RMSD` to `Desktop GUI and Logging`, `Application Configuration`, `Target Context Parsing`, `Ligand Identity and Receptor Prep`, `Interaction Analysis`, `PDBQT Preparation`?**
  _High betweenness centrality (0.156) - this node is a cross-community bridge._
- **Why does `SettingsDialog` connect `Application Configuration` to `Desktop GUI and Logging`?**
  _High betweenness centrality (0.117) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `SettingsDialog` (e.g. with `DockingApp` and `LogCollector`) actually correct?**
  _`SettingsDialog` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Native Ligand Preparation`, `Docking Summary Report`, `AutoDock Vina` to the rest of the system?**
  _15 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Docking Execution and RMSD` be split into smaller, more focused modules?**
  _Cohesion score 0.05628415300546448 - nodes in this community are weakly interconnected._
- **Should `Desktop GUI and Logging` be split into smaller, more focused modules?**
  _Cohesion score 0.062146892655367235 - nodes in this community are weakly interconnected._