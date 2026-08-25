# Graph Report - AMDP-Codex  (2026-08-09)

## Corpus Check
- 41 files · ~45,039 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 548 nodes · 900 edges · 27 communities (23 shown, 4 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 5 edges (avg confidence: 0.73)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- pipeline.py
- DockingApp
- SettingsDialog
- Knowledge Graph
- interaction_analysis.py
- hetatm_parser.py
- exporter.py
- What You Must Do When Invoked
- ligand_prep.py
- graphify reference: extra exports and benchmark
- graphify reference: query, path, explain
- graphify reference: add a URL and watch a folder
- graphify reference: commit hook and native CLAUDE.md integration
- graphify reference: incremental update and cluster-only
- graphify reference: GitHub clone and cross-repo merge
- graphify reference: transcribe video and audio
- AGENTS.md
- extraction-spec.md
- required
- properties
- docking.py
- required
- receptor_prep.py
- properties
- mmcif_parser.py
- parse_hetatm_records
- main.py

## God Nodes (most connected - your core abstractions)
1. `DockingApp` - 41 edges
2. `run_pipeline()` - 30 edges
3. `build_export()` - 22 edges
4. `SettingsDialog` - 20 edges
5. `ExportOptions` - 17 edges
6. `export_experiment()` - 15 edges
7. `parse_hetatm_records()` - 13 edges
8. `parse_interactions()` - 13 edges
9. `build_docking_target_context()` - 12 edges
10. `What You Must Do When Invoked` - 12 edges

## Surprising Connections (you probably didn't know these)
- `LogCollector` --uses--> `SettingsDialog`  [INFERRED]
  main.py → settings.py
- `DockingApp` --uses--> `SettingsDialog`  [INFERRED]
  main.py → settings.py
- `run_pipeline()` --calls--> `run_vina_docking()`  [EXTRACTED]
  pipeline.py → docking.py
- `run_pipeline()` --calls--> `calculate_redock_rmsd()`  [EXTRACTED]
  pipeline.py → docking.py
- `run_pipeline()` --calls--> `generate_summary_report()`  [EXTRACTED]
  pipeline.py → docking.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Extraction Pipeline** — _claude_skills_graphify_skill_ast_extraction, _claude_skills_graphify_skill_semantic_extraction, _claude_skills_graphify_skill_corpus_detection [INFERRED 0.95]
- **Knowledge Graph Exports** — _claude_skills_graphify_references_exports_neo4j_export, _claude_skills_graphify_references_exports_falkordb_export, _claude_skills_graphify_references_exports_wiki_export, _claude_skills_graphify_references_exports_mcp_server, _claude_skills_graphify_skill_html_visualization [INFERRED 0.85]
- **Graph Navigation** — _claude_skills_graphify_references_query_query, _claude_skills_graphify_references_query_path, _claude_skills_graphify_references_query_explain [INFERRED 0.95]

## Communities (27 total, 4 thin omitted)

### Community 0 - "pipeline.py"
Cohesion: 0.09
Nodes (37): _adaptive_exhaustiveness(), _compare_redock_atom_counts(), _detect_preserved_cofactor(), _download_and_process_sdf(), _extract_cid_from_filename(), _extract_sdf_name(), fetch_native_ligand_sdf(), _fetch_rcsb_chem_metadata() (+29 more)

### Community 1 - "DockingApp"
Cohesion: 0.06
Nodes (20): DockingApp, LogCollector, main(), _open_folder(), Any, Frame, Path, Tk (+12 more)

### Community 2 - "SettingsDialog"
Cohesion: 0.10
Nodes (21): Notebook, _initial_dir(), list_profiles(), _profile_path(), Frame, Path, Tk, Molecular Docking Pipeline ===================================== settings.py —… (+13 more)

### Community 3 - "Knowledge Graph"
Cohesion: 0.06
Nodes (38): Graphify Skill Trigger, File Watcher, URL Ingestion, Token Reduction Benchmark, FalkorDB Export, MCP Server, Neo4j Export, Wiki Export (+30 more)

### Community 4 - "interaction_analysis.py"
Cohesion: 0.16
Nodes (24): Element, _append_value(), _build_complex_pdb(), combine_interaction_csvs(), _find_plip_report(), _flatten_xml(), _format_score(), _is_docked_ligand_site() (+16 more)

### Community 5 - "hetatm_parser.py"
Cohesion: 0.14
Nodes (16): build_docking_target_context(), DockingTargetContext, _float_field(), iter_hetatm_atoms(), LigandAtomReference, Path, Molecular Docking Pipeline ========================== hetatm_parser.py — HETATM…, Yield one normalized dict per HETATM atom in *structure_file*. Format is chosen… (+8 more)

### Community 6 - "exporter.py"
Cohesion: 0.13
Nodes (43): _append(), build_export(), classify_hbond(), _clean(), _codebook(), _coordinates(), _experiment(), export_experiment() (+35 more)

### Community 7 - "What You Must Do When Invoked"
Cohesion: 0.08
Nodes (24): For /graphify add and --watch, For /graphify query, For the commit hook and native AGENTS.md integration, For --update and --cluster-only, /graphify, Honesty Rules, Interpreter guard for subcommands, Part A - Structural extraction for code files (+16 more)

### Community 8 - "ligand_prep.py"
Cohesion: 0.14
Nodes (24): _build_cmd(), ensure_3d_sdf(), _ensure_explicit_hs(), _meeko_hint(), prepare_ligand_pdbqt(), prepare_receptor_pdbqt(), CompletedProcess, Path (+16 more)

### Community 9 - "graphify reference: extra exports and benchmark"
Cohesion: 0.22
Nodes (8): graphify reference: extra exports and benchmark, Step 6b - Wiki (only if --wiki flag), Step 7 - Neo4j export (only if --neo4j or --neo4j-push flag), Step 7a - FalkorDB export (only if --falkordb or --falkordb-push flag), Step 7b - SVG export (only if --svg flag), Step 7c - GraphML export (only if --graphml flag), Step 7d - MCP server (only if --mcp flag), Step 8 - Token reduction benchmark (only if total_words > 5000)

### Community 10 - "graphify reference: query, path, explain"
Cohesion: 0.33
Nodes (5): For /graphify explain, For /graphify path, graphify reference: query, path, explain, Step 0 — Constrained query expansion (REQUIRED before traversal), Step 1 — Traversal

### Community 11 - "graphify reference: add a URL and watch a folder"
Cohesion: 0.50
Nodes (3): For /graphify add, For --watch, graphify reference: add a URL and watch a folder

### Community 12 - "graphify reference: commit hook and native CLAUDE.md integration"
Cohesion: 0.50
Nodes (3): For git commit hook, For native CLAUDE.md integration, graphify reference: commit hook and native CLAUDE.md integration

### Community 13 - "graphify reference: incremental update and cluster-only"
Cohesion: 0.50
Nodes (3): For --cluster-only, For --update (incremental re-extraction), graphify reference: incremental update and cluster-only

### Community 18 - "required"
Cohesion: 0.05
Nodes (45): best, counts, flags, quality, rank, ref, residues, scores (+37 more)

### Community 19 - "properties"
Cohesion: 0.04
Nodes (49): center, compact, contact_residues, contacts, engine, export_detail, full, HB_quality_rules (+41 more)

### Community 20 - "docking.py"
Cohesion: 0.14
Nodes (20): calculate_redock_rmsd(), generate_summary_report(), _optimal_rmsd(), _parse_pdbqt_model1_heavy_coords(), _parse_vina_scores(), Path, Molecular Docking Pipeline ===================================== docking.py —…, Calculate the RMSD between the crystal-pose native ligand and the top-ranked… (+12 more)

### Community 21 - "required"
Cohesion: 0.10
Nodes (20): codebook, experiment, global_analysis, ligands, provenance, ranking, receptor, reference (+12 more)

### Community 22 - "receptor_prep.py"
Cohesion: 0.16
Nodes (17): calculate_grid_box(), clean_receptor_with_pymol(), extract_native_ligand_coords(), CompletedProcess, Path, Molecular Docking Pipeline =====================================…, Return the heavy-atom crystal coordinates of the selected native ligand, parsed…, Compute the AutoDock Vina docking grid box from native-ligand heavy-atom… (+9 more)

### Community 23 - "properties"
Cohesion: 0.15
Nodes (13): outside_preferred_geometry, preferred_max_DA_A, preferred_min_angle_deg, properties, required, type, type, type (+5 more)

### Community 24 - "mmcif_parser.py"
Cohesion: 0.23
Nodes (11): _flush(), _normalize(), _normalize_row(), Molecular Docking Pipeline ========================== mmcif_parser.py —…, Yield every mmCIF token in the file, comments removed., Map an mmCIF '.' / '?' missing-value marker to an empty string., Emit normalized atom dicts for a completed _atom_site loop chunk., Map a raw `_atom_site` row to the pipeline's normalized atom dict. (+3 more)

### Community 25 - "parse_hetatm_records"
Cohesion: 0.28
Nodes (4): LigandInstanceKey, parse_hetatm_records(), Read a receptor structure file (PDB or mmCIF) and return one dict per unique…, Full structural identity of one selected ligand conformer.

### Community 26 - "main.py"
Cohesion: 0.33
Nodes (4): DockLLM export, build_chain_ligand_map(), Group candidate ligands by chain. Returns {chain_id: [list of candidate hetatm…, Molecular Docking Pipeline ===================================== main.py —…

## Knowledge Gaps
- **134 isolated node(s):** `$schema`, `$id`, `title`, `type`, `schema` (+129 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `DockingApp` connect `DockingApp` to `main.py`, `SettingsDialog`, `exporter.py`?**
  _High betweenness centrality (0.077) - this node is a cross-community bridge._
- **Why does `run_pipeline()` connect `pipeline.py` to `DockingApp`, `interaction_analysis.py`, `hetatm_parser.py`, `exporter.py`, `ligand_prep.py`, `docking.py`, `receptor_prep.py`, `parse_hetatm_records`, `main.py`?**
  _High betweenness centrality (0.062) - this node is a cross-community bridge._
- **Why does `SettingsDialog` connect `SettingsDialog` to `DockingApp`, `main.py`?**
  _High betweenness centrality (0.042) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `SettingsDialog` (e.g. with `DockingApp` and `LogCollector`) actually correct?**
  _`SettingsDialog` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `$schema`, `$id`, `title` to the rest of the system?**
  _134 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `pipeline.py` be split into smaller, more focused modules?**
  _Cohesion score 0.08677098150782361 - nodes in this community are weakly interconnected._
- **Should `DockingApp` be split into smaller, more focused modules?**
  _Cohesion score 0.060655737704918035 - nodes in this community are weakly interconnected._