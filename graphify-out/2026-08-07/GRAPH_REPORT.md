# Graph Report - Pro-Dock  (2026-08-07)

## Corpus Check
- 44 files · ~36,221 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 408 nodes · 644 edges · 18 communities (17 shown, 1 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 5 edges (avg confidence: 0.73)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- pipeline.py
- DockingApp
- SettingsDialog
- Graphify
- ligand_prep.py
- receptor_prep.py
- Knowledge Graph
- compress.py
- validate.py
- caveman-compress/README.md
- cavecrew/SKILL.md
- Caveman Help
- Caveman Compress
- caveman/SKILL.md
- caveman-commit
- caveman-review
- caveman-stats
- __init__.py

## God Nodes (most connected - your core abstractions)
1. `DockingApp` - 39 edges
2. `run_pipeline()` - 25 edges
3. `SettingsDialog` - 20 edges
4. `validate()` - 14 edges
5. `compress_file()` - 12 edges
6. `Knowledge Graph` - 11 edges
7. `prepare_ligand_pdbqt()` - 10 edges
8. `Graphify` - 10 edges
9. `detect_file_type()` - 9 edges
10. `run_vina_docking()` - 9 edges

## Surprising Connections (you probably didn't know these)
- `LogCollector` --uses--> `SettingsDialog`  [INFERRED]
  main.py → settings.py
- `DockingApp` --uses--> `SettingsDialog`  [INFERRED]
  main.py → settings.py
- `run_pipeline()` --calls--> `parse_hetatm_records()`  [EXTRACTED]
  pipeline.py → hetatm_parser.py
- `run_pipeline()` --calls--> `prepare_ligand_pdbqt()`  [EXTRACTED]
  pipeline.py → ligand_prep.py
- `run_pipeline()` --calls--> `prepare_receptor_pdbqt()`  [EXTRACTED]
  pipeline.py → ligand_prep.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Extraction Pipeline** — _claude_skills_graphify_skill_ast_extraction, _claude_skills_graphify_skill_semantic_extraction, _claude_skills_graphify_skill_corpus_detection [INFERRED 0.95]
- **Knowledge Graph Exports** — _claude_skills_graphify_references_exports_neo4j_export, _claude_skills_graphify_references_exports_falkordb_export, _claude_skills_graphify_references_exports_wiki_export, _claude_skills_graphify_references_exports_mcp_server, _claude_skills_graphify_skill_html_visualization [INFERRED 0.85]
- **Graph Navigation** — _claude_skills_graphify_references_query_query, _claude_skills_graphify_references_query_path, _claude_skills_graphify_references_query_explain [INFERRED 0.95]

## Communities (18 total, 1 thin omitted)

### Community 0 - "pipeline.py"
Cohesion: 0.06
Nodes (59): calculate_redock_rmsd(), generate_summary_report(), _optimal_rmsd(), _parse_pdb_heavy_coords(), _parse_pdbqt_model1_heavy_coords(), _parse_vina_scores(), Path, Molecular Docking Pipeline ===================================== docking.py —… (+51 more)

### Community 1 - "DockingApp"
Cohesion: 0.07
Nodes (16): DockingApp, LogCollector, main(), Any, Frame, Path, Tk, A dark card container: tk.Frame with a 1px highlight border + title. Returns a… (+8 more)

### Community 2 - "SettingsDialog"
Cohesion: 0.08
Nodes (31): build_chain_ligand_map(), parse_hetatm_records(), Path, Molecular Docking Pipeline ========================== hetatm_parser.py — HETATM…, Group candidate ligands by chain. Returns {chain_id: [list of candidate hetatm…, Read a PDB file and return one dict per unique HETATM residue group. Each dict…, _open_folder(), Molecular Docking Pipeline ===================================== main.py —… (+23 more)

### Community 3 - "Graphify"
Cohesion: 0.09
Nodes (24): Graphify Skill Trigger, File Watcher, URL Ingestion, Token Reduction Benchmark, Confidence Score Rubric, Hyperedge, Semantic Similarity, CLAUDE.md Integration (+16 more)

### Community 4 - "ligand_prep.py"
Cohesion: 0.14
Nodes (23): _build_cmd(), ensure_3d_sdf(), _ensure_explicit_hs(), _meeko_hint(), prepare_ligand_pdbqt(), prepare_receptor_pdbqt(), CompletedProcess, Path (+15 more)

### Community 5 - "receptor_prep.py"
Cohesion: 0.17
Nodes (17): calculate_grid_box(), clean_receptor_with_pymol(), extract_native_ligand_pdb(), _parse_heavy_atom_coords(), CompletedProcess, Path, Molecular Docking Pipeline =====================================…, Extract the native ligand from the original receptor PDB using pure Python… (+9 more)

### Community 6 - "Knowledge Graph"
Cohesion: 0.15
Nodes (14): FalkorDB Export, MCP Server, Neo4j Export, Wiki Export, Cross-Repository Merge, GitHub Clone, Explain, Path (+6 more)

### Community 7 - "compress.py"
Cohesion: 0.12
Nodes (27): main(), print_usage(), backup_dir_for(), build_compress_prompt(), build_fix_prompt(), call_claude(), compress_file(), is_sensitive_path() (+19 more)

### Community 8 - "validate.py"
Cohesion: 0.16
Nodes (22): benchmark_pair(), count_tokens(), main(), print_table(), Path, count_bullets(), extract_code_blocks(), extract_headings() (+14 more)

### Community 9 - "caveman-compress/README.md"
Cohesion: 0.09
Nodes (20): Before / After, Benchmarks, How It Work, <img src="../../docs/assets/dancing-rock.svg" width="20" height="20" alt="rock"/> Caveman (285 tokens), Install, 📄 Original (706 tokens), Part of Caveman, Security (+12 more)

### Community 10 - "cavecrew/SKILL.md"
Cohesion: 0.14
Nodes (12): cavecrew, Example chaining, How to invoke, Model overrides, See also, What it does, Auto-clarity (inherited), Chaining patterns (+4 more)

### Community 11 - "Caveman Help"
Cohesion: 0.14
Nodes (12): caveman-help, Example output, How to invoke, See also, What it does, Caveman Help, Configure Default Mode, Deactivate (+4 more)

### Community 12 - "Caveman Compress"
Cohesion: 0.17
Nodes (11): Boundaries, Caveman Compress, Compress, Compression Rules, Pattern, Preserve EXACTLY (never modify), Preserve Structure, Process (+3 more)

### Community 13 - "caveman/SKILL.md"
Cohesion: 0.17
Nodes (10): caveman, Example output, How to invoke, See also, What it does, Auto-Clarity, Boundaries, Intensity (+2 more)

### Community 14 - "caveman-commit"
Cohesion: 0.18
Nodes (9): caveman-commit, Example output, How to invoke, See also, What it does, Auto-Clarity, Boundaries, Examples (+1 more)

### Community 15 - "caveman-review"
Cohesion: 0.18
Nodes (9): caveman-review, Example output, How to invoke, See also, What it does, Auto-Clarity, Boundaries, Examples (+1 more)

### Community 16 - "caveman-stats"
Cohesion: 0.29
Nodes (5): caveman-stats, Example output, How to invoke, See also, What it does

## Knowledge Gaps
- **87 isolated node(s):** `What it does`, `How to invoke`, `Example chaining`, `Model overrides`, `See also` (+82 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `DockingApp` connect `DockingApp` to `SettingsDialog`?**
  _High betweenness centrality (0.081) - this node is a cross-community bridge._
- **Why does `run_pipeline()` connect `pipeline.py` to `DockingApp`, `SettingsDialog`, `ligand_prep.py`, `receptor_prep.py`?**
  _High betweenness centrality (0.064) - this node is a cross-community bridge._
- **Why does `SettingsDialog` connect `SettingsDialog` to `DockingApp`?**
  _High betweenness centrality (0.049) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `SettingsDialog` (e.g. with `DockingApp` and `LogCollector`) actually correct?**
  _`SettingsDialog` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `What it does`, `How to invoke`, `Example chaining` to the rest of the system?**
  _87 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `pipeline.py` be split into smaller, more focused modules?**
  _Cohesion score 0.05792349726775956 - nodes in this community are weakly interconnected._
- **Should `DockingApp` be split into smaller, more focused modules?**
  _Cohesion score 0.06936026936026936 - nodes in this community are weakly interconnected._