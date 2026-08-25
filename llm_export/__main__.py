"""Command-line DockLLM exporter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .exporter import ExportOptions, export_experiment


def main() -> int:
    parser = argparse.ArgumentParser(description="Export deterministic DockLLM JSON")
    parser.add_argument("--experiment", required=True, help="DockLLM source manifest JSON")
    parser.add_argument("--output", default="docking_experiment.llm.json")
    parser.add_argument("--detail", choices=("full", "standard", "compact"), default="standard")
    parser.add_argument("--pretty", action="store_true")
    parser.add_argument("--include-coordinates", action="store_true")
    parser.add_argument("--include-all-pose-interactions", action="store_true")
    args = parser.parse_args()
    context = json.loads(Path(args.experiment).read_text(encoding="utf-8"))
    export_experiment(context, args.output, ExportOptions(args.detail, args.pretty, args.include_coordinates, args.include_all_pose_interactions))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
