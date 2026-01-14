#!/usr/bin/env python3
"""Batch runner for Gen_ICT across existing partition files."""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Optional

BINARY_PATH = "./Gen_ICT"

CSV_HEADER = (
    "label,rep,solver_input,"
    "num_vertices,num_components,min_component_size,max_component_size,"
    "num_trees,elapsed_ms,delay_ms,elapsed_sort_ms,elapsed_stock_ms,status"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Gen_ICT across a suite of partitions and append timing/statistics to a CSV."
        )
    )
    parser.add_argument(
        "--inputs",
        nargs="*",
        default=["tests/*.txt"],
        help="Glob pattern(s) for deterministic partition files (default: tests/*.txt).",
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=1,
        help="How many timed repetitions per input (default: 1).",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=0,
        help="Untimed warmup runs per input to stabilize caches (default: 0).",
    )
    parser.add_argument(
        "--output",
        default="results.csv",
        help="CSV file to append measurements to (default: results.csv).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress messages and only report errors.",
    )
    parser.add_argument(
        "--time-limit",
        type=float,
        default=None,
        help="Optional wall-clock limit (seconds) passed to the solver.",
    )
    return parser.parse_args()


def discover_inputs(patterns: Iterable[str]) -> List[str]:
    """Expand the provided glob patterns into a sorted list of files."""
    seen: List[str] = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if not matches:
            sys.stderr.write(f"[run_bench] warning: no files match pattern '{pattern}'.\n")
        seen.extend(matches)
    # Preserve order but drop duplicates
    unique: List[str] = []
    for path in seen:
        if path not in unique:
            unique.append(path)
    return unique


def ensure_binary(binary: str) -> None:
    if not Path(binary).is_file():
        sys.stderr.write(f"[run_bench] error: solver binary '{binary}' not found.\n")
        sys.exit(1)
    if not os.access(binary, os.X_OK):
        sys.stderr.write(f"[run_bench] error: solver binary '{binary}' is not executable.\n")
        sys.exit(1)


def invoke_solver(binary: str, input_file: str, time_limit: Optional[float]) -> dict:
    """Run the solver and parse its CSV line into typed fields."""
    cmd = [binary, input_file]
    if time_limit is not None:
        cmd.append(str(time_limit))
    result = subprocess.run(
        cmd,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write(
            f"[run_bench] solver failed on {input_file} (exit {result.returncode}).\n"
        )
        if result.stderr:
            sys.stderr.write(result.stderr + "\n")
        sys.exit(result.returncode)

    stdout = result.stdout.strip()
    if not stdout:
        sys.stderr.write(
            f"[run_bench] solver produced no output for {input_file}.\n"
        )
        sys.exit(1)

    # Take the last non-empty line in case the binary prints extra info.
    for line in reversed(stdout.splitlines()):
        striped = line.strip()
        if striped:
            fields = [chunk.strip() for chunk in striped.split(",")]
            break
    else:  # no break
        sys.stderr.write(f"[run_bench] unable to parse solver output for {input_file}.\n")
        sys.exit(1)

    if len(fields) not in (9, 11):
        sys.stderr.write(
            f"[run_bench] unexpected column count ({len(fields)}) in line: '{striped}'.\n"
        )
        sys.exit(1)

    try:
        sort_ms = 0.0
        stock_ms = 0.0
        status = fields[8]
        if len(fields) == 11:
            sort_ms = float(fields[9])
            stock_ms = float(fields[10])
        return {
            "solver_input": fields[0],
            "num_vertices": int(fields[1]),
            "num_components": int(fields[2]),
            "min_component_size": int(fields[3]),
            "max_component_size": int(fields[4]),
            "num_trees": int(fields[5]),
            "elapsed_ms": float(fields[6]),
            "delay_ms": float(fields[7]),
            "elapsed_sort_ms": sort_ms,
            "elapsed_stock_ms": stock_ms,
            "status": status,
        }
    except ValueError as err:
        sys.stderr.write(
            f"[run_bench] failed to parse numeric fields from line '{striped}': {err}.\n"
        )
        sys.exit(1)


def run_repetitions(
    binary: str,
    input_path: str,
    label: str,
    reps: int,
    warmup: int,
    quiet: bool,
    time_limit: Optional[float],
) -> List[dict]:
    """Run warmups, then timed repetitions, returning parsed rows."""
    for _ in range(warmup):
        invoke_solver(binary, input_path, time_limit)
    rows: List[dict] = []
    for rep in range(1, reps + 1):
        stats = invoke_solver(binary, input_path, time_limit)
        row = {
            "label": label,
            "rep": rep,
            **stats,
        }
        rows.append(row)
        if not quiet:
            sys.stderr.write(
                f"[run_bench] {label} rep {rep}: {stats['elapsed_ms']:.3f} ms\n"
            )
    return rows


def append_rows(rows: Iterable[dict], output_path: str) -> None:
    if not rows:
        return
    write_header = not Path(output_path).exists()
    with open(output_path, "a", encoding="ascii") as handle:
        if write_header:
            handle.write(CSV_HEADER + "\n")
        for row in rows:
            line = (
                f"{row['label']},{row['rep']},"
                f"{row['solver_input']},{row['num_vertices']},{row['num_components']},"
                f"{row['min_component_size']},{row['max_component_size']},"
                f"{row['num_trees']},{row['elapsed_ms']:.3f},{row['delay_ms']},"
                f"{row['elapsed_sort_ms']:.6f},{row['elapsed_stock_ms']:.6f},{row['status']}"
            )
            handle.write(line + "\n")


def main() -> None:
    args = parse_args()
    ensure_binary(BINARY_PATH)

    all_rows: List[dict] = []

    input_files = discover_inputs(args.inputs)
    if input_files:
        for path in input_files:
            label = Path(path).name
            rows = run_repetitions(
                BINARY_PATH,
                path,
                label,
                reps=args.reps,
                warmup=args.warmup,
                quiet=args.quiet,
                time_limit=args.time_limit,
            )
            all_rows.extend(rows)
    elif not args.quiet:
        sys.stderr.write("[run_bench] no deterministic inputs selected.\n")

    append_rows(all_rows, args.output)
    if not args.quiet:
        sys.stderr.write(
            f"[run_bench] wrote {len(all_rows)} measurement(s) to {args.output}.\n"
        )


if __name__ == "__main__":
    main()
