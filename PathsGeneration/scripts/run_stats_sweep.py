#!/usr/bin/env python3
"""Sweep NUMBER_POSITION_AX1E3 (branching positions) / MAX_POSITION_KEEP_360 (max discretized position budget) / THRESHOLD_ANGLE using stats-only runs."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATTERN = re.compile(r"RESULT_SUMMARY\s+(.*)")
TREES_PATTERN = re.compile(r"Interconnection Trees:\s+(\d+)")
STATS_PATTERN = re.compile(r"STATS_SUMMARY\s+(.*)")
MAX_POSITION_WARNING = re.compile(r"Warning: MAX_POSITION_KEEP_360")

HEADER = [
    "dataset",
    "moc",
    "size_max",
    "number_branching_positions",
    "max_discretized_position",
    "threshold_angle_deg",
    "discretization_type",
    "step_grid_voxel",
    "distance_type",
    "path_boundary",
    "dynamic_path_limit",
    "status",
    "return_code",
    "results_count",
    "best_rmsd",
    "worst_rmsd",
    "min_path_length",
    "min_path_count",
    "min_path_best_rmsd",
    "max_path_length",
    "avg_path_length",
    "first_result_ms",
    "total_time_ms",
    "cpu_time_ms",
    "interconnection_trees",
    "branches",
    "without_intersections",
    "with_intersections",
    "max_position_warning",
    "mean_intervals_with",
    "mean_intervals_overall",
    "max_intervals",
    "mean_coverage_pct",
    "max_coverage_pct",
    "total_collisions",
    "circle_blocked",
    "boundary_allowed",
    "boundary_blocked",
    "build_flags",
    "log_file",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo-root", type=Path, default=Path("demos_onePath"), help="Folder with demo inputs")
    parser.add_argument("--names", nargs="*", help="Subset of demo names to keep")
    parser.add_argument("--binary", type=Path, default=Path("./bin/chemins.exe"), help="Compiled binary path")
    parser.add_argument("--size-min", type=int, default=5, help="Minimum value for -s")
    parser.add_argument("--size-max", type=int, default=15, help="Maximum value for -s (inclusive)")
    parser.add_argument(
        "--number-branching-positions",
        type=int,
        nargs="+",
        default=[3],
        dest="number_branching_positions",
        help="Values for NUMBER_POSITION_AX1E3 (default: 3)",
    )
    parser.add_argument(
        "--max-discretized-position",
        type=int,
        nargs="+",
        default=[12],
        dest="max_discretized_position",
        help="Values for MAX_POSITION_KEEP_360 (max discretized position budget, default: 12)",
    )
    parser.add_argument(
        "--threshold-angles",
        type=float,
        nargs="+",
        default=[15],
        help="Values (degrees) for THRESHOLD_ANGLE",
    )
    parser.add_argument(
        "--discretization-types",
        type=int,
        nargs="+",
        default=[1],
        dest="discretization_types",
        help="Values for DISCRITIZATION_TYPE (compile-time discretization mode)",
    )
    parser.add_argument(
        "--step-grid-voxels",
        type=float,
        nargs="+",
        default=[0.5],
        dest="step_grid_voxels",
        help="Values for STEP_GRID_VOXEL (grid spacing in angstroms)",
    )
    parser.add_argument(
        "--distance-types",
        nargs="+",
        default=["Euclidienne"],
        help="Distance backends to test (ENV DISTANCE_TYPE)",
    )
    parser.add_argument(
        "--path-boundary",
        choices=["auto", "on", "off"],
        default="auto",
        help="Force DIST_PATH_BOUNDARY on/off (via CAGE_PATH_BOUNDARY); 'auto' leaves the environment untouched",
    )
    parser.add_argument(
        "--dynamic-path-limit",
        choices=["auto", "on", "off"],
        nargs="+",
        default=["auto"],
        dest="dynamic_path_limit",
        help="Control the -l flag: 'on' passes -l 1, 'off' sends -l 0; repeat the option with both values to compare",
    )
    parser.add_argument("--max-results", type=int, default=1, help="Value for -r")
    parser.add_argument("--banned-edges", type=int, default=0, help="Value for -b (fixed to 0 by default)")
    parser.add_argument("--one-cage", type=int, default=0, help="Value for -t (fixed to 0 by default)")
    parser.add_argument("--timeout", type=int, default=120, help="Timeout per binary run (seconds)")
    parser.add_argument("--summary-csv", type=Path, default=Path("benchmark_runs/stats_sweep.csv"))
    parser.add_argument("--log-dir", type=Path, default=Path("benchmark_runs/logs"))
    parser.add_argument("--extra-cflags", default="", help="Additional flags appended to EXTRA_CFLAGS")
    parser.add_argument("--extra-args", nargs=argparse.REMAINDER, help="Extra CLI args forwarded to the binary")
    return parser.parse_args()


def discover_datasets(root: Path, selection: Optional[Iterable[str]]) -> List[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"Demo root {root} does not exist")
    candidates = sorted(p for p in root.iterdir() if p.is_dir())
    if selection:
        wanted = {name.strip() for name in selection}
        candidates = [p for p in candidates if p.name in wanted]
    return candidates


def guess_moc_number(dataset: Path) -> int:
    pattern = re.compile(r"_moc(\d+)\.mol2$")
    candidates: List[int] = []
    for mol in dataset.glob(f"{dataset.name}_moc*.mol2"):
        match = pattern.search(mol.name)
        if match:
            candidates.append(int(match.group(1)))
    if not candidates:
        raise FileNotFoundError(f"No moc file found in {dataset}")
    return min(candidates)


def run_make(command: Sequence[str]) -> None:
    result = subprocess.run(command, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)


def build_project(extra_flags: str) -> None:
    run_make(["make", "clean"])
    run_make(["make", f"EXTRA_CFLAGS={extra_flags}"])


def _parse_summary(text: str) -> Optional[Dict[str, str]]:
    tokens: Dict[str, str] = {}
    for line in text.splitlines():
        match = SUMMARY_PATTERN.search(line)
        if match:
            for chunk in match.group(1).split():
                if "=" not in chunk:
                    continue
                key, value = chunk.split("=", 1)
                tokens[key] = value
    return tokens or None


def _parse_stats_summary(text: str) -> Optional[Dict[str, str]]:
    for line in text.splitlines():
        match = STATS_PATTERN.search(line)
        if match:
            tokens: Dict[str, str] = {}
            for chunk in match.group(1).split():
                if "=" not in chunk:
                    continue
                key, value = chunk.split("=", 1)
                tokens[key] = value
            return tokens
    return None


def _parse_float(value: Optional[str]) -> Optional[float]:
    if value is None or value.upper() == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _extract_time_milliseconds(summary: Dict[str, str], key_ms: str, key_s: str) -> Optional[float]:
    """Return the time value in milliseconds, preferring explicit ms keys."""
    value_ms = summary.get(key_ms)
    if value_ms is not None:
        parsed_ms = _parse_float(value_ms)
        if parsed_ms is not None:
            return parsed_ms
    value_s = summary.get(key_s)
    parsed_s = _parse_float(value_s)
    if parsed_s is not None:
        return parsed_s * 1000.0
    return None


def _parse_int(value: Optional[str]) -> Optional[int]:
    if value is None or value.upper() == "NA":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def append_row(csv_path: Path, row: Dict[str, object]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = csv_path.exists()
    with csv_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def format_float(value: Optional[float], precision: str = ".6f") -> str:
    if value is None:
        return ""
    return f"{value:{precision}}"


def _ensure_text(data: Optional[object]) -> str:
    if data is None:
        return ""
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return str(data)


def _normalize_explicit_boundary(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    normalized = value.strip().lower()
    if not normalized or normalized == "auto":
        return None
    if normalized in {"0", "false", "f", "no", "n", "off"}:
        return "off"
    return "on"


def _boundary_state_from_env(env: Dict[str, str]) -> str:
    explicit = _normalize_explicit_boundary(env.get("CAGE_PATH_BOUNDARY"))
    if explicit:
        return explicit
    legacy = env.get("CAGE_DISABLE_PATH_BOUNDARY")
    if legacy is not None:
        legacy_norm = legacy.strip().lower()
        if legacy_norm and legacy_norm not in {"0", "false", "f", "no", "n", "off"}:
            return "off"
    return "auto"


def record_log(log_path: Path, command: Sequence[str], status: str, stdout: object, stderr: object) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    stdout_text = _ensure_text(stdout) or "<empty>\n"
    stderr_text = _ensure_text(stderr) or "<empty>\n"
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Command: {' '.join(command)}\n")
        handle.write(f"Status: {status}\n")
        handle.write("\n[stdout]\n")
        handle.write(stdout_text if stdout_text.endswith("\n") else stdout_text + "\n")
        handle.write("\n[stderr]\n")
        handle.write(stderr_text if stderr_text.endswith("\n") else stderr_text + "\n")


def run_case(
    binary: Path,
    dataset: Path,
    moc_number: int,
    size_max: int,
    distance_type: str,
    dynamic_path_limit: str,
    args: argparse.Namespace,
    build_flags: str,
    log_dir: Path,
) -> Dict[str, object]:
    env = os.environ.copy()
    env["CAGE_DISABLE_WRITE"] = "1"
    env["DISTANCE_TYPE"] = distance_type
    if args.path_boundary == "on":
        env["CAGE_PATH_BOUNDARY"] = "1"
        env.pop("CAGE_DISABLE_PATH_BOUNDARY", None)
    elif args.path_boundary == "off":
        env["CAGE_PATH_BOUNDARY"] = "0"
        env.pop("CAGE_DISABLE_PATH_BOUNDARY", None)
    path_boundary_state = _boundary_state_from_env(env)

    command = [
        str(binary),
        "-i",
        str(dataset.resolve()),
        "-n",
        str(moc_number),
        "-s",
        str(size_max),
        "-r",
        str(args.max_results),
        "-b",
        str(args.banned_edges),
        "-t",
        str(args.one_cage),
    ]
    if dynamic_path_limit == "on":
        command.extend(["-l", "1"])
    elif dynamic_path_limit == "off":
        command.extend(["-l", "0"])
    if args.extra_args:
        command.extend(args.extra_args)

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    safe_dist = re.sub(r"[^A-Za-z0-9._-]+", "_", distance_type)
    dyn_tag = dynamic_path_limit.capitalize()
    log_path = log_dir / f"{dataset.name}_s{size_max}_dist{safe_dist}_dyn{dyn_tag}_{timestamp}.log"

    status = "ok"
    return_code: Optional[int] = None
    stdout_text = ""
    stderr_text = ""

    def _merge_outputs(out: object, err: object) -> str:
        return _ensure_text(out) + "\n" + _ensure_text(err)

    try:
        completed = subprocess.run(
            command,
            env=env,
            capture_output=True,
            text=True,
            timeout=args.timeout,
            check=False,
        )
        stdout_text = completed.stdout
        stderr_text = completed.stderr
        return_code = completed.returncode
        if completed.returncode != 0:
            status = "error"
    except subprocess.TimeoutExpired as exc:
        stdout_text = exc.stdout or b""
        stderr_text = (exc.stderr or b"") + f"\nTimed out after {args.timeout}s".encode("utf-8")
        status = "timeout"

    record_log(log_path, command, status, stdout_text, stderr_text)

    merged_output = _merge_outputs(stdout_text, stderr_text)
    warning_present = bool(MAX_POSITION_WARNING.search(merged_output))

    summary = _parse_summary(merged_output)
    stats_summary = _parse_stats_summary(merged_output)
    trees = None
    tree_match = TREES_PATTERN.search(merged_output)
    if tree_match:
        trees = int(tree_match.group(1))

    if summary is None:
        if status == "ok":
            status = "missing-summary"
        results = None
        best = worst = min_path = max_path = avg_path = first = total = cpu = None
        min_path_count = None
        min_path_best = None
    else:
        results = _parse_int(summary.get("results"))
        best = _parse_float(summary.get("best_rmsd"))
        worst = _parse_float(summary.get("worst_rmsd"))
        min_path = _parse_int(summary.get("min_path"))
        min_path_count = _parse_int(summary.get("min_path_count"))
        min_path_best = _parse_float(summary.get("min_path_best_rmsd"))
        min_path_count = _parse_int(summary.get("min_path_count"))
        min_path_best = _parse_float(summary.get("min_path_best_rmsd"))
        max_path = _parse_int(summary.get("max_path"))
        avg_path = _parse_float(summary.get("avg_path"))
        first = _extract_time_milliseconds(summary, "first_result_ms", "first_result_s")
        total = _extract_time_milliseconds(summary, "total_time_ms", "total_time_s")
        cpu = _extract_time_milliseconds(summary, "cpu_time_ms", "cpu_time_s")

    if stats_summary is None:
        branches = without_inter = with_inter = max_inter = None
        mean_with = mean_all = mean_cov = max_cov = None
        collisions = blocked = None
        boundary_allowed = boundary_blocked = None
    else:
        branches = _parse_int(stats_summary.get("branches"))
        without_inter = _parse_int(stats_summary.get("without_intersections"))
        with_inter = _parse_int(stats_summary.get("with_intersections"))
        mean_with = _parse_float(stats_summary.get("mean_intervals_with"))
        mean_all = _parse_float(stats_summary.get("mean_intervals_overall"))
        max_inter = _parse_int(stats_summary.get("max_intervals"))
        mean_cov = _parse_float(stats_summary.get("mean_coverage_pct"))
        max_cov = _parse_float(stats_summary.get("max_coverage_pct"))
        collisions = _parse_int(stats_summary.get("total_collisions"))
        blocked = _parse_int(stats_summary.get("circle_blocked"))
        boundary_allowed = _parse_int(stats_summary.get("boundary_allowed"))
        boundary_blocked = _parse_int(stats_summary.get("boundary_blocked"))

    row = {
        "dataset": dataset.name,
        "moc": moc_number,
        "size_max": size_max,
        "number_branching_positions": "",
        "max_discretized_position": "",
        "threshold_angle_deg": "",
        "step_grid_voxel": "",
        "distance_type": distance_type,
        "path_boundary": path_boundary_state,
        "dynamic_path_limit": dynamic_path_limit,
        "status": status,
        "return_code": "" if return_code is None else return_code,
        "results_count": "" if results is None else results,
        "best_rmsd": format_float(best),
        "worst_rmsd": format_float(worst),
        "min_path_length": "" if min_path is None else min_path,
        "min_path_count": "" if min_path_count is None else min_path_count,
        "min_path_best_rmsd": format_float(min_path_best),
        "max_path_length": "" if max_path is None else max_path,
        "avg_path_length": format_float(avg_path, ".2f"),
        "first_result_ms": format_float(first, ".3f"),
        "total_time_ms": format_float(total, ".3f"),
        "cpu_time_ms": format_float(cpu, ".3f"),
        "interconnection_trees": "" if trees is None else trees,
        "branches": "" if branches is None else branches,
        "without_intersections": "" if without_inter is None else without_inter,
        "with_intersections": "" if with_inter is None else with_inter,
        "max_position_warning": "1" if warning_present else "0",
        "mean_intervals_with": format_float(mean_with, ".2f"),
        "mean_intervals_overall": format_float(mean_all, ".2f"),
        "max_intervals": "" if max_inter is None else max_inter,
        "mean_coverage_pct": format_float(mean_cov, ".2f"),
        "max_coverage_pct": format_float(max_cov, ".2f"),
        "total_collisions": "" if collisions is None else collisions,
        "circle_blocked": "" if blocked is None else blocked,
        "boundary_allowed": "" if boundary_allowed is None else boundary_allowed,
        "boundary_blocked": "" if boundary_blocked is None else boundary_blocked,
        "build_flags": build_flags.strip(),
        "log_file": str(log_path.resolve()),
    }
    return row


def main() -> int:
    args = parse_args()
    demo_root = args.demo_root if args.demo_root.is_absolute() else (REPO_ROOT / args.demo_root)
    binary_path = (args.binary if args.binary.is_absolute() else (REPO_ROOT / args.binary)).resolve()
    log_dir = args.log_dir if args.log_dir.is_absolute() else (REPO_ROOT / args.log_dir)
    csv_path = args.summary_csv if args.summary_csv.is_absolute() else (REPO_ROOT / args.summary_csv)
    log_dir.mkdir(parents=True, exist_ok=True)

    datasets = discover_datasets(demo_root, args.names)
    if not datasets:
        print("No demo directories found", file=sys.stderr)
        return 1

    size_values = list(range(args.size_min, args.size_max + 1))
    combos = []
    for num_pos in args.number_branching_positions:
        for max_keep in args.max_discretized_position:
            if max_keep <= 0:
                raise ValueError("Values from --max-discretized-position must be positive")
            for threshold in args.threshold_angles:
                threshold_rad = math.radians(threshold)
                for disc_type in args.discretization_types:
                    for step_voxel in args.step_grid_voxels:
                        if step_voxel <= 0:
                            raise ValueError("Values from --step-grid-voxels must be positive")
                        base_flags = (
                            f"-DENABLE_STATS -DNUMBER_POSITION_AX1E3={num_pos} -DNUMBER_POSITION_PATHS={num_pos} "
                            f"-DMAX_POSITION_KEEP_360={max_keep} -DTHRESHOLD_ANGLE={threshold_rad:.12f} "
                            f"-DDISCRITIZATION_TYPE={disc_type} -DSTEP_GRID_VOXEL={step_voxel:.12f}"
                        )
                        if args.extra_cflags:
                            base_flags = f"{base_flags} {args.extra_cflags}"
                        combos.append((num_pos, max_keep, threshold, disc_type, step_voxel, base_flags))

    had_error = False

    for num_pos, max_keep, threshold, disc_type, step_voxel, flags in combos:
        print(
            f"\n=== Build NBP={num_pos} MAX_DISCRETIZED={max_keep} THR={threshold} STEP={step_voxel} ==="
        )
        try:
            build_project(flags)
        except subprocess.CalledProcessError as exc:
            print(f"Build failed for flags '{flags}': return code {exc.returncode}", file=sys.stderr)
            had_error = True
            continue

        for dynamic_state in args.dynamic_path_limit:
            for distance in args.distance_types:
                for dataset in datasets:
                    try:
                        moc_number = guess_moc_number(dataset)
                    except FileNotFoundError as exc:
                        print(str(exc), file=sys.stderr)
                        had_error = True
                        continue
                    for size_max in size_values:
                        row = run_case(
                            binary_path,
                            dataset,
                            moc_number,
                            size_max,
                            distance,
                            dynamic_state,
                            args,
                            flags,
                            log_dir,
                        )
                        row["number_branching_positions"] = num_pos
                        row["max_discretized_position"] = max_keep
                        row["threshold_angle_deg"] = threshold
                        row["discretization_type"] = disc_type
                        row["step_grid_voxel"] = step_voxel
                        append_row(csv_path, row)
                        print(
                            f"{dataset.name} size={size_max} dist={distance} dyn={dynamic_state} disc={disc_type} step={step_voxel} status={row['status']} results={row['results_count']}"
                        )
                        if row["status"] != "ok":
                            had_error = True

    return 1 if had_error else 0


if __name__ == "__main__":
    raise SystemExit(main())
