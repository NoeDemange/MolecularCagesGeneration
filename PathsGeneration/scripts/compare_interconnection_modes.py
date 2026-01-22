#!/usr/bin/env python3
"""Compare interconnection tree execution modes on demos_realSubstrat datasets.

The script runs the CagePathGen binary on every requested substrate while toggling:

* The interconnection-tree strategy (`-g 0` on-the-fly vs `-g 1` store & sort)
* The banned-edges option (`-b 0/1`, configurable)

It prints a compact table and can optionally emit a CSV summary. By default,
`-t` is set to 0 so multiple cages per interconnection tree are permitted,
which makes it easier to observe the effect of banned edges.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATTERN = re.compile(r"RESULT_SUMMARY\s+(.*)")
TREES_PATTERN = re.compile(r"Interconnection Trees:\s+(\d+)")
COMPONENT_PATTERN = re.compile(r"COMPONENT_SUMMARY\s+(.*)")


@dataclass
class ModeConfig:
    name: str
    value: int


@dataclass
class DatasetTarget:
    path: Path
    moc_override: Optional[int]


MODES: Dict[str, ModeConfig] = {
    "on-the-fly": ModeConfig("on-the-fly", 0),
    "store-sort": ModeConfig("store-sort", 1),
}


CSV_FIELDS = [
    "dataset",
    "moc",
    "mode",
    "banned",
    "status",
    "return_code",
    "results",
    "best_rmsd",
    "worst_rmsd",
    "avg_path_length",
    "max_path_length",
    "time_first_result_ms",
    "total_time_ms",
    "interconnection_trees",
    "components",
    "min_component_vertices",
    "max_component_vertices",
    "min_atom_count",
    "max_atom_count",
    "results_dir",
]


def run_make(command: Sequence[str]) -> None:
    result = subprocess.run(command, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)


def build_project(extra_flags: str) -> None:
    run_make(["make", "clean"])
    run_make(["make", f"EXTRA_CFLAGS={extra_flags}"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--demo-root",
        type=Path,
        default=Path("demos_realSubstrat"),
        help="Directory containing the real substrate demos",
    )
    parser.add_argument(
        "--names",
        nargs="*",
        help="Optional subset of demo directories to run",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("./bin/chemins.exe"),
        help="Path to the compiled CagePathGen binary",
    )
    parser.add_argument(
        "--size-max",
        type=int,
        default=5,
        help="Value for -s (patterns per path)",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=20,
        help="Value for -r (results per dataset)",
    )
    parser.add_argument(
        "--one-cage",
        type=int,
        default=0,
        help="Value for -t (one cage per tree). Default 0 so multiple cages per tree are allowed.",
    )
    parser.add_argument(
        "--banned-values",
        type=int,
        nargs="+",
        default=[0, 1],
        help="Values to test for -b (banned edges toggle)",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=list(MODES.keys()),
        default=list(MODES.keys()),
        help="Interconnection-tree modes to compare",
    )
    parser.add_argument(
        "--default-moc",
        type=int,
        default=0,
        help="Valeur utilisée pour -n lorsqu'aucun fichier moc n'est trouvé (utile pour les entrées mono-fichier)",
    )
    parser.add_argument(
        "--distance-type",
        choices=["A*", "SSMTA*", "HYBRID", "Euclidienne"],
        help="Override DISTANCE_TYPE before running",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Timeout per run in seconds",
    )
    parser.add_argument(
        "--timeout-grace",
        type=int,
        default=15,
        help="Additional seconds to wait after sending SIGINT so the solver can flush stats",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path("results_compare_modes"),
        help="Répertoire racine pour stocker les sorties lorsqu'on conserve les résultats",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Fichier CSV optionnel pour ajouter les résumés de runs",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path("results_compare_modes/logs"),
        help="Répertoire dans lequel sauvegarder un log (commande, stdout, stderr) pour chaque exécution",
    )
    parser.add_argument(
        "--keep-results",
        action="store_true",
        help="Conserver les fichiers générés (.mol2, paramètres, temps). Si absent, le script active CAGE_DISABLE_WRITE",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Ne pas recompiler automatiquement le binaire avec -DENABLE_STATS avant les runs",
    )
    parser.add_argument(
        "--extra-cflags",
        default="",
        help="Flags additionnels passés à EXTRA_CFLAGS lors de la recompilation automatique",
    )
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments forwarded to the binary",
    )
    return parser.parse_args()


def _parse_dataset_spec(entry: str) -> tuple[str, Optional[int]]:
    if ":" not in entry:
        return entry, None
    name, moc_str = entry.split(":", 1)
    moc_str = moc_str.strip()
    if not moc_str:
        return name, None
    try:
        moc_value = int(moc_str)
    except ValueError:
        raise ValueError(f"Invalid moc number '{moc_str}' in '{entry}'") from None
    return name, moc_value


def discover_datasets(root: Path, wanted: Optional[Iterable[str]]) -> List[DatasetTarget]:
    targets: List[DatasetTarget] = []

    def _resolve_path(spec: str) -> Path:
        candidate = root / spec
        if candidate.exists():
            return candidate
        alt = Path(spec)
        if alt.exists():
            return alt
        raise FileNotFoundError(f"Dataset '{spec}' not found in {root} or as absolute path")

    if wanted:
        for raw_entry in wanted:
            entry = raw_entry.strip()
            if not entry:
                continue
            name, moc_override = _parse_dataset_spec(entry)
            path = _resolve_path(name)
            targets.append(DatasetTarget(path=path, moc_override=moc_override))
        return targets

    if root.is_file():
        targets.append(DatasetTarget(path=root, moc_override=None))
        return targets

    if not root.is_dir():
        raise FileNotFoundError(f"Demo root {root} does not exist")

    for candidate in sorted(p for p in root.iterdir() if p.is_dir()):
        targets.append(DatasetTarget(path=candidate, moc_override=None))
    return targets


def guess_moc_number(dataset: Path, moc_override: Optional[int], default_moc: int) -> int:
    if moc_override is not None:
        return moc_override
    if dataset.is_file():
        return default_moc
    pattern = re.compile(r"_moc(\d+)\.mol2$")
    numbers: List[int] = []
    for mol2 in dataset.glob(f"{dataset.name}_moc*.mol2"):
        match = pattern.search(mol2.name)
        if match:
            numbers.append(int(match.group(1)))
    if not numbers:
        return default_moc
    return min(numbers)


def parse_summary(payload: str) -> Dict[str, str]:
    for line in payload.splitlines():
        match = SUMMARY_PATTERN.search(line)
        if not match:
            continue
        tokens: Dict[str, str] = {}
        for chunk in match.group(1).split():
            if "=" not in chunk:
                continue
            key, value = chunk.split("=", 1)
            tokens[key] = value
        return tokens
    return {}


def parse_component_summary(payload: str) -> Dict[str, str]:
    for line in payload.splitlines():
        match = COMPONENT_PATTERN.search(line)
        if not match:
            continue
        tokens: Dict[str, str] = {}
        for chunk in match.group(1).split():
            if "=" not in chunk:
                continue
            key, value = chunk.split("=", 1)
            tokens[key] = value
        return tokens
    return {}


def extract_tree_count(payload: str) -> Optional[int]:
    match = TREES_PATTERN.search(payload)
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return None
    return None


def record_log(log_path: Path, command: Sequence[str], status: str, return_code: Optional[int], stdout_text: str, stderr_text: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    def _ensure_text(value: object) -> str:
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        if value is None:
            return ""
        return str(value)

    stdout_str = _ensure_text(stdout_text)
    stderr_str = _ensure_text(stderr_text)

    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Command: {' '.join(command)}\n")
        handle.write(f"Status: {status}\n")
        handle.write(f"Return code: {'' if return_code is None else return_code}\n")
        handle.write("\n[stdout]\n")
        handle.write(stdout_str if stdout_str.endswith("\n") else stdout_str + "\n")
        handle.write("\n[stderr]\n")
        handle.write(stderr_str if stderr_str.endswith("\n") else stderr_str + "\n")


def run_single(
    dataset: Path,
    moc_number: int,
    args: argparse.Namespace,
    mode: ModeConfig,
    banned_value: int,
) -> Dict[str, object]:
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    result_dir: Optional[Path] = None
    if args.keep_results:
        result_dir = args.results_root / dataset.name / f"mode_{mode.value}_b{banned_value}_{timestamp}"
        result_dir.mkdir(parents=True, exist_ok=True)

    log_dir = args.log_dir / dataset.name
    log_path = log_dir / f"mode_{mode.value}_b{banned_value}_{timestamp}.log"

    env = os.environ.copy()
    if args.keep_results and result_dir is not None:
        env["CAGE_RESULTS_DIR"] = str(result_dir.resolve())
    else:
        env["CAGE_DISABLE_WRITE"] = "1"
    if args.distance_type:
        env["DISTANCE_TYPE"] = args.distance_type

    command: List[str] = [
        str(args.binary),
        "-i",
        str(dataset.resolve()),
        "-n",
        str(moc_number),
        "-s",
        str(args.size_max),
        "-r",
        str(args.max_results),
        "-b",
        str(banned_value),
        "-t",
        str(args.one_cage),
        "-g",
        str(mode.value),
    ]
    if args.extra_args:
        command.extend(args.extra_args)

    status = "ok"
    return_code: Optional[int] = None
    stdout_text = ""
    stderr_text = ""
    start = time.perf_counter()

    def _ensure_text_local(value: object) -> str:
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        if value is None:
            return ""
        return str(value)

    process = subprocess.Popen(
        command,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    try:
        stdout_text, stderr_text = process.communicate(timeout=args.timeout)
        return_code = process.returncode
        if return_code != 0:
            status = "error"
    except subprocess.TimeoutExpired as exc:
        status = "timeout"
        stdout_chunks = [_ensure_text_local(exc.stdout)]
        stderr_chunks = [_ensure_text_local(exc.stderr), f"\nTimed out after {args.timeout}s; sending SIGINT for stats..."]
        try:
            process.send_signal(signal.SIGINT)
        except Exception as send_error:  # pragma: no cover - best-effort signal delivery
            stderr_chunks.append(f"\nFailed to send SIGINT: {send_error}")
        try:
            extra_stdout, extra_stderr = process.communicate(timeout=args.timeout_grace)
            stdout_chunks.append(_ensure_text_local(extra_stdout))
            stderr_chunks.append(_ensure_text_local(extra_stderr))
        except subprocess.TimeoutExpired:
            stderr_chunks.append(
                f"\nNo shutdown after {args.timeout_grace}s grace; sending SIGKILL."
            )
            process.kill()
            kill_stdout, kill_stderr = process.communicate()
            stdout_chunks.append(_ensure_text_local(kill_stdout))
            stderr_chunks.append(_ensure_text_local(kill_stderr))
        stdout_text = "".join(stdout_chunks)
        stderr_text = "".join(stderr_chunks)
        return_code = process.returncode

    stdout_text = _ensure_text_local(stdout_text)
    stderr_text = _ensure_text_local(stderr_text)

    total_time = time.perf_counter() - start
    record_log(log_path, command, status, return_code, stdout_text, stderr_text)
    merged_output = stdout_text + "\n" + stderr_text
    summary = parse_summary(merged_output)
    component_summary = parse_component_summary(merged_output)
    tree_count = extract_tree_count(merged_output)

    def _safe_float(*keys: str) -> Optional[float]:
        for key in keys:
            if not key:
                continue
            value = summary.get(key)
            if value is None or value.upper() == "NA":
                continue
            try:
                return float(value)
            except ValueError:
                continue
        return None

    def _safe_int_value(raw: Optional[str]) -> Optional[int]:
        if raw is None:
            return None
        if raw.upper() == "NA":
            return None
        try:
            return int(raw)
        except ValueError:
            return None

    def _safe_time_ms(key_ms: str, key_s: str) -> Optional[float]:
        value_ms = summary.get(key_ms)
        if value_ms is not None and value_ms.upper() != "NA":
            try:
                return float(value_ms)
            except ValueError:
                pass
        value_s = summary.get(key_s)
        if value_s is not None and value_s.upper() != "NA":
            try:
                return float(value_s) * 1000.0
            except ValueError:
                return None
        return None

    def _safe_int(*keys: str) -> Optional[int]:
        for key in keys:
            if not key:
                continue
            value = summary.get(key)
            parsed = _safe_int_value(value)
            if parsed is not None:
                return parsed
        return None

    def _component_int(key: str) -> Optional[int]:
        return _safe_int_value(component_summary.get(key))

    return {
        "dataset": dataset.name,
        "moc": moc_number,
        "mode": mode.name,
        "banned": banned_value,
        "status": status,
        "return_code": return_code,
        "results": _safe_int("results"),
        "best_rmsd": _safe_float("best_rmsd"),
        "worst_rmsd": _safe_float("worst_rmsd"),
        "avg_path_length": _safe_float("avg_path"),
        "max_path_length": _safe_int("max_path"),
        "time_first_result_ms": _safe_time_ms("first_result_ms", "first_result_s"),
        "total_time_ms": _safe_time_ms("total_time_ms", "total_time_s") or round(total_time * 1000.0, 3),
        "interconnection_trees": tree_count,
        "components": _component_int("components"),
        "min_component_vertices": _component_int("min_vertices"),
        "max_component_vertices": _component_int("max_vertices"),
        "min_atom_count": _safe_int("min_atoms"),
        "max_atom_count": _safe_int("max_atoms"),
        "results_dir": str(result_dir.resolve()) if result_dir else "(stats-only)",
    }


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows)


def print_table(rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    headers = [
        "dataset",
        "mode",
        "banned",
        "status",
        "results",
        "avg_path_length",
        "max_path_length",
        "interconnection_trees",
        "components",
        "min_component_vertices",
        "max_component_vertices",
        "min_atom_count",
        "max_atom_count",
        "time_first_result_ms",
        "total_time_ms",
    ]
    widths = {key: max(len(str(key)), max(len(str(row.get(key, ""))) for row in rows)) for key in headers}
    line = " | ".join(f"{key.upper():<{widths[key]}}" for key in headers)
    print(line)
    print("-" * len(line))
    for row in rows:
        print(
            " | ".join(
                f"{str(row.get(key, '')):<{widths[key]}}" for key in headers
            )
        )


def main() -> int:
    args = parse_args()
    binary_path = args.binary if args.binary.is_absolute() else (REPO_ROOT / args.binary)
    if not args.skip_build:
        extra_flags = "-DENABLE_STATS"
        if args.extra_cflags:
            extra_flags = f"{extra_flags} {args.extra_cflags}"
        print(f"Rebuilding project with EXTRA_CFLAGS=\"{extra_flags}\"")
        try:
            build_project(extra_flags)
        except subprocess.CalledProcessError as exc:
            print(f"Build failed (code {exc.returncode}). Check the compiler output above.", file=sys.stderr)
            return 1
    args.binary = binary_path.resolve()
    if not args.binary.is_file():
        print(f"Binary {args.binary} not found. Run 'make' first.", file=sys.stderr)
        return 1

    try:
        datasets = discover_datasets(args.demo_root, args.names)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    if args.keep_results:
        args.results_root.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, object]] = []
    for target in datasets:
        dataset = target.path
        print(f"\n>>> Dataset {dataset.name}")
        try:
            moc = guess_moc_number(dataset, target.moc_override, args.default_moc)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            continue
        for mode_name in args.modes:
            mode = MODES[mode_name]
            for banned in args.banned_values:
                print(f"  - mode={mode.name} (-g {mode.value}), banned={banned}")
                row = run_single(dataset, moc, args, mode, banned)
                all_rows.append(row)
                print(
                    f"    status={row['status']} return={row['return_code']} results={row['results']}\n"
                    f"    time={row['total_time_ms']}ms output={row['results_dir']}"
                )
                print(
                    "    "
                    f"trees={row['interconnection_trees']} components={row['components']} "
                    f"min_comp={row['min_component_vertices']} max_comp={row['max_component_vertices']} "
                    f"min_atoms={row['min_atom_count']} max_atoms={row['max_atom_count']}"
                )

    if args.summary_csv:
        write_csv(args.summary_csv, all_rows)

    print("\n=== Comparison Summary ===")
    print_table(all_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())