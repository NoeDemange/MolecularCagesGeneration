<h1 align="center"> 💻 Molecular Cage Guides ⚛️</h1>
<h3 align="center"> Paths Generation</h3>
<p>
</p>


This module produces the second half of the substrate-specific molecular cage construction guides. Starting from an initial partial cage (interaction points on the substrate), it explores and enumerates every feasible path connecting these interaction sites.

## Build

Compile the project with:
```sh
make
```

## Usage

### Demo
Run the ABINOS showcase with:
```sh
make demo
```

### Input format
The input must be a folder named *XX* that contains at least two files:

- XX.mol2 – the substrate description.
- XX_moc*n*.mol2 – each binding motif to connect (one file per motif, *n* is the index).

### Execution
Choose the distance metric (default: Euclidean):

```
export DISTANCE_TYPE="[A*, SSMTA*, Euclidienne, HYBRID]"
```

Reset to the default metric:

```
unset DISTANCE_TYPE
```

Launch the solver:
```
./bin/chemins.exe -i [Input folder] -n [Motif index] -s [Max motifs per path, default 5] -r [Max solutions, default 10] -b [Ban interconnection-tree edge, default 1 → true] -t [One solution per tree, default 1 → true] -p [Enable DIST_PATH_BOUNDARY, default 1 → true] -l [Dynamic path-length pruning, default 0 → false] -g [Interconnection-tree mode, default 1 → store+sort; 0 → on-the-fly]
```

Display inline help:
```sh
./bin/chemins.exe -h
```

Enforce or disable the `DIST_PATH_BOUNDARY` filter globally by exporting:

```
export CAGE_PATH_BOUNDARY=0   # disable pruning
# or
export CAGE_PATH_BOUNDARY=1   # enable pruning
```

This environment variable mirrors the `-p` flag.

To clamp the path length dynamically once a shorter solution exists for a given connection, pass `-l 1` (use `-l 0` to turn it off). This behaves like a branch-and-bound cap during backtracking.

### Cleaning build artifacts

Delete the executable and object files:
```sh
make clean
```
Delete those artifacts plus generated results:
```sh
make mrproper
```

### Visualization

Visualize the produced geometries with PyMOL (or any viewer):
```sh
make open_pymol
```

## Benchmarks for `demos_onePath`

The `demos_onePath/` folder contains lightweight cases that reproduce the path generation pipeline. After compiling (optionally with `make CFLAGS+=-DENABLE_STATS all` to collect statistics), iterate over the datasets using a plain shell loop or your own automation:

```sh
for dataset in demos_onePath/*/; do
    ./bin/chemins.exe -i "$dataset" -n 0 -s 5 -r 10
done
```

Adjust `-n`, `-s`, and `-r` according to how many `*_mocX.mol2` files each dataset holds. The sections below summarize the most relevant toggles when you batch experiments manually, with shell utilities, or via the helper Python scripts in `scripts/`.

### Custom output folders

The executable honors `CAGE_RESULTS_DIR`. Setting `CAGE_RESULTS_DIR=/tmp/cage_benchmark` writes every generated file (results, parameters, timings) into that directory, avoiding collisions between runs and simplifying metric collection.

### "Stats-only" mode (no files written)

- Export `CAGE_DISABLE_WRITE=1` to suppress `.mol2`, `*_parametres.txt`, and `*_time.txt` outputs. Instead, the program prints a single machine-readable line:
    ```text
    RESULT_SUMMARY results=4 best_rmsd=1.234500 worst_rmsd=1.560000 min_path=3 min_path_count=2 min_path_best_rmsd=1.234500 max_path=5 avg_path=4.10 min_atoms=210 min_atoms_count=1 min_atoms_best_rmsd=1.234500 max_atoms=248 avg_atoms=223.75 first_result_ms=312.000 total_time_ms=1874.000 cpu_time_ms=952.000
    ```
- This is ideal to compare builds, chain many inputs, or benchmark without stressing the filesystem.

### Sweeping NUMBER_POSITION_AX1E3 / MAX_POSITION_KEEP_360 / THRESHOLD_ANGLE

The script `scripts/run_stats_sweep.py` recompiles with various `EXTRA_CFLAGS`, launches every dataset in `demos_onePath`, and aggregates the `RESULT_SUMMARY` statistics into `benchmark_runs/stats_sweep.csv` (default). One typical invocation:

```sh
python3 scripts/run_stats_sweep.py \
    --number-branching-positions 3 4 \
    --max-discretized-position 12 24 \
    --threshold-angles 10 15 \
    --distance-types "Euclidienne" "A*" "HYBRID" \
    --size-min 5 --size-max 15 --max-results 1
```

Control the `DIST_PATH_BOUNDARY` filter during sweeps:

- `--path-boundary on` forces `CAGE_PATH_BOUNDARY=1`.
- `--path-boundary off` forces `CAGE_PATH_BOUNDARY=0`.
- `--path-boundary auto` (default) keeps whatever is already configured.

Probe the dynamic `-l` cutoff quickly:

- `--dynamic-path-limit on off` launches two batches (first `-l 1`, then `-l 0`).
- `--dynamic-path-limit on` always adds `-l 1`; `--dynamic-path-limit off` enforces `-l 0`.
- `--dynamic-path-limit auto` leaves the command untouched, relying on compile-time defaults.

Sweep the discretization mode:

- `--discretization-types 0 1` recompiles for each value consumed by the `DISCRITIZATION_TYPE` macro.
- Each value is forwarded through `-DDISCRITIZATION_TYPE=<value>` during compilation.

The CSV now includes `path_boundary`, `dynamic_path_limit`, `discretization_type`, `boundary_allowed`, `boundary_blocked`, `min_path_count`, and `min_path_best_rmsd` so you can track how many solutions hit the minimum length and the best RMSD among them. Additional columns (`best_min_path_length`, `best_max_path_length`, `best_avg_path_length`) expose the optimal length achieved for every path under the dynamic cutoff. Timings (`first_result_ms`, `total_time_ms`, `cpu_time_ms`) are stored directly in milliseconds for fine-grained comparisons.

The script recompiles automatically with `EXTRA_CFLAGS` (e.g., `-DENABLE_STATS -DNUMBER_POSITION_AX1E3=4 ...`). Supply extra definitions with `--extra-cflags "-DMY_FLAG=1"` when needed.

> Tip: you can also pass `EXTRA_CFLAGS` straight to `make`, for example `make EXTRA_CFLAGS="-DENABLE_STATS -DMAX_POSITION_KEEP_360=16"`, to trial specific discretization budgets by hand.

## Comparing interconnection-tree execution modes

The script `scripts/compare_interconnection_modes.py` evaluates both execution strategies (`-g 0` for on-the-fly expansion, `-g 1` for storing and sorting before `generatePaths`) on the `demos_realSubstrat` cases. It also explores several `-b` values to quantify the impact of banned edges without forcing `-t 1`.

Example run:

```sh
python3 scripts/compare_interconnection_modes.py \
    --names ABINOS YILLAG \
    --timeout 600 \
    --size-max 5 \
    --max-results 20 \
    --summary-csv results_compare_modes/summary.csv
```

Key features:

- Enables `CAGE_DISABLE_WRITE=1` by default to collect stats only (use `--keep-results` to generate files) and stores a full log per run under `--log-dir`.
- Creates `results_compare_modes/<substrate>/mode_<g>_b<banned>_<timestamp>/`, pointing `CAGE_RESULTS_DIR` there whenever outputs are preserved.
- Accepts targets as `NAME` or `NAME:moc`, honoring `--default-moc` if no `*_mocX.mol2` exists.
- Parses `RESULT_SUMMARY`, `COMPONENT_SUMMARY`, and the "Interconnection Trees" counter to populate the recap table (result counts, components, min/max sizes, min/max atoms, timings) and appends the same data to the CSV when `--summary-csv` is provided.
- Supports substrate filters (`--names`), forced distance types (`--distance-type`), and arbitrary CLI additions (`--extra-args ...`).

At the end of a session, the script prints a compact table (status, solutions, interconnection trees, components, min/max atoms, total runtime) and, when requested, appends the same metrics to a CSV suitable for plotting or regression. The CSV columns include `min_atom_count`, `max_atom_count`, `interconnection_trees`, `components`, `min_component_vertices`, and `max_component_vertices`, making it easier to track the generated structures in detail.
