## Gen_ICT: Interconnection Tree Enumeration

### Purpose

`Gen_ICT` enumerates every interconnection tree built from a partition. The executable prints one CSV line per run containing:

- `solver_input`, `num_vertices`, `num_components`, and the minimum/maximum component sizes
- `num_trees`, `elapsed_ms`, `delay_ms`
- `elapsed_sort_ms` and `elapsed_stock_ms` (weighted mode only)
- `status` (`ok` or `timeout`).

### Compilation

The `Makefile` exposes the `Gen_ICT` target.

```bash
# Fast build (default Makefile flags)
make Gen_ICT

# Clean and rebuild
make clean && make Gen_ICT
```

#### Enable the weighted build

Add `-DWEIGHTED` to the compilation flags to force the program to:

1. store every generated tree,
2. assign a weight to each tree edge,
3. sort the solutions at the end.

```bash
make clean
make CFLAGS="-O2 -DWEIGHTED" Gen_ICT
```

In this mode, the columns `elapsed_sort_ms` and `elapsed_stock_ms` report the time spent sorting and storing the weighted trees. Without `-DWEIGHTED`, both values remain `0.0`.

### Manual execution

```bash
./Gen_ICT tests/realYARZUN03_YOHPIU_YOHPOA.txt            # simple run
./Gen_ICT tests/realYARZUN03_YOHPIU_YOHPOA.txt 600        # stop after 600 s
```

When the time limit is reached, the program exits cleanly, keeps the counters already produced, and sets `status=timeout`.

### Automated benchmarks

The script `scripts/run_bench.py` runs `Gen_ICT` on a set of **existing** files.

```bash
python3 scripts/run_bench.py \
  --inputs tests_article/*.txt \
  --reps 1 \
  --time-limit 120 \
  --output results.csv
```

Available options:

- `--inputs`: one or more glob patterns to pick the partitions.
- `--reps`: measured repetitions per file (default `1`).
- `--warmup`: untimed warmup runs (default `0`).
- `--quiet`: hide progress logs.
- `--time-limit`: forwarded to `Gen_ICT`.
- `--output`: CSV file to append to (the header is added automatically).

Each row of `results.csv` now begins with the **file name** (`label`) followed by the repetition index. The remaining columns come directly from the `Gen_ICT` output, including `elapsed_sort_ms` and `elapsed_stock_ms`, which are non-zero only when the binary was compiled with `-DWEIGHTED`.

### Input format

1. First line: `N K` with `N` vertices and `K` components.
2. Next `N` lines: the component index (0-based) assigned to each vertex.

Example:

```
5 2
0
0
0
1
1
```

In this example, 5 vertices are split into 2 parts: the first three vertices belong to component `0`, the last two to component `1`.
