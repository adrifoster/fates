# Running a parameter sweep over `single_cohort`

`run_parameter_sweep.py` runs the `single_cohort` functional test once per
parameter file in a directory - built for the shade-tolerance
tolerance-inversion project's eventual Morris/LHC screening (see
`PARAMETER_SENSITIVITY.md` for which parameters actually affect this test's
output at all). It reuses the already-built executable across every draw
(parameters are read at runtime, not compiled in) and extracts each draw's
shade-tolerance metrics (`SingleCohort.calculate_all_metrics` - LCPplant,
LCPleaf, the growth-floor edge, LIE_LA, LIE_M, LUE) into one combined
dataset, without needing to keep every draw's multi-MB simulation output
around to do so.

## Prerequisite: build once

The executable is not rebuilt per draw - parameters are read from the JSON
file at runtime, so the same binary works for every draw. Build it once
before running a sweep:

```bash
conda activate fates_testing
cd testing
./run_functional_tests.py -t single_cohort
```

This leaves `SingleCohort_exe` at
`_build/testing/fates_single_cohort_ftest/SingleCohort_exe`, which is the
sweep script's default `--exe` location.

## Basic usage

```bash
cd testing/tests/functional/single_cohort
python run_parameter_sweep.py --param-dir /path/to/your/parameter/files/
```

That's sufficient if your parameter files are named uniquely (any
extension-having filenames work; the glob defaults to `*.json`) and you're
fine with the default output location (`_sweep_output/` at the repo root).

Full example, matching a real run:

```bash
python run_parameter_sweep.py \
    --param-dir ~/Documents/02_Projects/13_Shade_Tolerance/param_out \
    --pattern "trop_*.json" \
    --output-dir /path/to/sweep_output \
    --jobs 8
```

## CLI reference

| Flag | Default | Meaning |
|---|---|---|
| `--param-dir` | *(required unless `--combine-only`)* | Directory containing the parameter JSON files to sweep over |
| `--pattern` | `*.json` | Glob pattern for parameter files within `--param-dir` |
| `--output-dir` | `<repo_root>/_sweep_output` | Where all sweep output goes (flat, not one subdirectory per draw) - deliberately separate from `_run/`, the standard functional-test run directory, so a sweep never collides with a normal test run |
| `--exe` | `<repo_root>/_build/testing/fates_single_cohort_ftest/SingleCohort_exe` | Path to the already-built executable |
| `--jobs` | `1` | Number of draws to run concurrently. Safe to raise - every draw writes to its own uniquely-named file, so concurrent draws cannot collide |
| `--timeout` | `300` (seconds) | How long to allow a single draw before treating it as failed, rather than letting one hung/degenerate draw block the whole sweep |
| `--skip-existing` | off | Skip draws whose output file already exists - lets you resume a partially-completed sweep without redoing finished draws |
| `--combine-only` | off | Skip running any simulations entirely - just (re)build `sweep_metrics.nc` from whatever `<stem>_metrics.nc` files already exist in `--output-dir`. `--param-dir`/`--pattern`/`--exe` are ignored in this mode. For recovering when every draw already finished but the final combine step failed or was interrupted - see Gotchas |

## Output

Everything lands directly in `--output-dir` (flat - no per-draw
subdirectories), named from each parameter file's stem (filename without
extension):

| File | Contents |
|---|---|
| `<stem>.nc` | Raw reduced-output simulation file for that draw (see `FatesTestHistoryMod.F90`'s module header for exactly which variables `reduced_output` mode keeps) |
| `<stem>_metrics.nc` | That draw's shade-tolerance metrics (`SingleCohort.calculate_all_metrics`'s output - lcp_plant/lcp_leaf `(year, light_level)`, growth_floor_edge `(year,)`, lie_la/lie_m/lue `(light_level,)`, the last three from the final simulated year) |
| `<stem>_bracket_status.log` | `_zero_crossing` bracket-status summary for that draw (how many of the 250 LCP sweeps/10 growth-floor-edge sweeps cleanly bracketed vs. fell outside the swept range or hit a degenerate curve) |
| `<stem>_error.log` | Only written if the draw failed - captured stdout/stderr, truncated to ~20KB (head + tail) if the underlying error was unusually verbose (see Gotchas below) |
| `manifest.json` | One entry per parameter file: success/failure and, on failure, which log to check |
| `sweep_metrics.nc` | Every successful draw's metrics concatenated along a new `draw` dimension (coordinate = parameter file stem) - the one file you actually want for downstream analysis (GP emulation, history matching, etc.) |

Loading the combined result:

```python
import xarray as xr
metrics = xr.open_dataset("_sweep_output/sweep_metrics.nc")
metrics["lcp_plant"].isel(year=-1, light_level=-1)  # final-year LCPplant at full sun, per draw
```

If you also have a key file mapping draw name to the actual parameter
values used (as generated alongside the parameter files, e.g. `*_key.csv`),
join it against `metrics.draw` to relate emergent metrics back to inputs -
see the worked example below.

## Worked example

Sweeping 10 draws that perturb `fates_leaf_vcmax25top` at PFT index 0 (see
Gotchas - this has to be PFT index 0, not any other index) against a
`trop_key.csv` recording each draw's actual parameter values:

```python
import xarray as xr
import pandas as pd

metrics = xr.open_dataset("_sweep_output/sweep_metrics.nc")
key = pd.read_csv("param_out/trop_key.csv", index_col="ensemble")

summary = key.copy()
summary["lcp_plant_final"] = metrics["lcp_plant"].isel(year=-1, light_level=-1).to_pandas()
summary["lue_final"] = metrics["lue"].isel(light_level=-1).to_pandas()

print(summary.sort_values("fates_leaf_vcmax25top"))
print(summary["fates_leaf_vcmax25top"].corr(summary["lcp_plant_final"]))  # expect negative
print(summary["fates_leaf_vcmax25top"].corr(summary["lue_final"]))        # expect positive
```

In a real run of this, higher `vcmax25top` correlated with lower `LCPplant`
(r ≈ -0.70) and higher `LUE` (r ≈ 0.98) - both in the physically expected
direction, and a good sanity check that a new parameter set is actually
propagating through to the metrics rather than being silently ignored (see
the next section for why that check matters).

## Gotchas

- **PFT indexing.** `test_SingleCohort.F90` hardcodes `pft = 1`, a
  Fortran 1-indexed constant - this reads **JSON array index 0** (the
  first PFT) in every per-PFT parameter (`fates_leaf_vcmax25top`,
  `fates_leaf_slatop`, etc: arrays are 0-indexed in the JSON file, 1-indexed
  once loaded into the Fortran array). If a parameter-file generator
  perturbs any *other* index, every draw will read the same unperturbed
  default value for that parameter and the sweep will silently do nothing -
  no error, just identical output across every draw. **Always spot-check a
  couple of generated files** (diff against `parameter_files/fates_params_default.json`,
  confirm index 0 is what actually changed) before running a real sweep -
  this exact bug was caught only by inspection, not by any error message.
- **A malformed parameter file fails loudly, but verbosely.** FATES's JSON
  parser (`JSONParameterUtilsMod.F90`) has been observed to enter a runaway
  retry loop on malformed input (~145,000 repeated "could not find next
  tag" lines) before finally crashing with a Fortran runtime error, rather
  than failing fast. The sweep script still handles this correctly (that
  draw is marked failed, the rest of the sweep continues), and truncates
  the resulting log to ~20KB so it doesn't consume excess disk at sweep
  scale - but if you see a `_error.log` that looks like this, the fix is in
  whatever generates your parameter files, not in this script.
- **Output is safe to re-run even while a previous result is open elsewhere.**
  Every output file (raw per-draw netCDF, per-draw metrics, `sweep_metrics.nc`)
  is written atomically - to a same-directory `.tmp` file, then renamed into
  place - rather than overwritten in place. This was added after a real
  failure: netCDF4/HDF5 file locking made `sweep_metrics.nc` raise a
  `PermissionError` mid-write because a Jupyter kernel still had the
  previous run's file open (e.g. after following the "Loading the combined
  result" example above). With the atomic write, that kernel keeps reading
  its now-unlinked old data undisturbed, and the new write always succeeds
  - you no longer need to close notebooks/kernels before re-running a sweep.
  If the very last combine step still fails or gets interrupted for some
  other reason (every per-draw `<stem>_metrics.nc` already exists, but
  `sweep_metrics.nc` doesn't, or is stale), use `--combine-only` to rebuild
  just that final file without re-running any simulations:
  `python run_parameter_sweep.py --output-dir <same dir> --combine-only`.
- **`reduced_output` mode omits most diagnostic output on purpose.** See
  `PARAMETER_SENSITIVITY.md` for the full parameter-usage trace and
  `FatesTestHistoryMod.F90`'s module header for exactly which variables
  survive - if a metric you need isn't in `calculate_all_metrics`, it's
  likely because the underlying variable isn't written in reduced mode at
  all, not a computation gap.
- **Metrics are final-year-only for LIE_LA/LIE_M/LUE.** `calculate_all_metrics`
  reports these three at the last simulated year by default. The
  underlying `calculate_lie_la`/`calculate_lie_m`/`calculate_lue` in
  `single_cohort_test.py` each accept an explicit `year` argument if you
  need a different year's snapshot or want to build your own ontogenetic
  trajectory across the sweep - `calculate_all_metrics` doesn't expose that
  as a sweep-wide option, so call those directly per draw if needed.
