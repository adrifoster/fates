"""
Notebook-friendly helpers for analyzing single_cohort parameter-sweep output.

Complements run_parameter_sweep.py (which runs the sweep) and
single_cohort_test.py's SingleCohort.calculate_all_metrics (which defines
the actual metrics): these functions are for interactive prototyping against
an existing sweep's --output-dir - reading the raw per-draw simulation
files, recomputing metrics directly from them (rather than the
already-extracted <stem>_metrics.nc files - useful if you want to try a
different `year` or inspect fields calculate_all_metrics doesn't bundle),
and joining in a sweep's parameter "key" file (e.g. trop_key.csv) so metrics
can be related back to the parameter values that produced them.

Typical notebook usage:
    import sys
    sys.path.insert(0, "/path/to/testing")
    sys.path.insert(0, "/path/to/testing/tests/functional/single_cohort")
    from sweep_analysis import calculate_metrics_from_raw, load_key, attach_key

    metrics = calculate_metrics_from_raw("/path/to/_sweep_output")
    key = load_key("/path/to/param_out")
    metrics = attach_key(metrics, key)

    # e.g. final-year LCPplant at full sun vs. the key's vcmax25top column
    metrics.plot.scatter(x="fates_leaf_vcmax25top", y="lcp_plant")
"""
import contextlib
import io
import sys
from pathlib import Path

import pandas as pd
import xarray as xr

# both insertions are needed: single_cohort_test.py itself (this directory),
# and the `framework` package it imports from (testing/, three levels up) -
# see single_cohort_test.py's own `from framework...` imports
_THIS_DIR = Path(__file__).resolve().parent
_TESTING_DIR = _THIS_DIR.parents[2]
sys.path.insert(0, str(_TESTING_DIR))
sys.path.insert(0, str(_THIS_DIR))
from single_cohort_test import SingleCohort  # noqa: E402  (path setup must come first)


def _is_raw_output(path: Path) -> bool:
    """True if `path` looks like a raw per-draw simulation output file
    (<stem>.nc, written by run_parameter_sweep.py's run_one_draw) rather
    than a derived sweep artifact (<stem>_metrics.nc, sweep_metrics.nc, or a
    leftover .tmp from an interrupted atomic write)

    Args:
        path (Path): candidate file

    Returns:
        bool: whether this looks like a raw simulation output file
    """
    return (
        path.suffix == ".nc"
        and not path.name.endswith("_metrics.nc")
        and path.name != "sweep_metrics.nc"
    )


def find_raw_outputs(output_dir, pattern: str = "*.nc") -> dict:
    """Finds every raw per-draw simulation output file in a sweep's
    --output-dir, keyed by draw name (the parameter file's stem - see
    run_parameter_sweep.py's output-naming convention)

    Args:
        output_dir (str | Path): a run_parameter_sweep.py --output-dir
        pattern (str): glob pattern to search with (default: "*.nc")

    Returns:
        dict[str, Path]: {draw_name: path}, sorted by draw name
    """
    output_dir = Path(output_dir)
    files = sorted(f for f in output_dir.glob(pattern) if _is_raw_output(f))
    if not files:
        raise FileNotFoundError(
            f"No raw per-draw output files found in {output_dir} (pattern '{pattern}')"
        )
    return {f.stem: f for f in files}


def load_raw_outputs(output_dir, pattern: str = "*.nc") -> dict:
    """Opens every raw per-draw simulation output file in a sweep's
    --output-dir, keyed by draw name - for interactively inspecting a few
    specific draws' full (reduced_output) fields in a notebook (e.g.
    plotting a compensation-point trajectory for one draw the way
    single_cohort_test.py's own plotting methods do, via
    SingleCohort.plot_compensation_point_trajectory(dataset, False, None)).

    Opens lazily (xarray only reads array data on access), but each returned
    Dataset keeps its own file handle open until you call .close() on it.
    Fine for pulling out a handful of draws to inspect closely, but for
    computing something over EVERY draw in a large sweep, use
    calculate_metrics_from_raw instead - it opens and closes one file at a
    time and never holds hundreds of file handles open simultaneously (worth
    caring about here specifically: an open file handle blocking a later
    write is exactly what broke run_parameter_sweep.py's sweep_metrics.nc
    write earlier in this project)

    Args:
        output_dir (str | Path): a run_parameter_sweep.py --output-dir
        pattern (str): glob pattern to search with (default: "*.nc")

    Returns:
        dict[str, xr.Dataset]: {draw_name: dataset}
    """
    return {
        name: xr.open_dataset(path)
        for name, path in find_raw_outputs(output_dir, pattern).items()
    }


def calculate_metrics_from_raw(
    output_dir, pattern: str = "*.nc", year: int = None, verbose: bool = False,
) -> xr.Dataset:
    """Recomputes shade-tolerance metrics (SingleCohort.calculate_all_metrics)
    directly from every raw per-draw simulation output file in a sweep's
    --output-dir, combined into one Dataset along a new "draw" dimension -
    the same shape of result run_parameter_sweep.py's own sweep_metrics.nc
    holds, but computed fresh from the full per-draw files rather than read
    back from the already-extracted <stem>_metrics.nc files. Useful for
    prototyping a different `year` for the Onoda metrics (LIE_LA/LIE_M/LUE
    default to the final simulated year - see calculate_all_metrics) without
    changing run_parameter_sweep.py or re-running any simulations.

    Opens and closes one raw file at a time (never holds more than one file
    handle open at once) - safe to run over a sweep with hundreds of draws.

    Args:
        output_dir (str | Path): a run_parameter_sweep.py --output-dir
        pattern (str): glob pattern to search with (default: "*.nc")
        year (int): 1-indexed simulated year for LIE_LA/LIE_M/LUE (default:
            None, the final simulated year - see calculate_all_metrics)
        verbose (bool): if True, let calculate_all_metrics' per-draw
            bracket-status summary (see _print_bracket_status) print to
            stdout as usual - useful when checking one or a few draws
            closely, but produces several lines of output per draw, per
            metric, so defaults to suppressed here across a whole sweep
            (matching run_parameter_sweep.py's own extract_metrics)

    Returns:
        xr.Dataset: lcp_plant/lcp_leaf (draw, year, light_level),
            growth_floor_edge (draw, year), lie_la/lie_m/lue
            (draw, light_level)
    """
    raw_files = find_raw_outputs(output_dir, pattern)

    all_metrics = []
    for draw_name, path in raw_files.items():
        with xr.open_dataset(path) as data:
            if verbose:
                metrics = SingleCohort.calculate_all_metrics(data, year=year)
            else:
                with contextlib.redirect_stdout(io.StringIO()):
                    metrics = SingleCohort.calculate_all_metrics(data, year=year)
        all_metrics.append(metrics.expand_dims(draw=[draw_name]))

    return xr.concat(all_metrics, dim="draw")


def load_key(param_dir, filename: str = None) -> pd.DataFrame:
    """Loads a sweep's parameter "key" CSV (e.g. trop_key.csv - the
    ensemble-to-parameter-value lookup a generator script writes alongside
    its parameter files), indexed by draw name.

    Note: the values in this file are whatever the generator wrote - in the
    sweep this has been used with so far, they are normalized [0,1]
    design-space coordinates, NOT the actual physical values substituted
    into each draw's JSON (confirmed by direct inspection: a key entry of
    0.176 for fates_leaf_slatop corresponded to an actual substituted value
    of 0.0206, not 0.176 m2/gC). Treat these as design-space coordinates for
    correlation/ranking purposes - fine for the sanity-check style analysis
    this module is meant for, but not literal parameter values - unless you
    know your specific generator's key file records physical values instead.

    Args:
        param_dir (str | Path): directory containing the key CSV (typically
            the same --param-dir passed to run_parameter_sweep.py)
        filename (str): explicit key filename; if None, looks for exactly
            one "*_key.csv" file in param_dir

    Returns:
        pd.DataFrame: key values, indexed by draw/ensemble name
    """
    param_dir = Path(param_dir)
    if filename is not None:
        path = param_dir / filename
    else:
        matches = sorted(param_dir.glob("*_key.csv"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected exactly one *_key.csv in {param_dir}, found "
                f"{len(matches)}: {matches}. Pass filename= explicitly."
            )
        path = matches[0]

    key = pd.read_csv(path)
    index_col = "ensemble" if "ensemble" in key.columns else key.columns[0]
    return key.set_index(index_col)


def attach_key(metrics: xr.Dataset, key: pd.DataFrame) -> xr.Dataset:
    """Attaches a sweep's parameter key values (see load_key) to a metrics
    Dataset (see calculate_metrics_from_raw, or xr.open_dataset on
    run_parameter_sweep.py's sweep_metrics.nc) as non-dimensional
    coordinates alongside the existing "draw" dimension - e.g.
    metrics.fates_leaf_vcmax25top becomes directly usable for
    plotting/selecting/correlating against any metric, without conflating
    input parameter values with computed output metrics (which stay in
    metrics.data_vars, untouched).

    Aligns by draw name: a draw present in `metrics` but missing from `key`
    (or vice versa) gets NaN on whichever side is missing, rather than
    raising - so a draw that failed the sweep (and so has no metrics)
    simply carries NaN parameter values rather than disappearing, and a key
    entry for a draw that isn't in `metrics` is silently unused.

    Args:
        metrics (xr.Dataset): output of calculate_metrics_from_raw, or
            xr.open_dataset(".../sweep_metrics.nc")
        key (pd.DataFrame): output of load_key

    Returns:
        xr.Dataset: `metrics` with one new coordinate per column of `key`
    """
    key_aligned = key.reindex(metrics.draw.values)
    coords = {col: ("draw", key_aligned[col].values) for col in key_aligned.columns}
    return metrics.assign_coords(coords)
