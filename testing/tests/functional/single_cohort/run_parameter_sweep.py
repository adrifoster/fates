"""
Parameter sweep driver for the single_cohort functional test.

Runs the already-built SingleCohort_exe once per parameter file in a
directory, each in reduced_output mode with its own uniquely-named output
file (see test_SingleCohort.F90's "COMMAND LINE" header comment) - all
draws share one flat --output-dir rather than each getting its own
subdirectory, since the executable's 3rd command-line argument lets this
script assign each draw a unique filename directly. Extracts each
successful draw's shade-tolerance metrics (see single_cohort_test.py's
SingleCohort.calculate_all_metrics) into a per-draw <stem>_metrics.nc, plus
one combined sweep_metrics.nc across all successful draws, indexed by
parameter file name.

Does NOT rebuild FATES - parameters are read at runtime from the JSON file,
so the executable only needs to be built once for the entire sweep. Build it
first (e.g. `./run_functional_tests.py -t single_cohort` from testing/) if
it doesn't already exist at the default location.

Every output file (raw per-draw netCDF, per-draw metrics, combined
sweep_metrics.nc) is written atomically - to a same-directory .tmp file,
then renamed into place (see _atomic_to_netcdf and run_one_draw) - so
re-running a sweep is safe even while a previous run's output file is open
elsewhere (e.g. loaded into a Jupyter notebook via xr.open_dataset): the
open reader keeps reading its now-unlinked old data undisturbed, and the
write never contends for a lock on a file another process has open.

Output naming, all directly under --output-dir (<stem> = the parameter
file's name without extension):
    <stem>.nc                  raw reduced-output simulation file
    <stem>_metrics.nc          this draw's shade-tolerance metrics
    <stem>_bracket_status.log  _zero_crossing bracket-status summary
    <stem>_error.log           only written if the draw failed
    manifest.json              per-draw success/failure record
    sweep_metrics.nc           every successful draw's metrics, concatenated
                               along a new "draw" dimension (coord = <stem>)

Usage:
    python run_parameter_sweep.py --param-dir /path/to/params/

    python run_parameter_sweep.py --param-dir /path/to/params/ \\
        --pattern "draw_*.json" --output-dir /path/to/sweep_output \\
        --jobs 8

    # rebuild sweep_metrics.nc from existing <stem>_metrics.nc files only,
    # without re-running any simulations (e.g. after the combine step
    # itself failed or was interrupted, but every draw already finished)
    python run_parameter_sweep.py --output-dir /path/to/sweep_output --combine-only
"""
import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import xarray as xr

# both insertions are needed: single_cohort_test.py itself (this directory),
# and the `framework` package it imports from (testing/, three levels up) -
# see single_cohort_test.py's own `from framework...` imports
_THIS_DIR = Path(__file__).resolve().parent
_TESTING_DIR = _THIS_DIR.parents[2]
_REPO_ROOT = _THIS_DIR.parents[3]
sys.path.insert(0, str(_TESTING_DIR))
sys.path.insert(0, str(_THIS_DIR))
from single_cohort_test import SingleCohort  # noqa: E402  (path setup must come first)

_DEFAULT_EXE = _REPO_ROOT / "_build" / "testing" / "fates_single_cohort_ftest" / "SingleCohort_exe"
_DEFAULT_OUTPUT_DIR = _REPO_ROOT / "_sweep_output"

# JSONParameterUtilsMod.F90's JSON parser has been observed to enter a
# runaway retry loop on malformed input (JSONFindTag repeating "could not
# find next tag" ~145,000 times before finally crashing) rather than
# failing fast - cap how much of a failed draw's captured stdout/stderr
# gets written to disk so one degenerate parameter file can't blow up
# error-log size (and, at sweep scale, total disk usage) the way it did in
# testing here
_MAX_LOG_CHARS = 20_000


def _truncate_log(text: str) -> str:
    """Keeps the head and tail of `text`, dropping the middle if it exceeds
    _MAX_LOG_CHARS - the crash signature is almost always at the very start
    (the first error) or the very end (the final traceback/runtime error),
    so both are worth keeping even when a runaway loop bloats the middle

    Args:
        text (str): captured stdout/stderr to (possibly) truncate

    Returns:
        str: `text` unchanged if short enough, otherwise head + a marker +
            tail
    """
    if len(text) <= _MAX_LOG_CHARS:
        return text
    half = _MAX_LOG_CHARS // 2
    omitted = len(text) - _MAX_LOG_CHARS
    return (
        f"{text[:half]}\n"
        f"... [{omitted} characters omitted - see _MAX_LOG_CHARS in this script] ...\n"
        f"{text[-half:]}"
    )


def _atomic_to_netcdf(dataset: xr.Dataset, path: Path) -> None:
    """Writes `dataset` to `path` atomically: to a same-directory temp file
    first, then os.replace (rename) into place. A process that already has
    `path` open (e.g. a Jupyter kernel that loaded it via xr.open_dataset -
    the exact scenario that broke this script the first time: netCDF4/HDF5's
    file locking made a plain in-place to_netcdf raise a PermissionError
    while another process held the file open) keeps reading its now-unlinked
    old inode undisturbed, and this write never contends for a lock on a
    file another process has open, since it's writing to a brand-new path
    until the very last, near-instantaneous rename step

    Args:
        dataset (xr.Dataset): dataset to write
        path (Path): destination path to atomically replace
    """
    tmp_path = path.with_name(path.name + ".tmp")
    dataset.to_netcdf(tmp_path)
    os.replace(tmp_path, path)


def run_one_draw(exe_path: Path, param_file: Path, out_path: Path, timeout: float) -> tuple[bool, str]:
    """Runs one parameter draw in reduced_output mode, writing to a
    caller-chosen output path (test_SingleCohort.F90's optional 3rd
    command-line argument) rather than relying on the driver's default
    filename - this is what lets every draw share one flat --output-dir
    instead of needing a directory each. Invokes the executable by absolute
    path rather than copying it per draw (unlike FunctionalTest.run's
    single-test convention) - copying a multi-MB binary for every one of
    potentially thousands of draws is pure waste when they all use the
    identical, already-built executable.

    The executable itself is told to write to a temp path, then renamed into
    place with os.replace only on success (same atomicity rationale as
    _atomic_to_netcdf) - this means out_path only ever appears/changes as a
    complete, successful result: a crash, timeout, or failed run never
    leaves a partial/corrupt file at out_path, and a process with a
    previous run's out_path already open is never disturbed by a rerun

    Args:
        exe_path (Path): path to the already-built SingleCohort_exe
        param_file (Path): parameter file for this draw
        out_path (Path): where this draw should write its reduced-output
            netCDF - must be unique across draws sharing an --output-dir
        timeout (float): seconds to allow before killing a stuck run - this
            driver's whole design premise is "seconds per run" (see
            test_SingleCohort.F90's header comment), so a hang here is
            anomalous (e.g. a degenerate LHC draw sending the Ci-bisection
            solver into a pathological loop) and should fail that one draw
            rather than block the rest of the sweep indefinitely

    Returns:
        (bool, str): whether the run succeeded, and combined stdout/stderr
            (or a timeout message)
    """
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    try:
        result = subprocess.run(
            [str(exe_path), str(param_file.resolve()), "reduced_output", str(tmp_path)],
            cwd=out_path.parent, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        log = (exc.stdout or "") + (exc.stderr or "")
        tmp_path.unlink(missing_ok=True)
        return False, f"TIMED OUT after {timeout}s\n{log}"

    if result.returncode != 0:
        tmp_path.unlink(missing_ok=True)
        return False, result.stdout + result.stderr

    os.replace(tmp_path, out_path)
    return True, result.stdout + result.stderr


def extract_metrics(out_path: Path) -> tuple[xr.Dataset, str]:
    """Opens a draw's reduced-output file and computes its shade-tolerance
    metrics. calculate_all_metrics prints a bracket-status summary per
    metric (see _print_bracket_status) - useful for one interactive run,
    but noisy across a whole sweep - so it's captured here instead of left
    on stdout, and returned for the caller to log per-draw rather than
    discard: LHC/Morris draws are expected to hit degenerate/out-of-range
    curves far more often than the single default-parameter run this was
    built against did, and that's exactly the signal this status summary
    was designed to surface

    Args:
        out_path (Path): this draw's reduced-output netCDF (see run_one_draw)

    Returns:
        (xarray Dataset, str): this draw's metrics (see
            SingleCohort.calculate_all_metrics), and the captured
            bracket-status log text
    """
    data = xr.open_dataset(out_path)
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        metrics = SingleCohort.calculate_all_metrics(data)
    data.close()
    return metrics, captured.getvalue()


def combine_metrics(output_dir: Path) -> Path | None:
    """(Re)builds sweep_metrics.nc from whatever <stem>_metrics.nc files
    already exist in output_dir, without re-running any simulations or
    re-extracting metrics from the raw per-draw netCDFs - for recovering
    from a crash/interruption after the per-draw work already finished (see
    --combine-only), or for folding in draws added since the last combine.
    Excludes sweep_metrics.nc itself from the glob - its name also happens
    to end in "_metrics.nc"

    Args:
        output_dir (Path): directory containing <stem>_metrics.nc files

    Returns:
        Path | None: path to the written sweep_metrics.nc, or None if no
            per-draw metrics files were found
    """
    metrics_files = sorted(
        f for f in output_dir.glob("*_metrics.nc") if f.name != "sweep_metrics.nc"
    )
    if not metrics_files:
        return None

    all_metrics = []
    for f in metrics_files:
        stem = f.name[: -len("_metrics.nc")]
        data = xr.open_dataset(f)
        all_metrics.append(data.expand_dims(draw=[stem]))
        data.close()

    combined = xr.concat(all_metrics, dim="draw")
    combined_path = output_dir / "sweep_metrics.nc"
    _atomic_to_netcdf(combined, combined_path)
    return combined_path


def main():
    """Parses arguments, runs every parameter file in --param-dir through
    SingleCohort_exe in reduced_output mode (one flat --output-dir, unique
    filenames per draw), and writes per-draw plus combined shade-tolerance
    metrics
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--param-dir", type=Path,
        help="Directory containing the parameter JSON files to sweep over. "
             "Required unless --combine-only is given",
    )
    parser.add_argument(
        "--pattern", default="*.json",
        help="Glob pattern for parameter files within --param-dir (default: *.json)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=_DEFAULT_OUTPUT_DIR,
        help="Directory to write all sweep output to (flat, not one "
             "subdirectory per draw) - kept separate from the standard "
             "functional-test run directory (_run/) so a sweep never "
             f"collides with a normal test run (default: {_DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--exe", type=Path, default=_DEFAULT_EXE,
        help=f"Path to the already-built SingleCohort_exe (default: {_DEFAULT_EXE})",
    )
    parser.add_argument(
        "--jobs", type=int, default=1,
        help="Number of draws to run concurrently (default: 1, sequential). "
             "Safe to raise: every draw writes to its own uniquely-named "
             "file, so concurrent draws cannot collide",
    )
    parser.add_argument(
        "--timeout", type=float, default=300.0,
        help="Seconds to allow a single draw before treating it as failed "
             "(default: 300)",
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip draws whose output file already exists - resume a "
             "partially-completed sweep without redoing finished draws",
    )
    parser.add_argument(
        "--combine-only", action="store_true",
        help="Skip running any simulations entirely - just (re)build "
             "sweep_metrics.nc from whatever <stem>_metrics.nc files "
             "already exist in --output-dir. For recovering after the "
             "per-draw work already finished but the final combine step "
             "failed or was interrupted (--param-dir/--pattern/--exe are "
             "ignored in this mode)",
    )
    args = parser.parse_args()

    if args.combine_only:
        combined_path = combine_metrics(args.output_dir)
        if combined_path is None:
            raise FileNotFoundError(
                f"No <stem>_metrics.nc files found in {args.output_dir} - "
                "nothing to combine"
            )
        print(f"Combined metrics written to {combined_path}")
        return

    if args.param_dir is None:
        raise ValueError("--param-dir is required unless --combine-only is given")

    if not args.exe.exists():
        raise FileNotFoundError(
            f"SingleCohort_exe not found at {args.exe} - build it first, e.g. "
            "`./run_functional_tests.py -t single_cohort` from testing/ "
            "(parameters are read at runtime, so the executable only needs "
            "to be built once for the whole sweep)."
        )

    param_files = sorted(args.param_dir.glob(args.pattern))
    if not param_files:
        raise FileNotFoundError(
            f"No parameter files matching '{args.pattern}' found in {args.param_dir}"
        )
    stems = [pf.stem for pf in param_files]
    if len(set(stems)) != len(stems):
        raise ValueError(
            "Two or more parameter files share the same stem (filename "
            "without extension) - every draw's output filename is derived "
            "from its parameter file's stem, so these would collide. "
            "Rename the parameter files so every stem is unique."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Sweeping {len(param_files)} parameter files -> {args.output_dir}")

    def process_one(param_file: Path):
        out_path = args.output_dir / f"{param_file.stem}.nc"
        if args.skip_existing and out_path.exists():
            return param_file, True, "skipped (existing output)"
        return (param_file, *run_one_draw(args.exe, param_file, out_path, args.timeout))

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(process_one, pf): pf for pf in param_files}
        for i, future in enumerate(as_completed(futures), start=1):
            param_file, ok, log = future.result()
            print(f"[{i}/{len(param_files)}] {param_file.name}: {'ok' if ok else 'FAILED'}")
            results.append((param_file, ok, log))

    # process in the original sorted order regardless of completion order,
    # so the manifest/combined dataset are deterministic run to run
    results.sort(key=lambda r: param_files.index(r[0]))

    manifest = []
    all_metrics = []
    for param_file, ok, log in results:
        out_path = args.output_dir / f"{param_file.stem}.nc"
        entry = {"param_file": param_file.name, "output_file": str(out_path), "success": ok}

        if not ok:
            (args.output_dir / f"{param_file.stem}_error.log").write_text(_truncate_log(log))
            entry["error"] = f"run failed - see {param_file.stem}_error.log"
            manifest.append(entry)
            continue

        try:
            metrics, bracket_log = extract_metrics(out_path)
        except Exception as exc:  # noqa: BLE001 - one bad draw must not kill the sweep
            entry["success"] = False
            entry["error"] = f"metrics extraction failed: {exc}"
            manifest.append(entry)
            print(f"FAILED (metrics): {param_file.name} - {exc}")
            continue

        (args.output_dir / f"{param_file.stem}_bracket_status.log").write_text(bracket_log)
        _atomic_to_netcdf(metrics, args.output_dir / f"{param_file.stem}_metrics.nc")
        all_metrics.append(metrics.expand_dims(draw=[param_file.stem]))
        manifest.append(entry)

    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    n_ok = sum(entry["success"] for entry in manifest)
    print(f"\n{n_ok}/{len(param_files)} draws succeeded. Manifest: {manifest_path}")

    if all_metrics:
        combined = xr.concat(all_metrics, dim="draw")
        combined_path = args.output_dir / "sweep_metrics.nc"
        _atomic_to_netcdf(combined, combined_path)
        print(f"Combined metrics for {len(all_metrics)} draws written to {combined_path}")


if __name__ == "__main__":
    main()
