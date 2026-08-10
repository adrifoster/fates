"""
Fits a &fates_test_site namelist file (see FatesTestSiteMod.F90) from a
directory of DATM-format hourly forcing files, the same way the driver's
built-in BCI, Panama defaults were derived from
~/Documents/02_Projects/06_FATES/BCI_datm.

Method (developed by refitting the BCI data itself against
FatesTestSiteMod.F90's documented defaults):
    - latitude_deg is read directly off the data's LATIXY coordinate - it is
      a site descriptor, not a fitted quantity.
    - t_annual_mean/t_annual_amp/t_annual_peak_doy (the annual TBOT cycle)
      and t_diurnal_amp/t_diurnal_peak_hour (the diurnal TBOT cycle) are all
      recovered from a single joint ordinary-least-squares fit of hourly
      TBOT against two harmonics (annual + diurnal):

          TBOT ~ mean + A_ann*cos(w_ann*doy) + B_ann*sin(w_ann*doy)
                      + A_diu*cos(w_diu*hr)  + B_diu*sin(w_diu*hr)

      with w_ann = 2*pi/365, w_diu = pi/12 - matching the exact functional
      form FatesTestEnvironmentMod.F90's SetHour evaluates at runtime. This
      is linear in its 5 coefficients (no nonlinear solver, no initial-guess
      sensitivity); amplitude/peak are then recovered from each harmonic's
      (A, B) coefficient pair via hypot/atan2. day-of-year (doy) is computed
      on a fixed 365-day calendar (Feb 29 dropped if present) to match the
      driver's own hardcoded 365-day year, and hour-of-day (hr) is converted
      from the data's timestamp to local solar hour via the site's longitude
      (UTC hour + longitude/15, wrapped to [0, 24)) - matching the "local
      solar hour" FatesTestLightEnvMod.F90's coszen calculation already
      assumes.

      Refitting BCI_datm this way reproduces t_annual_mean/t_annual_amp/
      t_annual_peak_doy exactly (299.117/0.6527/120.555 vs. the documented
      299.1172/0.6527/120.555) but lands ~0.3 hours off the documented
      t_diurnal_peak_hour (13.74 vs. 13.450) and slightly higher on
      t_diurnal_amp (2.271 vs. 2.2547). That gap traces to a floating-point
      decoding artifact in the BCI files themselves: about a third of their
      "hourly" cftime timestamps decode as e.g. 04:59:59.9996 rather than
      exactly 05:00:00, and naively truncating those to an integer hour (an
      easy mistake, since it silently mislabels a third of the records by an
      hour) happens to land within ~0.04 hours of the documented BCI values.
      This function does the technically-correct thing instead - it keeps
      each timestamp's actual sub-hour remainder (which also fits BCI's own
      data marginally better: RMSE 1.497 vs. 1.510 K) - so a from-scratch
      fit against a clean new-site dataset should be trusted over an exact
      match to BCI's diurnal terms specifically.
    - relative_humidity is deliberately NOT fit from data - BCI's default
      (0.70) is a fixed simplifying assumption, not fit to BCI's actual
      relative humidity (which averages ~89%). This script requires it as an
      explicit --relative-humidity argument for the same reason: it is a
      judgment call about the driver's fixed-RH simplification, not a
      quantity a regression can respond for you. If the input data includes
      an RH variable, its mean/median are printed as a reference point only.

Usage:
    python fit_site_namelist.py --datm-dir /path/to/SITE_datm \\
        --relative-humidity 0.70 --output site_namelist.nml

    # custom file glob and TBOT/RH/lat/lon variable names, for DATM-like
    # data that doesn't use CESM DATM's usual names
    python fit_site_namelist.py --datm-dir /path/to/SITE_datm \\
        --relative-humidity 0.75 --pattern "*.nc" \\
        --tbot-var TBOT --rh-var RH --lat-var LATIXY --lon-var LONGXY
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

# fixed 365-day (no-leap) month lengths - matches the FATES functional-test
# driver's own hardcoded days_per_year=365 (test_SingleCohort.F90), which
# has no leap-day handling at all. Day-of-year is computed from this table
# directly (month/day fields), rather than trusting the input data's own
# possibly-Gregorian dayofyear, so a leap-calendar input site fits onto
# exactly the same 365-day cycle the driver will actually run
_NOLEAP_MONTH_LENGTHS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
_NOLEAP_MONTH_OFFSETS = np.cumsum([0] + _NOLEAP_MONTH_LENGTHS[:-1])

_W_ANN = 2.0 * np.pi / 365.0  # annual angular frequency [rad/day]
_W_DIU = np.pi / 12.0         # diurnal angular frequency [rad/hour]

_NAMELIST_VAR_ORDER = [
    "latitude_deg", "t_annual_mean", "t_annual_amp", "t_annual_peak_doy",
    "t_diurnal_amp", "t_diurnal_peak_hour", "relative_humidity",
]


def _noleap_doy(times) -> np.ndarray:
    """Computes day-of-year (1-365) from an array of datetime-like objects
    (cftime or numpy datetime64), on a fixed 365-day calendar - see
    _NOLEAP_MONTH_LENGTHS. Timestamps falling on Feb 29 have no representation
    in this fixed calendar and are dropped (see the caller, fit_site).

    Args:
        times: array-like of objects exposing .month and .day (cftime
            datetimes, or pandas/numpy datetime64 converted via pandas)

    Returns:
        (np.ndarray, np.ndarray): (doy, keep_mask) - doy is only valid where
            keep_mask is True (i.e. not a dropped Feb 29)
    """
    months = np.array([t.month for t in times])
    days = np.array([t.day for t in times])
    keep_mask = ~((months == 2) & (days == 29))
    doy = _NOLEAP_MONTH_OFFSETS[months - 1] + days
    return doy.astype(float), keep_mask


def _local_solar_hour(times, longitude_deg: float) -> np.ndarray:
    """Converts timestamps (assumed UTC, the standard CESM DATM convention)
    to local solar hour-of-day [0, 24), via a simple longitude/15 offset -
    matching FatesTestLightEnvMod.F90's coszen calculation, which takes
    hour_of_day as local solar hour with no further correction. This omits
    the equation-of-time correction (which shifts true solar noon by up to
    ~16 minutes across the year) - on the BCI data, adding it changed
    t_diurnal_peak_hour by only ~0.008 hours (under 30 seconds), well
    inside this single-harmonic fit's own residual, so the added complexity
    isn't worth it here. Keeps each timestamp's actual minute/second
    remainder (not just the truncated hour) - some DATM files have hourly
    timestamps that decode with floating-point jitter (e.g. 04:59:59.9996
    instead of exactly 05:00:00); truncating those would silently mislabel
    the affected records by a full hour.

    Args:
        times: array-like of objects exposing .hour, .minute, .second
        longitude_deg (float): site longitude, 0-360 or -180-180 [deg E]

    Returns:
        np.ndarray: local solar hour-of-day [0, 24)
    """
    lon_signed = ((longitude_deg + 180.0) % 360.0) - 180.0
    offset_hours = lon_signed / 15.0
    hour_utc = np.array([t.hour + t.minute / 60.0 + t.second / 3600.0 for t in times])
    return (hour_utc + offset_hours) % 24.0


def fit_site(
    datm_dir: Path,
    pattern: str = "*.nc",
    tbot_var: str = "TBOT",
    rh_var: str = "RH",
    lat_var: str = "LATIXY",
    lon_var: str = "LONGXY",
) -> dict:
    """Fits the annual/diurnal TBOT-cycle namelist values from a directory of
    DATM-format hourly forcing files (one point, split across any number of
    files - e.g. one per month, matching BCI_datm's layout).

    Args:
        datm_dir (Path): directory containing the site's DATM netCDF files
        pattern (str, optional): glob pattern for files within datm_dir.
            Defaults to "*.nc".
        tbot_var (str, optional): name of the hourly temperature variable
            [K]. Defaults to "TBOT" (CESM DATM's convention).
        rh_var (str, optional): name of the hourly relative humidity
            variable [%], used only to print a reference value - not fit.
            Defaults to "RH". Set to None/"" if the data has no RH variable.
        lat_var (str, optional): name of the time-invariant latitude
            coordinate [deg N]. Defaults to "LATIXY".
        lon_var (str, optional): name of the time-invariant longitude
            coordinate [deg E, 0-360 or -180-180]. Defaults to "LONGXY".

    Returns:
        dict: fit results - the five fitted namelist values plus
            latitude_deg, n_hours (records used), n_dropped_feb29,
            rmse_k (fit residual RMSE [K]), and (if rh_var data was found)
            rh_mean_pct/rh_median_pct for reference
    """
    files = sorted(Path(datm_dir).glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching '{pattern}' found in {datm_dir}")

    all_times, all_tbot, all_rh = [], [], []
    lat = lon = None
    for f in files:
        with xr.open_dataset(f) as ds:
            if tbot_var not in ds:
                raise KeyError(f"{f}: no '{tbot_var}' variable (see --tbot-var)")
            site_lat = float(np.asarray(ds[lat_var]).reshape(-1)[0])
            site_lon = float(np.asarray(ds[lon_var]).reshape(-1)[0])
            if lat is None:
                lat, lon = site_lat, site_lon
            elif abs(site_lat - lat) > 1e-6 or abs(site_lon - lon) > 1e-6:
                raise ValueError(
                    f"{f}: {lat_var}/{lon_var} ({site_lat}, {site_lon}) disagree with "
                    f"the first file's ({lat}, {lon}) - is this one site's data?"
                )

            all_times.append(np.asarray(ds["time"].values))
            all_tbot.append(np.asarray(ds[tbot_var]).reshape(len(ds["time"]), -1)[:, 0])
            if rh_var and rh_var in ds:
                all_rh.append(np.asarray(ds[rh_var]).reshape(len(ds["time"]), -1)[:, 0])

    times = np.concatenate(all_times)
    tbot = np.concatenate(all_tbot).astype(float)

    doy, keep_mask = _noleap_doy(times)
    hour = _local_solar_hour(times, lon)
    finite_mask = keep_mask & np.isfinite(tbot)
    n_dropped_feb29 = int((~keep_mask).sum())

    doy, hour, tbot_fit = doy[finite_mask], hour[finite_mask], tbot[finite_mask]

    design = np.column_stack([
        np.ones_like(doy),
        np.cos(_W_ANN * doy), np.sin(_W_ANN * doy),
        np.cos(_W_DIU * hour), np.sin(_W_DIU * hour),
    ])
    coeffs, _, _, _ = np.linalg.lstsq(design, tbot_fit, rcond=None)
    c0, c1, c2, c3, c4 = coeffs
    resid = tbot_fit - design @ coeffs

    result = {
        "latitude_deg": lat,
        "t_annual_mean": float(c0),
        "t_annual_amp": float(np.hypot(c1, c2)),
        "t_annual_peak_doy": float((np.arctan2(c2, c1) / _W_ANN) % 365.0),
        "t_diurnal_amp": float(np.hypot(c3, c4)),
        "t_diurnal_peak_hour": float((np.arctan2(c4, c3) / _W_DIU) % 24.0),
        "n_hours": int(finite_mask.sum()),
        "n_dropped_feb29": n_dropped_feb29,
        "rmse_k": float(np.sqrt(np.mean(resid**2))),
    }

    if all_rh:
        rh = np.concatenate(all_rh).astype(float)
        rh = rh[np.isfinite(rh)]
        if rh.size:
            result["rh_mean_pct"] = float(np.mean(rh))
            result["rh_median_pct"] = float(np.median(rh))

    return result


def write_namelist(values: dict, relative_humidity: float, output: Path, source_desc: str) -> None:
    """Writes a &fates_test_site namelist file (FatesTestSiteMod.F90's
    ReadSiteNamelist format) from fit_site's results plus a chosen
    relative_humidity.

    Args:
        values (dict): fit_site's return value
        relative_humidity (float): chosen relative_humidity [0-1] - not fit,
            see this module's header comment
        output (Path): path to write the namelist file to
        source_desc (str): human-readable description of the input data,
            recorded in the file's header comment for provenance
    """
    nml_values = {**values, "relative_humidity": relative_humidity}
    lines = [
        "! Generated by fit_site_namelist.py"
        f" on {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}",
        f"! Source data: {source_desc}",
        f"! Fit used {values['n_hours']} hourly TBOT records"
        + (f" ({values['n_dropped_feb29']} Feb-29 records dropped)"
           if values["n_dropped_feb29"] else "")
        + f", residual RMSE {values['rmse_k']:.4f} K",
    ]
    if "rh_mean_pct" in values:
        lines.append(
            f"! Input data's actual RH: mean {values['rh_mean_pct']:.1f}%, "
            f"median {values['rh_median_pct']:.1f}% (relative_humidity below "
            "was NOT fit to this - see this script's header comment)"
        )
    lines.append("&fates_test_site")
    for var in _NAMELIST_VAR_ORDER:
        lines.append(f"  {var:<19s} = {nml_values[var]!r}")
    lines.append("/")

    Path(output).write_text("\n".join(lines) + "\n")


def main():
    """Parses arguments, fits the annual/diurnal TBOT-cycle namelist values
    from --datm-dir, and writes the resulting &fates_test_site namelist file
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--datm-dir", type=Path, required=True,
        help="Directory containing the site's DATM-format hourly netCDF "
             "files (e.g. one per month, matching BCI_datm's layout)",
    )
    parser.add_argument(
        "--relative-humidity", type=float, required=True,
        help="Fixed canopy-air relative humidity to write to the namelist "
             "[0-1]. NOT fit from data - a deliberate choice (see this "
             "script's header comment). If --datm-dir has an RH variable, "
             "its mean/median are printed as a reference point",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("site_namelist.nml"),
        help="Path to write the &fates_test_site namelist file to "
             "(default: site_namelist.nml in the current directory)",
    )
    parser.add_argument(
        "--pattern", default="*.nc",
        help="Glob pattern for files within --datm-dir (default: *.nc)",
    )
    parser.add_argument("--tbot-var", default="TBOT", help="Hourly temperature variable name [K] (default: TBOT)")
    parser.add_argument("--rh-var", default="RH", help="Hourly relative humidity variable name [%%], for reference only (default: RH)")
    parser.add_argument("--lat-var", default="LATIXY", help="Latitude coordinate name (default: LATIXY)")
    parser.add_argument("--lon-var", default="LONGXY", help="Longitude coordinate name (default: LONGXY)")
    args = parser.parse_args()

    if not 0.0 < args.relative_humidity <= 1.0:
        raise ValueError(
            f"--relative-humidity must be in (0, 1] (got {args.relative_humidity}) "
            "- it's a fraction, e.g. 0.70 for 70% RH, not a percentage"
        )

    values = fit_site(
        args.datm_dir, pattern=args.pattern, tbot_var=args.tbot_var,
        rh_var=args.rh_var, lat_var=args.lat_var, lon_var=args.lon_var,
    )

    print(f"Fit from {values['n_hours']} hourly TBOT records "
          f"(RMSE {values['rmse_k']:.4f} K):")
    for var in _NAMELIST_VAR_ORDER[:-1]:
        print(f"  {var:<19s} = {values[var]}")
    print(f"  relative_humidity   = {args.relative_humidity}  (chosen, not fit)")
    if "rh_mean_pct" in values:
        print(f"  [reference only: input data's actual RH mean/median = "
              f"{values['rh_mean_pct']:.1f}% / {values['rh_median_pct']:.1f}%]")

    write_namelist(
        values, args.relative_humidity, args.output,
        source_desc=f"{args.datm_dir} ({args.pattern})",
    )
    print(f"\nNamelist written to {args.output}")


if __name__ == "__main__":
    main()
