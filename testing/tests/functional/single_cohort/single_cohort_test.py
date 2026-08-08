"""
Concrete class for running the single_cohort functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from framework.utils.plotting import get_color_palette
from framework.functional_test import FunctionalTest


class SingleCohort(FunctionalTest):
    """Single-cohort light-sweep test class"""

    name = "single_cohort"

    # netCDF layers above a light level/year's actual nv are filled with this
    # sentinel (fates_unset_r8, see FatesConstantsMod.F90) rather than left
    # undefined - mask it out before plotting the per-leaf-layer profile
    _FILL_VALUE = -1.0e36

    # test_SingleCohort.F90's reduced_output mode writes to this filename
    # instead of self.out_file (see its "COMMAND LINE" header comment) - a
    # distinct name, not one shared with the full-output file, specifically
    # so a stale file left over from a manual reduced_output run can never be
    # silently mistaken for the full-output file plot_output requires. Keep
    # this in sync with test_SingleCohort.F90's out_file assignment
    _REDUCED_OUT_FILE = "single_cohort_out_reduced.nc"

    # a cohort with less leaf area than this is treated as effectively
    # leafless/dead for compensation-point purposes (see _state_on_diagnostic_days)
    _ALIVE_TREELAI_THRESHOLD = 1.0e-3

    # _zero_crossing bracket status codes - distinguishes "no crossing because
    # y never changes sign" from "too few finite samples to say", so a
    # degenerate or out-of-sweep-range light-response curve (e.g. a light
    # level extreme enough to kill the cohort, or a compensation point outside
    # the swept PPFD range) can be diagnosed rather than reported as the same
    # silent NaN as every other failure mode
    _BRACKET_OK = 0          # sign change found inside the sweep
    _BRACKET_ALL_POS = 1     # y > 0 everywhere: true crossing is below x.min()
    _BRACKET_ALL_NEG = 2     # y < 0 everywhere: true crossing is above x.max()
    _BRACKET_DEGENERATE = 3  # fewer than 2 finite samples, or a zero-slope bracket
    _BRACKET_LABELS = {
        _BRACKET_OK: "bracketed",
        _BRACKET_ALL_POS: "y > 0 everywhere (crossing below swept range)",
        _BRACKET_ALL_NEG: "y < 0 everywhere (crossing above swept range)",
        _BRACKET_DEGENERATE: "degenerate (too few finite points / zero-slope bracket)",
    }

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots output from the single-cohort test. Requires the full-output
        file (self.out_file) - four of the plots below (plot_light_response,
        plot_daily_net_carbon, plot_gpp_light_response, and the per-leaf-layer
        profile inputs) need variables reduced_output mode omits, so this
        deliberately does not try to run against a reduced-output file.

        Args:
            run_dir (str): run directory
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        out_path = os.path.join(run_dir, self.out_file)
        if not os.path.exists(out_path):
            reduced_path = os.path.join(run_dir, self._REDUCED_OUT_FILE)
            hint = (
                f" A reduced_output file exists instead ({reduced_path}), but plot_output "
                f"needs the full diagnostic suite - re-run without the reduced_output "
                f"command-line argument to regenerate {self.out_file}."
                if os.path.exists(reduced_path)
                else ""
            )
            raise FileNotFoundError(
                f"Expected full-output file not found: {out_path}.{hint}"
            )

        light_dat = xr.open_dataset(out_path)

        self.plot_light_response(light_dat, save_figs, plot_dir)
        self.plot_growth_trajectory(light_dat, save_figs, plot_dir)
        self.plot_daily_net_carbon(light_dat, save_figs, plot_dir)
        self.plot_gpp_light_response(light_dat, save_figs, plot_dir)
        self.plot_compensation_point_trajectory(light_dat, save_figs, plot_dir)
        self.plot_compensation_point_vs_size(light_dat, save_figs, plot_dir)
        self.plot_light_response_curves(light_dat, save_figs, plot_dir)
        self.plot_lie_sterck_trajectory(light_dat, save_figs, plot_dir)
        self.plot_onoda_efficiency(light_dat, save_figs, plot_dir)

    @staticmethod
    def _style_axis(axis):
        """Applies the shared minimalist axis styling used across these plots

        Args:
            axis (matplotlib.axes.Axes): axis to style
        """
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(bottom=False, left=False)
        axis.set_axisbelow(True)
        axis.grid(axis="y", lw=0.5, alpha=0.3, linestyle="--")

    @classmethod
    def _zero_crossing(cls, x, y):
        """Interpolates the first negative-to-positive zero-crossing of
        y(x), linear in x (not log(x)): near a light compensation point, net
        assimilation is close to linear in raw PPFD (the initial slope is
        quantum yield, A ~ phi*I - Rd), not in log(PPFD), so interpolating
        directly in x is the better-motivated choice right at the crossing,
        even though the sweep itself is log-spaced. Falls back to any sign
        change (not just negative-to-positive) if the first pass finds none,
        in case the curve is non-monotonic

        Args:
            x (np.ndarray): sample points, increasing
            y (np.ndarray): y(x) at each sample point, same length as x

        Returns:
            (float, int): interpolated x where y crosses zero (or np.nan),
                and one of the _BRACKET_* status codes explaining why, if not
        """
        x = np.asarray(x)
        y = np.asarray(y)
        finite = np.isfinite(y)
        if finite.sum() < 2:
            return np.nan, cls._BRACKET_DEGENERATE

        x = x[finite]
        y = y[finite]

        if np.all(y > 0):
            return np.nan, cls._BRACKET_ALL_POS
        if np.all(y < 0):
            return np.nan, cls._BRACKET_ALL_NEG

        sign_change = np.where((y[:-1] < 0) & (y[1:] >= 0))[0]
        if sign_change.size == 0:
            sign_change = np.where(np.diff(np.sign(y)) != 0)[0]
            if sign_change.size == 0:
                return np.nan, cls._BRACKET_DEGENERATE

        i = sign_change[0]
        if y[i + 1] == y[i]:
            return np.nan, cls._BRACKET_DEGENERATE

        crossing = x[i] + (0.0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
        return float(crossing), cls._BRACKET_OK

    @classmethod
    def _print_bracket_status(cls, label, status):
        """Prints a one-line-per-outcome summary of how many entries in
        `status` fell into each _BRACKET_* category from _zero_crossing

        Args:
            label (str): name to print, e.g. "LCPplant"
            status (np.ndarray): array of _BRACKET_* codes
        """
        total = status.size
        print(f"{label} bracket status ({total} sweeps):")
        for code, count in zip(*np.unique(status, return_counts=True)):
            print(f"  {cls._BRACKET_LABELS.get(int(code), code):<55s} {count:5d}  ({100 * count / total:.1f}%)")

    @classmethod
    def _compensation_points(cls, data: xr.Dataset):
        """Computes LCPplant and LCPleaf (see _zero_crossing) for every
        (year, light_level) in the light-response sweep, printing a bracket
        status summary for each

        Args:
            data (xarray Dataset): the light environment dataset

        Returns:
            (np.ndarray, np.ndarray): lcp_plant, lcp_leaf, each shaped
                (year, light_level), with unbracketed entries left as np.nan
        """
        years = data.year.values
        light_levels = data.light_level.values
        ppfd = data.ppfd.values
        net_plant = data["gross_assim"] - data["total_resp"]
        leaf_anet = data["leaf_anet"]

        lcp_plant = np.full((len(years), len(light_levels)), np.nan)
        lcp_leaf = np.full((len(years), len(light_levels)), np.nan)
        status_plant = np.zeros((len(years), len(light_levels)), dtype=int)
        status_leaf = np.zeros((len(years), len(light_levels)), dtype=int)

        for iyear, year in enumerate(years):
            for ilight, light_frac in enumerate(light_levels):
                lcp_plant[iyear, ilight], status_plant[iyear, ilight] = cls._zero_crossing(
                    ppfd, net_plant.sel(year=year, light_level=light_frac).values
                )
                lcp_leaf[iyear, ilight], status_leaf[iyear, ilight] = cls._zero_crossing(
                    ppfd, leaf_anet.sel(year=year, light_level=light_frac).values
                )

        cls._print_bracket_status("LCPplant", status_plant)
        cls._print_bracket_status("LCPleaf", status_leaf)

        return lcp_plant, lcp_leaf

    @staticmethod
    def _state_on_diagnostic_days(data: xr.Dataset):
        """Selects daily cohort state on the days the annual light-response
        diagnostics were captured (day 1 of each simulated year, solar noon -
        see test_SingleCohort.F90's LightResponseSweep/LeafNetAssimSweep
        header comments), so a size axis (e.g. total leaf area) can be lined
        up against the (year, light_level)-indexed compensation points from
        _compensation_points. days_per_year is derived from the data rather
        than hardcoded, since it is a compile-time constant on the Fortran side

        Args:
            data (xarray Dataset): the light environment dataset

        Returns:
            xarray Dataset: data's daily variables, indexed by (year, light_level)
                instead of (time, light_level)
        """
        days_per_year = data.sizes["time"] // data.sizes["year"]
        day_idx = (data.year.values - 1) * days_per_year  # 0-based time index of day 1
        return data.isel(time=xr.DataArray(day_idx, dims="year", coords={"year": data.year}))

    @classmethod
    def _growth_floor_edge(cls, data: xr.Dataset):
        """Computes the annual dbh-growth-floor edge for each simulated year:
        the light fraction below which DailyPRTAllometricCarbon's phase-3
        stature-growth integrator is never invoked at all
        (PRTAllometricCarbonMod.F90's if_carbon_increment gate, guarding on
        carbon_balance > calloc_abs_error), so annual dbh increment is
        identically 0.0 rather than merely small - this is the
        tolerance_inversion_plan's original "light fraction at which annual
        dbh increment = 0" metric.

        Unlike LCPplant/LCPleaf, this is not a smooth zero-crossing: it is the
        edge of a hard plateau (a binary gate, not a clamp), so annual dbh
        increment is exactly 0.0 across a whole low-light interval rather than
        approaching zero asymptotically. _zero_crossing still locates the edge
        correctly despite being written for smooth crossings: at the
        plateau's last point y is exactly 0.0, so the interpolation formula
        x[i] + (0-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]) collapses to x[i] exactly
        - the last light fraction still inside the zero-growth plateau, not
        some point between it and the first light fraction with positive
        growth. Interpolated in log(light_frac), matching the sweep's own log
        spacing (there is no physiological argument for one over the other
        here, unlike the PPFD interpolation in _zero_crossing's own
        docstring, since this is a hard threshold rather than a physiological
        curve)

        Args:
            data (xarray Dataset): the light environment dataset

        Returns:
            np.ndarray: growth-floor edge per year, shaped (year,), in units
                of incident light fraction [0-1]. np.nan for a year whose
                entire swept light-fraction range fell on one side of the
                floor (see _zero_crossing's ALL_POS/ALL_NEG status codes)
        """
        light_frac = data.light_level.values
        log_light_frac = np.log(light_frac)
        years = data.year.values
        days_per_year = data.sizes["time"] // data.sizes["year"]

        dbh = data["dbh"]
        edge = np.full(len(years), np.nan)
        status = np.zeros(len(years), dtype=int)
        for iyear, year in enumerate(years):
            end_idx = int(year * days_per_year - 1)
            start_idx = max(int((year - 1) * days_per_year - 1), 0)
            annual_increment = (dbh.isel(time=end_idx) - dbh.isel(time=start_idx)).values
            crossing, status[iyear] = cls._zero_crossing(log_light_frac, annual_increment)
            if status[iyear] == cls._BRACKET_OK:
                edge[iyear] = np.exp(crossing)

        cls._print_bracket_status("Growth-floor edge", status)
        return edge

    # ==========================================================================
    # Shade-tolerance metric calculations - public API, meant to be called
    # directly (e.g. from an eventual Morris/LHC post-processing script) as
    # well as by the plotting methods below, which build on these rather than
    # recomputing the same math inline.
    # ==========================================================================

    @classmethod
    def calculate_compensation_points(cls, data: xr.Dataset) -> xr.Dataset:
        """Whole-plant (LCPplant) and leaf-level (LCPleaf) light compensation
        points, per (year, light_level) - see _compensation_points/
        _zero_crossing for the underlying computation. Entries where the
        cohort had already died back (treelai <= _ALIVE_TREELAI_THRESHOLD on
        that year's diagnostic day - see _state_on_diagnostic_days) are set
        to NaN in addition to any _zero_crossing bracket-status failures, so
        a downstream consumer (e.g. an LHC/Morris post-processing script)
        does not have to separately reproduce this masking to avoid feeding a
        post-mortality "phantom" compensation point into a GP emulator or
        history-matching step

        Args:
            data (xarray Dataset): the light environment dataset (reduced_output
                or full - both carry everything this needs)

        Returns:
            xarray Dataset with two (year, light_level) variables: lcp_plant,
                lcp_leaf, both umol m-2 s-1, NaN where dead or unbracketed
        """
        lcp_plant, lcp_leaf = cls._compensation_points(data)

        alive = (
            cls._state_on_diagnostic_days(data)["treelai"] > cls._ALIVE_TREELAI_THRESHOLD
        ).transpose("year", "light_level").values
        lcp_plant = np.where(alive, lcp_plant, np.nan)
        lcp_leaf = np.where(alive, lcp_leaf, np.nan)

        coords = {"year": data.year.values, "light_level": data.light_level.values}
        return xr.Dataset(
            {
                "lcp_plant": (
                    ("year", "light_level"), lcp_plant,
                    {"units": "umol m-2 s-1",
                     "long_name": "whole-plant light compensation point (Sterck et al. 2013)"},
                ),
                "lcp_leaf": (
                    ("year", "light_level"), lcp_leaf,
                    {"units": "umol m-2 s-1",
                     "long_name": "leaf-level light compensation point (Sterck et al. 2013)"},
                ),
            },
            coords=coords,
        )

    @classmethod
    def calculate_growth_floor_edge(cls, data: xr.Dataset) -> xr.DataArray:
        """Annual dbh-growth-floor edge (see _growth_floor_edge) - the
        tolerance_inversion_plan's original "annual dbh increment = 0" metric,
        wrapped as a labeled DataArray for use outside the plotting code

        Args:
            data (xarray Dataset): the light environment dataset

        Returns:
            xarray DataArray, dims (year,), units fraction of full sun. NaN
                for a year whose entire swept light-fraction range fell on
                one side of the floor
        """
        edge = cls._growth_floor_edge(data)
        return xr.DataArray(
            edge, dims=("year",), coords={"year": data.year.values},
            name="growth_floor_edge",
            attrs={
                "units": "fraction of full sun",
                "long_name": "annual dbh-growth-floor edge (tolerance_inversion_plan's original 'dbh increment = 0' metric)",
            },
        )

    @staticmethod
    def calculate_lie(data: xr.Dataset) -> xr.DataArray:
        """Sterck et al. (2013) light interception efficiency - a direct
        passthrough of the driver's own light_intercept_eff output (already
        computed in Fortran, see test_SingleCohort.F90); included here purely
        for a uniform calculate_* metrics API, not because any computation
        happens on the Python side

        Args:
            data (xarray Dataset): the light environment dataset

        Returns:
            xarray DataArray, dims (time, light_level), dimensionless [0-1]
        """
        return data["light_intercept_eff"]

    @staticmethod
    def _year_time_slice(data: xr.Dataset, year: int):
        """0-based time-index slice spanning the given simulated year's own
        days (1-indexed, matching data.year) - used by the Onoda metrics
        below to average over exactly one year without hardcoding
        days_per_year (a compile-time constant on the Fortran side)

        Args:
            data (xarray Dataset): the light environment dataset
            year (int): 1-indexed simulated year

        Returns:
            slice: 0-based time-index slice for that year's days
        """
        days_per_year = data.sizes["time"] // data.sizes["year"]
        return slice((year - 1) * days_per_year, year * days_per_year)

    @classmethod
    def calculate_lie_la(cls, data: xr.Dataset, year: int = None) -> xr.DataArray:
        """Onoda et al. (2013)'s LIE_LA (Phi/LA): light interception
        efficiency per unit leaf area, averaged over the requested simulated
        year (default: the final year) at each swept light level. light_frac
        stands in for a tree's competitive light environment/social position,
        in place of Onoda's original height gradient across co-occurring
        trees - see test_SingleCohort.F90's header comment

        Args:
            data (xarray Dataset): the light environment dataset
            year (int): 1-indexed simulated year to average over (default:
                the final simulated year)

        Returns:
            xarray DataArray, dims (light_level,), units J m-2 leaf day-1
        """
        if year is None:
            year = int(data.year.values[-1])
        time_slice = cls._year_time_slice(data, year)

        phi = data["daily_absorbed_par_indiv"]
        leaf_area = data["treelai"] * data["crown_area"] / data["n"]
        lie_la = (phi.isel(time=time_slice) / leaf_area.isel(time=time_slice)).mean(dim="time")

        lie_la.name = "lie_la"
        lie_la.attrs = {
            "units": "J m-2 leaf day-1",
            "long_name": f"light interception efficiency per unit leaf area (Onoda et al. 2013), year {year} mean",
        }
        return lie_la

    @classmethod
    def calculate_lie_m(cls, data: xr.Dataset, year: int = None) -> xr.DataArray:
        """Onoda et al. (2013)'s LIE_M (Phi/M): light interception efficiency
        per unit total biomass, averaged over the requested simulated year
        (default: the final year) at each swept light level. M sums all five
        carbon pools (leaf/fine root/sapwood/structure/storage) - FATES has
        no clean above-ground-only split for struct_c the way Onoda's M
        excludes fine roots, so this is a total-biomass analogue, not a
        literal reproduction of their above-ground-only M

        Args:
            data (xarray Dataset): the light environment dataset
            year (int): 1-indexed simulated year to average over (default:
                the final simulated year)

        Returns:
            xarray DataArray, dims (light_level,), units J kgC-1 day-1
        """
        if year is None:
            year = int(data.year.values[-1])
        time_slice = cls._year_time_slice(data, year)

        phi = data["daily_absorbed_par_indiv"]
        biomass = data["leaf_c"] + data["fnrt_c"] + data["sapw_c"] + data["struct_c"] + data["storage_c"]
        lie_m = (phi.isel(time=time_slice) / biomass.isel(time=time_slice)).mean(dim="time")

        lie_m.name = "lie_m"
        lie_m.attrs = {
            "units": "J kgC-1 day-1",
            "long_name": f"light interception efficiency per unit biomass (Onoda et al. 2013), year {year} mean",
        }
        return lie_m

    @classmethod
    def calculate_lue(cls, data: xr.Dataset, year: int = None) -> xr.DataArray:
        """Onoda et al. (2013)'s LUE: biomass increment over the requested
        simulated year (default: the final year) divided by that year's
        cumulative absorbed PAR, at each swept light level. Can be negative
        (net biomass loss under carbon starvation at the lowest light
        levels) - see calculate_lie_m for the biomass-pool caveat, which
        applies here too. Unlike calculate_lie_la/_lie_m (an average within
        one year), this differences across the year boundary (end of year
        `year` minus end of year `year - 1`), matching how a biomass
        increment must be measured

        Args:
            data (xarray Dataset): the light environment dataset
            year (int): 1-indexed simulated year (default: the final
                simulated year)

        Returns:
            xarray DataArray, dims (light_level,), units kgC J-1. NaN where
                that year's cumulative absorbed PAR was zero
        """
        if year is None:
            year = int(data.year.values[-1])
        days_per_year = data.sizes["time"] // data.sizes["year"]
        end_idx = year * days_per_year - 1
        start_idx = max((year - 1) * days_per_year - 1, 0)

        biomass = data["leaf_c"] + data["fnrt_c"] + data["sapw_c"] + data["struct_c"] + data["storage_c"]
        phi = data["daily_absorbed_par_indiv"]

        delta_biomass = biomass.isel(time=end_idx) - biomass.isel(time=start_idx)
        sum_phi = phi.isel(time=slice(start_idx, end_idx)).sum(dim="time")
        lue = xr.where(sum_phi > 0, delta_biomass / sum_phi, np.nan)

        lue.name = "lue"
        lue.attrs = {
            "units": "kgC J-1",
            "long_name": f"light use efficiency (Onoda et al. 2013), year {year}",
        }
        return lue

    @classmethod
    def calculate_all_metrics(cls, data: xr.Dataset, year: int = None) -> xr.Dataset:
        """Bundles every shade-tolerance metric this class computes into one
        small Dataset: LCPplant/LCPleaf (year, light_level), the growth-floor
        edge (year,), and LIE_LA/LIE_M/LUE (light_level,) at the requested
        year (default: the final simulated year). Meant to be the actual
        per-draw artifact an eventual Morris/LHC parameter sweep saves out,
        rather than keeping every draw's full (even reduced_output)
        simulation netCDF around - this Dataset is a few KB; the
        reduced_output simulation file itself is several MB. Sterck's LIE
        (light_intercept_eff) is intentionally left out of this bundle: it
        is already a small (time, light_level) array in the source file
        itself (see calculate_lie), and does not reduce to a single scalar
        or per-year value the way the others do without an arbitrary choice
        this function does not want to make silently

        Args:
            data (xarray Dataset): the light environment dataset
            year (int): 1-indexed simulated year for LIE_LA/LIE_M/LUE
                (default: None, the final simulated year - see
                calculate_lie_la/_lie_m/_lue). LCPplant/LCPleaf/the
                growth-floor edge are unaffected by this - they already
                report every simulated year, not just one

        Returns:
            xarray Dataset with lcp_plant/lcp_leaf (year, light_level),
                growth_floor_edge (year,), and lie_la/lie_m/lue
                (light_level,) - the last three from the requested year
        """
        metrics = cls.calculate_compensation_points(data)
        metrics = metrics.assign(
            growth_floor_edge=cls.calculate_growth_floor_edge(data),
            lie_la=cls.calculate_lie_la(data, year=year),
            lie_m=cls.calculate_lie_m(data, year=year),
            lue=cls.calculate_lue(data, year=year),
        )
        return metrics

    @classmethod
    def plot_light_response(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots total absorbed PAR (year 1, solar noon) against incident light
        fraction, one line per leaf layer, colored by layer depth (nlevleaf can
        run past 20, beyond get_color_palette's cap, so this uses a continuous
        colormap + colorbar instead of a discrete per-layer legend). Layers above
        a given light level's nv are filled with _FILL_VALUE (unoccupied that
        year) and masked out here

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        year1 = data.isel(year=0)
        total_par = (year1["parsun_z"] + year1["parsha_z"]).where(
            year1["parsun_z"] > cls._FILL_VALUE
        )

        leaf_layers = data.nlevleaf.values
        cmap = plt.get_cmap("viridis")
        norm = mcolors.Normalize(vmin=leaf_layers.min(), vmax=leaf_layers.max())

        fig, axis = plt.subplots(figsize=(7, 5))
        for layer in leaf_layers:
            layer_par = total_par.sel(nlevleaf=layer)
            if layer_par.isnull().all():
                continue
            axis.plot(
                data.light_level.values,
                layer_par.values,
                lw=1.5,
                color=cmap(norm(layer)),
            )

        SingleCohort._style_axis(axis)
        axis.set_xscale("log")
        axis.set_ylim(bottom=0.0)
        axis.set_xlabel("Incident light (fraction of full sun)", fontsize=11)
        axis.set_ylabel("Absorbed PAR at solar noon (W m$^{-2}$ ground)", fontsize=11)
        axis.set_title("Light response by leaf layer (year 1)", fontsize=11)

        mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array([])
        cbar = fig.colorbar(mappable, ax=axis, pad=0.02)
        cbar.set_label("Leaf layer (1 = top of crown)", fontsize=10)

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_light_response.png"))

    @staticmethod
    def plot_growth_trajectory(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots dbh through time, one line per swept light level, colored by
        incident light fraction

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_levels = data.light_level.values

        cmap = plt.get_cmap("viridis")
        norm = mcolors.LogNorm(vmin=light_levels.min(), vmax=light_levels.max())

        fig, axis = plt.subplots(figsize=(7, 5))
        for light_frac in light_levels:
            axis.plot(
                data.time.values,
                data["dbh"].sel(light_level=light_frac).values,
                lw=1.5,
                color=cmap(norm(light_frac)),
            )

        SingleCohort._style_axis(axis)
        axis.set_xlim(left=0.0)
        axis.set_ylim(bottom=0.0)
        axis.set_xlabel("Day", fontsize=11)
        axis.set_ylabel("dbh (cm)", fontsize=11)
        axis.set_title("Growth trajectory across the light sweep", fontsize=11)

        mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array([])
        cbar = fig.colorbar(mappable, ax=axis, pad=0.02)
        cbar.set_label("Incident light (fraction of full sun)", fontsize=10)

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_growth_trajectory.png"))

    @staticmethod
    def plot_daily_net_carbon(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots final-day daily net carbon (GPP - leaf dark resp - nonleaf MR)
        against incident light fraction - the whole-plant light-response curve,
        showing where the light compensation point (net carbon = 0) falls

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_frac = data.light_level.values
        net_c = data["daily_net_c"].isel(time=-1).values

        colors = get_color_palette(4)

        fig, axis = plt.subplots(figsize=(7, 5))
        axis.fill_between(
            light_frac, net_c, 0.0, where=(net_c >= 0), color=colors[0],
            alpha=0.15, interpolate=True,
        )
        axis.fill_between(
            light_frac, net_c, 0.0, where=(net_c < 0), color=colors[2],
            alpha=0.15, interpolate=True,
        )
        axis.axhline(0.0, color="0.4", lw=1, linestyle="--")
        axis.plot(light_frac, net_c, lw=2, marker="o", markersize=4, color=colors[0])

        SingleCohort._style_axis(axis)
        axis.set_xscale("log")
        axis.set_xlabel("Incident light (fraction of full sun)", fontsize=11)
        axis.set_ylabel("Daily net carbon (kgC indiv$^{-1}$ day$^{-1}$)", fontsize=11)
        axis.set_title(
            "Daily net carbon balance, final day (light compensation point)",
            fontsize=11,
        )
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_daily_net_carbon.png"))

    @staticmethod
    def plot_gpp_light_response(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots final-day daily GPP against incident light fraction on log-log
        axes - a cheap sanity check: in the light-limited regime GPP should scale
        roughly linearly with absorbed light (slope ~1 on log-log axes), so a
        dashed slope-1 reference line (anchored at the dimmest level) makes
        deviations from that easy to spot at a glance

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_frac = data.light_level.values
        gpp = data["daily_gpp"].isel(time=-1).values
        reference = gpp[0] * (light_frac / light_frac[0])

        colors = get_color_palette(2)

        fig, axis = plt.subplots(figsize=(7, 5))
        axis.plot(
            light_frac, reference, lw=1, linestyle="--", color="0.4",
            label="Linear in light (slope 1)",
        )
        axis.plot(
            light_frac, gpp, lw=2, marker="o", markersize=4, color=colors[0],
            label="GPP",
        )

        SingleCohort._style_axis(axis)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Incident light (fraction of full sun)", fontsize=11)
        axis.set_ylabel("Daily GPP (kgC indiv$^{-1}$ day$^{-1}$)", fontsize=11)
        axis.set_title("GPP light response, final day (log-log)", fontsize=11)
        axis.legend(loc="upper left", frameon=False)
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_gpp_log_log.png"))

    @classmethod
    def plot_compensation_point_trajectory(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots whole-plant (LCPplant) and leaf-level (LCPleaf) light
        compensation points, and the annual dbh-growth-floor edge, against
        year. LCPplant/LCPleaf are one line per swept light level, colored by
        incident light fraction: LCPplant is the PPFD where
        gross_assim - total_resp crosses zero (net of self-shading and
        non-leaf maintenance respiration); LCPleaf is the PPFD where leaf_anet
        crosses zero (single canopy layer, no self-shading) - see
        test_SingleCohort.F90's LightResponseSweep/LeafNetAssimSweep. The
        growth-floor edge (see _growth_floor_edge) is a single line in
        light_frac units, not PPFD, since it is itself a threshold within the
        light_frac sweep dimension rather than a per-light-level quantity -
        it is the tolerance_inversion_plan's original "annual dbh increment
        = 0" metric, and sits strictly above LCPplant since
        DailyPRTAllometricCarbon's turnover-replacement and storage-refill
        phases (1-2) claim the first share of any positive carbon balance
        before stature growth (phase 3) sees any. Together these three show
        the ontogenetic rise/fall in compensation point reported by Sterck et
        al. (2013), the LCPleaf < LCPplant separation driven by self-shading
        plus non-leaf respiration, and the width of the "treading water" zone
        between LCPplant and the growth-floor edge

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_levels = data.light_level.values
        years = data.year.values

        compensation = cls.calculate_compensation_points(data)
        lcp_plant = compensation["lcp_plant"].values  # (year, light_level), already alive+bracket masked (NaN)
        lcp_leaf = compensation["lcp_leaf"].values
        growth_floor_edge = cls.calculate_growth_floor_edge(data).values

        cmap = plt.get_cmap("viridis")
        norm = mcolors.LogNorm(vmin=light_levels.min(), vmax=light_levels.max())

        fig, (ax_plant, ax_leaf, ax_floor) = plt.subplots(3, 1, figsize=(7, 11), sharex=True)
        for ilight, light_frac in enumerate(light_levels):
            color = cmap(norm(light_frac))
            mask_plant = np.isfinite(lcp_plant[:, ilight])
            ax_plant.plot(
                years[mask_plant], lcp_plant[mask_plant, ilight], lw=1.5, marker="o", markersize=3, color=color,
            )
            mask_leaf = np.isfinite(lcp_leaf[:, ilight])
            ax_leaf.plot(
                years[mask_leaf], lcp_leaf[mask_leaf, ilight], lw=1.5, marker="o", markersize=3, color=color,
            )
        ax_floor.plot(years, growth_floor_edge, lw=1.5, marker="o", markersize=4, color="0.2")

        cls._style_axis(ax_plant)
        cls._style_axis(ax_leaf)
        cls._style_axis(ax_floor)
        ax_plant.set_ylabel("LCP$_{plant}$ ($\\mu$mol m$^{-2}$ s$^{-1}$)", fontsize=11)
        ax_leaf.set_ylabel("LCP$_{leaf}$ ($\\mu$mol m$^{-2}$ s$^{-1}$)", fontsize=11)
        ax_floor.set_ylabel("Growth-floor edge\n(fraction of full sun)", fontsize=11)
        ax_floor.set_yscale("log")
        ax_floor.set_xlabel("Year", fontsize=11)
        ax_plant.set_title(
            "Whole-plant vs leaf-level LCP (Sterck et al. 2013) and the\n"
            "dbh-growth-floor edge (tolerance_inversion_plan's original metric)",
            fontsize=11,
        )

        mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array([])
        cbar = fig.colorbar(mappable, ax=[ax_plant, ax_leaf, ax_floor], pad=0.02)
        cbar.set_label("Incident light (fraction of full sun)", fontsize=10)

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_compensation_point_trajectory.png"))

    @classmethod
    def plot_compensation_point_vs_size(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots whole-plant (LCPplant) and leaf-level (LCPleaf) light
        compensation points against total leaf area, one trajectory per
        swept light level, colored by incident light fraction - the Sterck
        et al. (2013) Fig. 2b analogue. Plotting against size rather than
        year (see plot_compensation_point_trajectory) matters here because
        light levels reach a given size at very different rates: two light
        levels compared at the same year are not at the same developmental
        stage, whereas comparing at the same leaf area is what Sterck's own
        figure does and is the more faithful reproduction of their point

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_levels = data.light_level.values

        compensation = cls.calculate_compensation_points(data)
        lcp_plant = compensation["lcp_plant"].values  # (year, light_level), already alive+bracket masked (NaN)
        lcp_leaf = compensation["lcp_leaf"].values

        state = cls._state_on_diagnostic_days(data)
        # per-individual, matching Sterck's per-seedling leaf area - not just
        # treelai*crown_area, since cohort%n can drop below 1 via carbon-
        # starvation mortality over the run (see plot_onoda_efficiency's
        # leaf_area, computed the same way)
        total_leaf_area = (state["treelai"] * state["crown_area"] / state["n"]).transpose(
            "year", "light_level"
        ).values

        cmap = plt.get_cmap("viridis")
        norm = mcolors.LogNorm(vmin=light_levels.min(), vmax=light_levels.max())

        fig, (ax_plant, ax_leaf) = plt.subplots(1, 2, figsize=(12, 5))
        for ilight, light_frac in enumerate(light_levels):
            color = cmap(norm(light_frac))
            mask_plant = np.isfinite(lcp_plant[:, ilight])
            if mask_plant.sum() >= 2:
                ax_plant.plot(
                    total_leaf_area[mask_plant, ilight], lcp_plant[mask_plant, ilight],
                    lw=1.5, marker="o", markersize=3, color=color,
                )
            mask_leaf = np.isfinite(lcp_leaf[:, ilight])
            if mask_leaf.sum() >= 2:
                ax_leaf.plot(
                    total_leaf_area[mask_leaf, ilight], lcp_leaf[mask_leaf, ilight],
                    lw=1.5, marker="o", markersize=3, color=color,
                )

        cls._style_axis(ax_plant)
        cls._style_axis(ax_leaf)
        ax_plant.set_xscale("log")
        ax_leaf.set_xscale("log")
        ax_plant.set_xlabel("Total leaf area (m$^2$)", fontsize=11)
        ax_leaf.set_xlabel("Total leaf area (m$^2$)", fontsize=11)
        ax_plant.set_ylabel("LCP$_{plant}$ ($\\mu$mol m$^{-2}$ s$^{-1}$)", fontsize=11)
        ax_leaf.set_ylabel("LCP$_{leaf}$ ($\\mu$mol m$^{-2}$ s$^{-1}$)", fontsize=11)
        fig.suptitle(
            "Light compensation point vs. plant size (Sterck et al. 2013, Fig. 2b analogue)",
            fontsize=11,
        )

        mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array([])
        cbar = fig.colorbar(mappable, ax=[ax_plant, ax_leaf], pad=0.02)
        cbar.set_label("Incident light (fraction of full sun)", fontsize=10)

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_compensation_point_vs_size.png"))

    @classmethod
    def plot_light_response_curves(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None,
                                    light_level: float = 1.0):
        """Plots net whole-plant assimilation (gross_assim - total_resp)
        against swept PPFD, one curve per simulated year, for a single growth
        light level - the ontogenetic shift of the zero crossing is the
        LCPplant change with age visualized directly, complementing
        plot_compensation_point_trajectory/_vs_size's extracted scalar with
        the underlying curve shape (useful for spotting non-monotonic or
        otherwise degenerate curves that a bracket-status code alone can be
        harder to interpret at a glance)

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
            light_level (float): growth light level to select (nearest match)
        """
        sel = data.sel(light_level=light_level, method="nearest")
        net_plant = sel["gross_assim"] - sel["total_resp"]
        years = sel.year.values
        ppfd = sel.ppfd.values

        cmap = plt.get_cmap("plasma")

        fig, axis = plt.subplots(figsize=(7, 5))
        axis.axhline(0.0, color="0.4", lw=1, linestyle="--")
        for iyear, year in enumerate(years):
            axis.plot(
                ppfd, net_plant.sel(year=year).values,
                lw=1.2, color=cmap(iyear / max(len(years) - 1, 1)),
                label=f"year {year}" if year in (years[0], years[-1]) else None,
            )

        cls._style_axis(axis)
        axis.set_xscale("log")
        axis.set_xlabel("PPFD ($\\mu$mol m$^{-2}$ s$^{-1}$)", fontsize=11)
        axis.set_ylabel("Net whole-plant assimilation (kgC indiv$^{-1}$ s$^{-1}$)", fontsize=11)
        axis.set_title(
            f"Whole-plant light-response curves, growth light = {float(sel.light_level):.3g}",
            fontsize=11,
        )
        axis.legend(loc="lower right", frameon=False, fontsize=9)
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_light_response_curves.png"))

    @classmethod
    def plot_lie_sterck_trajectory(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots light interception efficiency (Sterck et al. 2013) against
        day, one line per swept light level, colored by incident light
        fraction - shows the self-shading efficiency declining as each
        trajectory's cohort grows, and its ceiling approaching 1.0 (not the
        PFT's bare PAR absorptance) for small, minimally self-shaded cohorts

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_levels = data.light_level.values

        cmap = plt.get_cmap("viridis")
        norm = mcolors.LogNorm(vmin=light_levels.min(), vmax=light_levels.max())

        fig, axis = plt.subplots(figsize=(7, 5))
        for light_frac in light_levels:
            axis.plot(
                data.time.values,
                data["light_intercept_eff"].sel(light_level=light_frac).values,
                lw=1.5,
                color=cmap(norm(light_frac)),
            )

        cls._style_axis(axis)
        axis.set_ylim(0.0, 1.0)
        axis.set_xlim(left=0.0)
        axis.set_xlabel("Day", fontsize=11)
        axis.set_ylabel("Light interception efficiency (-)", fontsize=11)
        axis.set_title("Self-shading efficiency through ontogeny (Sterck et al. 2013)", fontsize=11)

        mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array([])
        cbar = fig.colorbar(mappable, ax=axis, pad=0.02)
        cbar.set_label("Incident light (fraction of full sun)", fontsize=10)

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_lie_trajectory.png"))

    @classmethod
    def plot_onoda_efficiency(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots light interception efficiency per unit leaf area (LIE_LA,
        averaged over the final simulated year), per unit total biomass
        (LIE_M, same averaging) and light use efficiency (LUE, the final
        year's total biomass increment divided by that year's cumulative
        absorbed PAR) against incident light fraction (Onoda et al. 2013).
        light_frac stands in here for a tree's competitive light environment/
        social position, in place of Onoda's original height gradient across
        co-occurring trees - see test_SingleCohort.F90's header comment.
        M (total biomass) sums all five carbon pools - FATES has no clean
        above-ground-only split for struct_c the way Onoda's M excludes fine
        roots, so this is a total-biomass analogue, not a literal
        reproduction of their above-ground-only M. LUE can be negative
        (net biomass loss under carbon starvation at the lowest light
        levels), so unlike LIE_LA/LIE_M it is not log-scaled

        Args:
            data (xarray Dataset): the light environment dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        light_frac = data.light_level.values

        lie_la = cls.calculate_lie_la(data).values
        lie_m = cls.calculate_lie_m(data).values
        lue = cls.calculate_lue(data).values

        colors = get_color_palette(3)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        axes[0].plot(light_frac, lie_la, lw=2, marker="o", markersize=4, color=colors[0])
        axes[0].set_xscale("log")
        axes[0].set_yscale("log")
        axes[0].set_ylabel("LIE$_{LA}$ (J m$^{-2}$ leaf day$^{-1}$)", fontsize=11)
        axes[0].set_title("Light interception\nper unit leaf area", fontsize=11)

        axes[1].plot(light_frac, lie_m, lw=2, marker="o", markersize=4, color=colors[1])
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_ylabel("LIE$_{M}$ (J kgC$^{-1}$ day$^{-1}$)", fontsize=11)
        axes[1].set_title("Light interception\nper unit biomass", fontsize=11)

        axes[2].axhline(0.0, color="0.4", lw=1, linestyle="--")
        axes[2].fill_between(
            light_frac, lue, 0.0, where=(lue >= 0), color=colors[2], alpha=0.15, interpolate=True,
        )
        axes[2].fill_between(
            light_frac, lue, 0.0, where=(lue < 0), color=colors[2], alpha=0.15, interpolate=True,
        )
        axes[2].plot(light_frac, lue, lw=2, marker="o", markersize=4, color=colors[2])
        axes[2].set_xscale("log")
        axes[2].set_ylabel("LUE (kgC J$^{-1}$)", fontsize=11)
        axes[2].set_title("Light use efficiency", fontsize=11)

        for axis in axes:
            cls._style_axis(axis)
            axis.set_xlabel("Incident light (fraction of full sun)", fontsize=11)

        fig.suptitle(
            "Light interception & use efficiency, final year (Onoda et al. 2013)",
            fontsize=11,
        )
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "single_cohort_onoda_efficiency.png"))
