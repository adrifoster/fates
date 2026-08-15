"""
Concrete class for running the single_cohort functional test for FATES.
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from framework.utils.plotting import get_color_palette
from framework.functional_test import FunctionalTest


class SingleCohort(FunctionalTest):
    """Single-cohort light-sweep test class"""

    name = "single_cohort"

    # netCDF layers above a light level/year's actual nv are filled with this value
    _FILL_VALUE = -1.0e36

    # a cohort with less leaf area than this is treated as effectively
    # leafless/dead for compensation-point purposes
    _ALIVE_TREELAI_THRESHOLD = 1.0e-3

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots output from the single-cohort test.

        Args:
            run_dir (str): run directory
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        out_path = os.path.join(run_dir, self.out_file)
        if not os.path.exists(out_path):
            raise FileNotFoundError(
                f"Output file not found: {out_path}"
            )

        light_dat = xr.open_dataset(out_path)

        self.plot_light_response(light_dat, save_figs, plot_dir)
        self.plot_growth_trajectory(light_dat, save_figs, plot_dir)
        self.plot_daily_net_carbon(light_dat, save_figs, plot_dir)
        self.plot_gpp_light_response(light_dat, save_figs, plot_dir)
        self.plot_light_response_curves(light_dat, save_figs, plot_dir)

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
        final_year = data.isel(year=-1)
        total_par = (final_year["parsun_z"] + final_year["parsha_z"]).where(
            final_year["parsun_z"] != cls._FILL_VALUE
        )

        leaf_layers = data.nlevleaf.values[np.any(~np.isnan(total_par), axis=0)]
        cmap = plt.get_cmap("viridis_r")
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
        axis.set_title("Light response by leaf layer (final year)", fontsize=11)

        mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array([])
        cbar = fig.colorbar(mappable, ax=axis, pad=0.02)
        cbar.set_label("Leaf layer (1 = to of canopy)", fontsize=10)

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
        axes

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
    def plot_light_response_curves(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None,
                                    light_level: float = 1.0):
        """Plots net whole-plant assimilation (gross_assim - total_resp)
        against PPFD, one curve per simulated year, for a single growth
        light level

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

