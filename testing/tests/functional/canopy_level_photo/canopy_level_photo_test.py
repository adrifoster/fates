"""
Concrete class for running a CanopyLevelPhoto functional tests for FATES.
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from framework.functional_test import FunctionalTest


class CanopyLevelPhoto(FunctionalTest):
    """CanopyLevelPhoto test class
    """
    name = "canopy_level_photo"

    # matches FatesConstantsMod's t_water_freeze_k_1atm, used only to convert
    # the temp sweep's K output to degC for plotting
    _T_FREEZE_K = 273.15

    # matches FatesTestEnvironmentMod's sea_level_press [Pa]
    _CAN_PRESS_PA = 101325.0

    # matches FatesTestEnvironmentMod's default_co2_molfrac [ppm]
    _DEFAULT_CO2_PPM = 380.0

    # (swept dimension, x-axis label, plot-title fragment) for each of the
    # sweeps
    _SWEEPS = [
        ("par", "Incident PPFD at canopy top ($\\mu$mol m$^{-2}$ s$^{-1}$)", "PAR"),
        ("co2", "CO$_2$ (ppm)", "CO$_2$"),
        ("vpd", "Vapor pressure deficit (kPa)", "VPD"),
        ("temp", "Leaf temperature ($^{\\circ}$C)", "leaf temperature"),
        ("soilfrac", "Soil water content (fraction of saturation)", "soil water content"),
    ]

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        data = xr.open_dataset(os.path.join(run_dir, self.out_file))

        for dim, xlabel, title in self._SWEEPS:
            self.plot_sweep(data, dim, xlabel, title, save_figs, plot_dir)

        self.plot_profile(data, save_figs, plot_dir)

                          
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
    def plot_sweep(cls, data: xr.Dataset, dim: str, xlabel: str, title: str,
                    save_fig: bool, plot_dir: str = None):
        """Plots gross/net photosynthesis, stomatal conductance, and
        intracellular CO2 against one swept variable, for the single PFT
        the test was run with

        Args:
            data (xarray Dataset): the leaf-level photosynthesis dataset
            dim (str): swept dimension name (e.g. "par" for
                anet_bypar/agross_bypar/gs_bypar/ci_bypar)
            xlabel (str): x-axis label
            title (str): title fragment, e.g. "PAR" -> "... vs. PAR"
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        x = data[dim].values
        if dim == "temp":
            x = x - cls._T_FREEZE_K
        elif dim == "vpd":
            x = x / 1000.0  # Pa -> kPa
        elif dim == "co2":
            x = x / cls._CAN_PRESS_PA * 1.0e6  # Pa -> ppm
            
        lai_vals = data["lai"].values

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        panels = [
            (axes[0], f"canopy_agross_by{dim}",
             "Canopy gross photosynthesis\n($\\mu$molC m$^{-2}$ ground s$^{-1}$)"),
            (axes[1], f"canopy_anet_by{dim}",
             "Canopy net photosynthesis\n($\\mu$molC m$^{-2}$ ground s$^{-1}$)"),
        ]

        for axis, varname, ylabel in panels:
            for ilai, lai in enumerate(lai_vals):
                axis.plot(x, data[varname].isel(lai=ilai).values, lw=1.2,
                          label=f"LAI = {lai:g}")
            if varname.startswith("canopy_anet"):
                axis.axhline(0.0, lw=0.6, color="0.4", linestyle=":")
            cls._style_axis(axis)
            axis.set_xlabel(xlabel, fontsize=10)
            axis.set_ylabel(ylabel, fontsize=10)
            axis.legend(frameon=False, fontsize=9)

        fig.suptitle(f"Canopy-level photosynthesis vs. {title}", fontsize=12)
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, f"canopy_level_photo_{dim}.png"))


    @classmethod
    def plot_profile(cls, data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plots the within-canopy vertical profile at the reference condition,
        one line per prescribed LAI

        Layers beyond a given LAI's occupied layer count (nv) are _FillValue and
        so drop out of the plot automatically - the driver registers the fill
        value, and xarray masks it on read.

        Args:
            data (xarray Dataset): the canopy-level photosynthesis dataset
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        lai_vals = data["lai"].values
        layer = data["layer"].values

        panels = [
            ("parsun_z", "Absorbed PAR, sunlit\n(W m$^{-2}$ ground)"),
            ("parsha_z", "Absorbed PAR, shaded\n(W m$^{-2}$ ground)"),
            ("nscaler_z", "Nitrogen-scaling factor (-)"),
            ("anet_z", "Net photosynthesis\n($\\mu$molC m$^{-2}$ leaf s$^{-1}$)"),
        ]

        fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), sharey=True)
        for axis, (varname, xlabel) in zip(axes, panels):
            for ilai, lai in enumerate(lai_vals):
                axis.plot(data[varname].isel(lai=ilai).values, layer, lw=1.2,
                          marker="o", markersize=3, label=f"LAI = {lai:g}")
            cls._style_axis(axis)
            axis.set_xlabel(xlabel, fontsize=10)

        axes[0].set_ylabel("Leaf layer (1 = canopy top)", fontsize=10)
        axes[0].invert_yaxis()
        axes[0].legend(frameon=False, fontsize=9)

        fig.suptitle("Within-canopy profile at the reference condition", fontsize=12)
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "canopy_level_photo_profile.png"))


