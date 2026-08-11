"""
Concrete class for running the canopy_level_photo functional test for FATES.
"""
import os
import xarray as xr
import matplotlib.pyplot as plt
from framework.functional_test import FunctionalTest
from framework.utils.plotting import style_axis


class CanopyLevelPhoto(FunctionalTest):
    """Canopy-level photosynthesis sensitivity sweep test class

    The canopy companion to leaf_level_photo: the same five sweeps (PAR, CO2,
    VPD, leaf temperature, soil water content) over the same ranges, but run at
    each prescribed leaf area index and integrated to a canopy flux per unit
    ground area. Each sweep therefore plots as one line per LAI, so the LAI
    spread at a given x value is the canopy-structure effect the leaf-level
    test cannot show.
    """

    name = "canopy_level_photo"

    # matches FatesConstantsMod's t_water_freeze_k_1atm, used only to convert
    # the temp sweep's K output to degC for plotting
    _T_FREEZE_K = 273.15

    # (swept dimension, x-axis label, plot-title fragment) - variable names in
    # the output follow the <metric>_by<dim> convention the driver writes
    _SWEEPS = [
        ("par", "Incident PPFD at canopy top ($\\mu$mol m$^{-2}$ s$^{-1}$)", "PAR"),
        ("co2", "CO$_2$ partial pressure (Pa)", "CO$_2$"),
        ("vpd", "Vapor pressure deficit (kPa)", "VPD"),
        ("temp", "Leaf temperature ($^{\\circ}$C)", "leaf temperature"),
        ("soilfrac", "Soil water content (fraction of saturation)", "soil water content"),
    ]

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots output from the canopy_level_photo test - one 2-panel figure
        (canopy gross and net photosynthesis, one line per prescribed LAI) per
        sweep, plus the within-canopy profile at the reference condition

        Args:
            run_dir (str): run directory
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        data = xr.open_dataset(os.path.join(run_dir, self.out_file))

        for dim, xlabel, title in self._SWEEPS:
            self.plot_sweep(data, dim, xlabel, title, save_figs, plot_dir)

        self.plot_profile(data, save_figs, plot_dir)

    @classmethod
    def _swept_x(cls, data: xr.Dataset, dim: str):
        """Returns the swept coordinate in the units its axis label states

        Args:
            data (xarray Dataset): the canopy-level photosynthesis dataset
            dim (str): swept dimension name
        """
        x = data[dim].values
        if dim == "temp":
            return x - cls._T_FREEZE_K
        if dim == "vpd":
            return x / 1000.0  # Pa -> kPa
        return x

    @classmethod
    def plot_sweep(cls, data: xr.Dataset, dim: str, xlabel: str, title: str,
                   save_fig: bool, plot_dir: str = None):
        """Plots canopy gross and net photosynthesis against one swept
        variable, one line per prescribed LAI

        Args:
            data (xarray Dataset): the canopy-level photosynthesis dataset
            dim (str): swept dimension name (e.g. "par" for
                canopy_anet_bypar/canopy_agross_bypar)
            xlabel (str): x-axis label
            title (str): title fragment, e.g. "PAR" -> "... vs. PAR"
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        x = cls._swept_x(data, dim)
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
            # net photosynthesis goes negative in the dark and at full stomatal
            # closure, so mark the compensation line rather than leaving the
            # sign change to be read off the tick labels
            if varname.startswith("canopy_anet"):
                axis.axhline(0.0, lw=0.6, color="0.4", linestyle=":")
            style_axis(axis)
            axis.set_xlabel(xlabel, fontsize=10)
            axis.set_ylabel(ylabel, fontsize=10)
            axis.legend(frameon=False, fontsize=9)

        fig.suptitle(f"Canopy-level photosynthesis vs. {title}", fontsize=12)
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, f"canopy_level_photo_{dim}.png"))

    @staticmethod
    def plot_profile(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
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
            style_axis(axis)
            axis.set_xlabel(xlabel, fontsize=10)

        axes[0].set_ylabel("Leaf layer (1 = canopy top)", fontsize=10)
        axes[0].invert_yaxis()
        axes[0].legend(frameon=False, fontsize=9)

        fig.suptitle("Within-canopy profile at the reference condition", fontsize=12)
        fig.tight_layout()

        if save_fig:
            fig.savefig(os.path.join(plot_dir, "canopy_level_photo_profile.png"))
