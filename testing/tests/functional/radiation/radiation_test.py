"""
Concrete class for running the radiation functional tests for FATES.
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
from framework.utils.general import round_up
from framework.utils.plotting import blank_plot, get_color_palette
from framework.functional_test import FunctionalTest

BAND_NAMES = {1: "Visible", 2: "Near-Infrared"}


class Radiation(FunctionalTest):
    """Radiation test class"""

    name = "radiation"

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots all radiation plots

        Args:
            run_dir (str): run directory
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """

        radiation_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        for band in radiation_dat.band.values:
            self.plot_intensity_profile(radiation_dat, int(band), save_figs, plot_dir)
            self.plot_absorbed_radiation(radiation_dat, int(band), save_figs, plot_dir)
            self.plot_sun_fraction(radiation_dat, int(band), save_figs, plot_dir)

    @staticmethod
    def plot_intensity_profile(data: xr.Dataset, band: int, save_fig: bool, plot_dir: str = None):
        """Plot normalized beam/diffuse radiation intensity vs. depth into the canopy

        Args:
            data (xarray Dataset): the radiation dataset
            band (int): shortwave band to plot (1=visible, 2=near-infrared)
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        band_dat = data.sel(band=band)
        max_vai = round_up(float(data.vai.max()))

        blank_plot(1.0, 0.0, max_vai, 0.0, draw_horizontal_lines=False)

        colors = get_color_palette(3)
        plt.plot(band_dat.r_beam.values, data.vai.values, lw=2, color=colors[0], label="beam")
        plt.plot(band_dat.r_diff_dn.values, data.vai.values, lw=2, color=colors[1], label="diffuse down")
        plt.plot(band_dat.r_diff_up.values, data.vai.values, lw=2, color=colors[2], label="diffuse up")

        plt.gca().invert_yaxis()
        plt.xlabel("Normalized radiation intensity [-]", fontsize=11)
        plt.ylabel("Vegetation area index [m$^2$ m$^{-2}$]", fontsize=11)
        plt.title(f"Radiation intensity profile ({BAND_NAMES[band]})", fontsize=11)
        plt.legend(loc="lower right")

        if save_fig:
            fig_name = os.path.join(plot_dir, f"radiation_intensity_profile_band{band}.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_absorbed_radiation(data: xr.Dataset, band: int, save_fig: bool, plot_dir: str = None):
        """Plot absorbed leaf and stem radiation vs. depth into the canopy

        Args:
            data (xarray Dataset): the radiation dataset
            band (int): shortwave band to plot (1=visible, 2=near-infrared)
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        band_dat = data.sel(band=band)
        total_abs = band_dat.rb_abs_leaf + band_dat.rd_abs_leaf + band_dat.r_abs_stem
        max_abs = round_up(float(total_abs.max()), decimals=3)
        max_vai = round_up(float(data.vai.max()))

        blank_plot(max_abs, 0.0, max_vai, 0.0, draw_horizontal_lines=False)

        plt.plot(total_abs.values, data.vai.values, lw=2, color=get_color_palette(1)[0])

        plt.gca().invert_yaxis()
        plt.xlabel("Absorbed radiation [W m$^{-2}$]", fontsize=11)
        plt.ylabel("Vegetation area index [m$^2$ m$^{-2}$]", fontsize=11)
        plt.title(f"Absorbed radiation profile ({BAND_NAMES[band]})", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, f"radiation_absorbed_profile_band{band}.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_sun_fraction(data: xr.Dataset, band: int, save_fig: bool, plot_dir: str = None):
        """Plot sunlit fraction of leaf area vs. depth into the canopy

        Args:
            data (xarray Dataset): the radiation dataset
            band (int): shortwave band to plot (1=visible, 2=near-infrared)
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        band_dat = data.sel(band=band)
        max_vai = round_up(float(data.vai.max()))

        blank_plot(1.0, 0.0, max_vai, 0.0, draw_horizontal_lines=False)

        plt.plot(band_dat.leaf_sun_frac.values, data.vai.values, lw=2, color=get_color_palette(1)[0])

        plt.gca().invert_yaxis()
        plt.xlabel("Sunlit fraction of leaf area [-]", fontsize=11)
        plt.ylabel("Vegetation area index [m$^2$ m$^{-2}$]", fontsize=11)
        plt.title(f"Sunlit leaf fraction profile ({BAND_NAMES[band]})", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, f"radiation_sun_frac_profile_band{band}.png")
            plt.savefig(fig_name)
