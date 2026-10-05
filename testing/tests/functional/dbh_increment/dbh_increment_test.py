"""
Concrete class for running a DbhIncrement functional tests for FATES.
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from framework.utils.general import round_up
from framework.utils.plotting import blank_plot, get_color_palette
from framework.functional_test import FunctionalTest


class DbhIncrement(FunctionalTest):
    """DbhIncrement test class
    """
    name = "dbh_increment"

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots dbh increment output

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """

        dbh_incr_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        self.plot_dbh_increment(dbh_incr_dat, save_figs, plot_dir)

    @staticmethod
    def plot_dbh_increment(data: xr.Dataset, save_fig: bool, plot_dir: str):
        """Plot annual dbh increment over time

        Args:
            data (xarray DataSet): the dbh increment dataset
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        max_year = data.year.values.max()
        max_incr = round_up(data.dbh_incr.values.max(), decimals=1)

        blank_plot(max_year, 0.0, max_incr, 0.0, draw_horizontal_lines=False)

        colors = get_color_palette(1)
        plt.plot(data.year.values, data.dbh_incr.values, lw=2, color=colors[0])

        plt.xlabel("Year", fontsize=11)
        plt.ylabel("DBH increment (cm yr$^{-1}$)", fontsize=11)
        plt.title("Simulated annual DBH increment for input parameter file", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, "dbh_increment_plot.png")
            plt.savefig(fig_name)
