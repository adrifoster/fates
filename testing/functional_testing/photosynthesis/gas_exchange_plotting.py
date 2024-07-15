"""Utility functions for gas exchange functional tests
"""
import os
import math
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def plot_gas_exchange_dat(run_dir, out_file, save_figs, plot_dir):
    """Reads in and plots gas exchange test output

    Args:
        run_dir (str): run directory
        out_file (str): output file
        save_figs (bool): whether or not to save the figures
        plot_dir (str): plot directory
    """
    
    ## observations from Bernacchi et al. 2001
    kc_vals = [40.49]
    ko_vals = [27840.0]
    cpt_vals = [4.275]
    temp_vals = np.repeat(25 + 273.15, len(kc_vals))

    # read in quadratic data
    gas_exchange_dat = xr.open_dataset(os.path.join(run_dir, out_file))
    
    plot_dict = {
      'mm_co2': {
        'varname': 'Michaelis-Menten Constant for CO$_2$',
        'units': 'Pa',
        'temp_obs': temp_vals,
        'var_obs': kc_vals,
        },
      'mm_o2': {
        'varname': 'Michaelis-Menten Constant for O$_2$',
        'units': 'Pa',
        'temp_obs': temp_vals,
        'var_obs': ko_vals,
        },
      'co2_comp_pt': {
        'varname': 'CO$_2$ Compensation Point',
        'units': 'Pa',
        'temp_obs': temp_vals,
        'var_obs': cpt_vals,
        },
      'cf': {
        'varname': 'Conversion factor between molar form and velocity form of conductance and resistance',
        'units': 'umol m$^{-3}$',
        'temp_obs': None,
        'var_obs': None,
        },
      'bl_conductance': {
        'varname': 'Leaf Boundary Layer Conductance',
        'units': 'umol H$_2$O m$^{-2}$ s$^{-1}$',
        'temp_obs': None,
        'var_obs': None,
        },
      'air_vpress_constrained': {
        'varname': 'Vapor Pressure of Air, Constrained',
        'units': 'Pa',
        'temp_obs': None,
        'var_obs': None,
        },
    }
    
    # plot output
    for plot, attributes in plot_dict.items():
        plot_gas_var(gas_exchange_dat, plot, attributes['temp_obs'],
                     attributes['var_obs'], 
                     f"{attributes['varname']} ({attributes['units']})", plot_dir, save_figs)
    
def plot_gas_var(dat, var, temp_obs, var_obs, var_label, plot_dir, save_figs):
    """Plots and individual output variable from the gas exchange test

    Args:
        dat (xarray DataSet): the output DataSet from the test
        var (str): variable name
        temp_obs (float): temperature observation
        var_obs (float): variable observation
        var_label (str): label to use for variable (include units)
        plot_dir (str): plotting directory
        save_figs (bool): whether or not to save figures
    """
    
    # color palette:
    colors = ["#2B2F42", "#8D99AE", "#EF233C"]
    
    fig, ax = plt.subplots(figsize=(6, 5))
    dat[var].plot(c=colors[0], linewidth=2, label='FATES')
    
    if temp_obs is not None:
        ax.scatter(temp_obs, var_obs,
            label="Bernacchi et al. (2001)",
            color=colors[2]
        )

    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    ax.spines["bottom"].set_bounds(min(dat.temperat), max(dat.temperat))

    ax.set_xlabel("Vegetation Temperature (K)")
    ax.set_ylabel(var_label)
    ax.legend()
   
    if save_figs:
        fig_name = os.path.join(plot_dir, f"gas_exch_{var}_test.png")
        plt.savefig(fig_name)