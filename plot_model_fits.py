import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
import sys
from pipes_utils import *

# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

# datafiles = ["phil_model_1", "phil_model_2", "phil_model_3", "phil_model_4",
#              "phil_model_5", "phil_model_6", "phil_model_7", "phil_model_8",
#              "phil_model_9", "phil_model_10"]

datafiles = ["phil_model_4"]

exponential1 = {}
exponential1["age"] = (7.5, 12.5)   # Gyr
exponential1["tau"] = (0.5, 2.0)    # Gyr
exponential1["massformed"] = (5.0, 12.5)   # log_10(M*/M_solar)
exponential1["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

exponential2 = {}
exponential2["age"] = (0.0, 3.5)   # Gyr, lifetime of F type stars
exponential2["tau"] = (0.1, 2.0)    # Gyr
exponential2["massformed"] = (0.0, 12.5)   # log_10(M*/M_solar)
exponential2["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

dblplaw = {}
dblplaw["tau"] = (0.0, 4.5) # Do not let the peak occur beyond ~6.5 Gyr
dblplaw["alpha"] = (5.0, 10.0) # Formation must begin and end reasonably fast
dblplaw["beta"] = (5.0, 10.0) # Formation must begin and end reasonably fast
dblplaw["massformed"] = (5.0, 12.5)
dblplaw["metallicity"] = (0.0, 2.5)

delayed = {}                   # Delayed Tau model t*e^-(t/tau)
delayed["age"] = (7.5, 12.5)   # Time since SF began: Gyr
delayed["tau"] = (0.1, 2.0)    # Timescale of decrease: Gyr
delayed["massformed"] = (5.0, 12.5)
delayed["metallicity"] = (0.0, 2.5)

lognormal = {}
lognormal["tmax"] = (1.0, 6.0)
lognormal["fwhm"] = (0.1, 3.0)
lognormal["massformed"] = (5.0, 12.5)
lognormal["metallicity"] = (0.0, 2.5)

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0.0, 2.0)

nebular = {}
nebular["logU"] = (-4.0,-2.0)

chi_squ_vals = {}

for filename in datafiles:

    # Get the galaxy object and a priori model compenents dictionary.
    galaxy, model_components = import_spectrum(filename)
    # Calculate redshift constraints.
    z_low, z_high = model_components["redshift"] - 0.001,  model_components["redshift"] + 0.001

    try:
        f = open("pipes/posterior/exponential_burst_r4/"+filename+".h5")
        f.close()
        fit_instructions = {
        "redshift"     : (z_low, z_high), # z varies tight_layout around z_obs.
        "t_bc"         : (0.013, 0.021),  # Constraints from Murray 2011.
        "veldisp"      : (50.0, 450.0),   # Constrained by Faber-Jackson. TODO: Lookup Minkowski 1962!
        "exponential1" : exponential1,    # Add the exp SFH component.
        "exponential2" : exponential2,    # Add the burst component.
        "dust"         : dust,
        "nebular"      : nebular
        }
        # Do a fit with both an old exponential component and recent burst.
        fit = pipes.fit(galaxy, fit_instructions, run="exponential_burst_r4")
        fit.fit(verbose=False)
        # Create a dictionary for storying posterior sample distribution widths.
        chi_squ_vals["exponential_burst_r4"] = chi_squared(galaxy, fit)
        # Select the fit with lowest chi-ssquared value and plot it.
        plt.tight_layout()
        fig = fit.plot_sfh_posterior(save=True, show=False)
    except IOError:
        print("Run not exponential_burst_r4 not accessible for "+filename)

    try:
        f = open("pipes/posterior/dblplaw_burst_r4/"+filename+".h5")
        f.close()
        fit_instructions = {
        "redshift"     : (z_low, z_high), # z varies tight_layout around z_obs.
        "t_bc"         : (0.013, 0.021),  # Constraints from Murray 2011.
        "veldisp"      : (50.0, 450.0),   # Constrained by Faber-Jackson. TODO: Lookup Minkowski 1962!
        "dblplaw"      : dblplaw,         # Add the dblplaw SFH component.
        "exponential2" : exponential2,    # Add the burst component.
        "dust"         : dust,
        "nebular"      : nebular
        }
        # Do a fit with both an old delayed component and recent burst.
        fit = pipes.fit(galaxy, fit_instructions, run="dblplaw_burst_r4")
        fit.fit(verbose=False)
        # Create a dictionary for storying posterior sample distribution widths.
        chi_squ_vals["dblplaw_burst_r4"] = chi_squared(galaxy, fit)
        # Select the fit with lowest chi-squared value and plot it.
        plt.tight_layout()
        fig = fit.plot_sfh_posterior(save=True, show=False)
    except IOError:
        print("Run not dblplaw_burst_r4 not accessible for "+filename)

    try:
        f = open("pipes/posterior/delayed_burst_r4/"+filename+".h5")
        f.close()
        fit_instructions = {
        "redshift"     : (z_low, z_high), # z varies tight_layout around z_obs.
        "t_bc"         : (0.013, 0.021),  # Constraints from Murray 2011.
        "veldisp"      : (50.0, 450.0),   # Constrained by Faber-Jackson. TODO: Lookup Minkowski 1962!
        "delayed"      : delayed,         # Add the delayed SFH component.
        "exponential2" : exponential2,    # Add the burst component.
        "dust"         : dust,
        "nebular"      : nebular
        }
        # Do a fit with both an old delayed component and recent burst.
        fit = pipes.fit(galaxy, fit_instructions, run="delayed_burst_r4")
        fit.fit(verbose=False)
        # Create a dictionary for storying posterior sample distribution widths.
        chi_squ_vals["delayed_burst_r4"] = chi_squared(galaxy, fit)
        # Select the fit with lowest chi-squared value and plot it.
        plt.tight_layout()
        fig = fit.plot_sfh_posterior(save=True, show=False)
    except IOError:
        print("Run not delayed_burst_r4 not accessible for "+filename)

    try:
        f = open("pipes/posterior/lognormal_burst_r4/"+filename+".h5")
        f.close()
        fit_instructions = {
        "redshift"     : (z_low, z_high), # z varies tight_layout around z_obs.
        "t_bc"         : (0.013, 0.021),  # Constraints from Murray 2011.
        "veldisp"      : (50.0, 450.0),   # Constrained by Faber-Jackson. TODO: Lookup Minkowski 1962!
        "lognormal"    : lognormal,       # Add the lognormal SFH component.
        "exponential2" : exponential2,    # Add the burst component.
        "dust"         : dust,
        "nebular"      : nebular
        }
        # Do a fit with both an old delayed component and recent burst.
        fit = pipes.fit(galaxy, fit_instructions, run="lognormal_burst_r4")
        fit.fit(verbose=False)
        # Create a dictionary for storying posterior sample distribution widths.
        chi_squ_vals["lognormal_burst_r4"] = chi_squared(galaxy, fit)
        # Select the fit with lowest chi-squared value and plot it.
        plt.tight_layout()
        fig = fit.plot_sfh_posterior(save=True, show=False)

    except IOError:
        print("Run not lognormal_burst_r4 not accessible for "+filename)


    # # Reload all the saved fits.

    # fit = pipes.fit(galaxy, fit_instructions, run="exponential_burst_r4")
    # fit.fit(verbose=False)
    # chi_squ_vals = {"exponential_burst_r4" : chi_squared(galaxy, fit)}
    #
    # fit = pipes.fit(galaxy, fit_instructions, run="dblplaw_burst_r4")
    # fit.fit(verbose=True)
    # chi_squ_vals["dblplaw_burst_r4"] = chi_squared(galaxy, fit)
    #
    # fit = pipes.fit(galaxy, fit_instructions, run="delayed_burst_r4")
    # fit.fit(verbose=False)
    # chi_squ_vals["delayed_burst_r4"] = chi_squared(galaxy, fit)
    #
    # fit = pipes.fit(galaxy, fit_instructions, run="lognormal_burst_r4")
    # fit.fit(verbose=False)
    # chi_squ_vals["lognormal_burst_r4"] = chi_squared(galaxy, fit)

    # Get the functional form with the lowest chi-squared value.
    best_func = min(chi_squ_vals, key=lambda k: chi_squ_vals[k])
    # Select the fit with lowest chi-squared value and plot it.
    fit = pipes.fit(galaxy, fit_instructions, run=best_func)

    # fit.plot_spectrum_posterior(show=True, save=False)  # Shows the input and fitted spectrum/photometry
    # fit.plot_sfh_posterior(show=True, save=False)       # Shows the fitted star-formation history
    # fit.plot_1d_posterior(show=True, save=False)        # Shows 1d posterior probability distributions
    # fit.plot_corner(show=True, save=False)              # Shows 1d and 2d posterior probability distributions

    print_posterior(fit)
    plot_corner(fit, names=["exponential2:massformed", "lognormal:fwhm", "veldisp"])
    
    # import matplotlib as mpl
    # fig = plt.figure(figsize=(12, 7))
    # gs = mpl.gridspec.GridSpec(7, 4, hspace=3., wspace=0.1)
    #
    # ax1 = plt.subplot(gs[:4, :])
    #
    # labels = ["lognormal:fwhm", "exponential2:massformed"]
    #
    # post_quantities = dict(zip(labels, [fit.posterior.samples[l] for l in labels]))
    #
    # axes = []
    # for i in range(4):
    #     axes.append(plt.subplot(gs[4:, i]))
    #     pipes.plotting.plot_corner(post_quantities[labels[i]], axes[-1], smooth=True, label=labels[i])
    #
    # plt.show()
