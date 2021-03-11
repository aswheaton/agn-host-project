import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
import sys
from pipes_utils import *

# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

datafiles = [sys.argv[1]]

exponential1 = {}
exponential1["age"] = (7.5, 12.5)   # Gyr
exponential1["tau"] = (0.5, 2.0)    # Gyr
exponential1["massformed"] = (5.0, 15.0)   # log_10(M*/M_solar)
exponential1["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

exponential2 = {}
exponential2["age"] = (0.0, 3.5)   # Gyr, lifetime of F type stars
exponential2["tau"] = (0.1, 2.0)    # Gyr
exponential2["massformed"] = (0.0, 15.0)   # log_10(M*/M_solar)
exponential2["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

dblplaw = {}
dblplaw["tau"] = (0.0, 4.5) # Do not let the peak occur beyond ~6.5 Gyr
dblplaw["alpha"] = (5.0, 10.0) # Formation must begin and end reasonably fast
dblplaw["beta"] = (5.0, 10.0) # Formation must begin and end reasonably fast
dblplaw["massformed"] = (5.0, 15.0)
dblplaw["metallicity"] = (0.0, 2.5)

delayed = {}                   # Delayed Tau model t*e^-(t/tau)
delayed["age"] = (7.5, 12.5)           # Time since SF began: Gyr
delayed["tau"] = (0.1, 2.0)           # Timescale of decrease: Gyr
delayed["massformed"] = (5.0, 15.0)
delayed["metallicity"] = (0.0, 2.5)

lognormal = {}
lognormal["tmax"] = (1.0, 6.0 )
lognormal["fwhm"] = (0.1, 3.0)
lognormal["massformed"] = (5.0, 15.0)
lognormal["metallicity"] = (0.0, 2.5)

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0.0, 2.0)

nebular = {}
nebular["logU"] = (-4.0,-2.0)

for filename in datafiles:

    print("Running initial exponential fit for {}...".format(filename)),
    # Get the galaxy object and a priori model compenents dictionary.
    galaxy, model_components = import_spectrum(filename)

    # Calculate redshift constraints.
    z_low, z_high = model_components["redshift"] - 0.001,  model_components["redshift"] + 0.001

    # Create (or reset) the fit instructions dictionary.
    fit_instructions = {
    "redshift"     : (z_low, z_high),   # Obs. redshift from 0-10.
    "t_bc"         : (0.005, 2.0),
    "veldisp"      : (1.0, 1000.0),
    "veldisp_prior": "log_10",
    "exponential1" : exponential1, # Add the exp SFH component.
    "exponential2" : exponential2, # Add the burst component.
    "dust"         : dust,
    "nebular"      : nebular
    }

    # Do an initial fit with only an exponential compontent, over a large parameter space.
    fit = pipes.fit(galaxy, fit_instructions, run="exponential_burst_final")
    fit.fit(verbose=False)

    # Create a dictionary for storying posterior sample distribution widths.
    chi_squ_vals = {"exponential_burst_final" : chi_squared(galaxy, fit)}

    fit_instructions.pop("exponential1")
    fit_instructions["dblplaw"] = dblplaw
    fit = pipes.fit(galaxy, fit_instructions, run="dblplaw_burst_final")
    fit.fit(verbose=True)
    chi_squ_vals["dblplaw_burst_final"] = chi_squared(galaxy, fit)

    fit_instructions.pop("dblplaw", None)
    fit_instructions["delayed"] = delayed
    fit = pipes.fit(galaxy, fit_instructions, run="delayed_burst_final")
    fit.fit(verbose=False)
    chi_squ_vals["delayed_burst_final"] = chi_squared(galaxy, fit)

    fit_instructions.pop("delayed", None)
    fit_instructions["lognormal"] = lognormal
    fit = pipes.fit(galaxy, fit_instructions, run="lognormal_burst_final")
    fit.fit(verbose=False)
    chi_squ_vals["lognormal_burst_final"] = chi_squared(galaxy, fit)

    # Get the functional form with the lowest chi-squared value.
    # best_func = min(chi_squ_vals, key=lambda k: chi_squ_vals[k])
    # Select the fit with lowest chi-squared value and plot it.
    # fit = pipes.fit(galaxy, fit_instructions, run=best_func)
    # plt.tight_layout()
    # fig = fit.plot_sfh_posterior(save=True, show=False)

    print ('parameter     median     16th percentile     84th percentile')
    for key in fit.posterior.samples.keys():
        print(key+": ", np.median(fit.posterior.samples[key]), np.percentile(fit.posterior.samples[key], 16), np.percentile(fit.posterior.samples[key], 84))
