import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
import sys
from pipes_utils import *

from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

datafiles = [sys.argv[1]]
# datafiles = ["phil_model_02", "phil_model_03", "phil_model_04"]
# datafiles = ["phil_model_05", "phil_model_06", "phil_model_07"]
# datafiles = ["phil_model_08", "phil_model_09", "phil_model_10"]

exponential1 = {}
exponential1["age"] = (3.5, 10.0)   # Gyr
exponential1["tau"] = (0.1, 2.0)    # Gyr
exponential1["massformed"] = (0.0, 15.0)   # log_10(M*/M_solar)
exponential1["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

exponential2 = {}
exponential2["age"] = (0.0, 10.0)   # Gyr
exponential2["tau"] = (0.1, 2.0)    # Gyr
exponential2["massformed"] = (0.0, 15.0)   # log_10(M*/M_solar)
exponential2["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

dblplaw = {}
dblplaw["tau"] = (0.0, 15.0)
dblplaw["alpha"] = (0.0, 10.0)
dblplaw["beta"] = (0.0, 10.0)
dblplaw["massformed"] = (0.0, 15)
dblplaw["metallicity"] = (0.0, 2.5)

delayed = {}                   # Delayed Tau model t*e^-(t/tau)
delayed["age"] = (3.5, 10.0)           # Time since SF began: Gyr
delayed["tau"] = (0.1, 2.0)           # Timescale of decrease: Gyr
delayed["massformed"] = (0.0, 15.0)
delayed["metallicity"] = (0.0, 2.5)

lognormal = {}
lognormal["tmax"] = (3.5, 10.5)
lognormal["fwhm"] = (0.1, 5.0)
lognormal["massformed"] = (0.0, 15.0)
lognormal["metallicity"] = (0.0, 2.5)

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0.0, 1.0)

for filename in datafiles:

    print("Running initial exponential fit for {}...".format(filename)),

    # Create (or reset) the fit instructions dictionary.
    fit_instructions = {
    "redshift"     : (0.0, 0.1),   # Obs. redshift from 0-10.
    "t_bc"         : (0.005, 0.015),
    "veldisp"      : (150.0, 200.0),
    "exponential1" : exponential1, # Add the exp SFH component.
    "exponential2" : exponential2, # Add the burst component.
    "dust"         : dust
    }

    # Do an initial fit with only an exponential compontent.
    galaxy, model_components = import_spectrum(filename)
    fit = pipes.fit(galaxy, fit_instructions, run="r2_exponential_burst")
    fit.fit(verbose=False)

    # Create a dictionary for storying posterior sample distribution widths.
    chi_squ_vals = {"r2_exponential_burst" : chi_squared(galaxy, fit)}

    fit_instructions.pop("exponential1")
    fit_instructions["dblplaw"] = dblplaw
    fit = pipes.fit(galaxy, fit_instructions, run="r2_dblplaw_burst")
    fit.fit(verbose=True)
    chi_squ_vals["r2_dblplaw_burst"] = chi_squared(galaxy, fit)

    fit_instructions.pop("dblplaw", None)
    fit_instructions["delayed"] = delayed
    fit = pipes.fit(galaxy, fit_instructions, run="r2_delayed_burst")
    fit.fit(verbose=False)
    chi_squ_vals["r2_delayed_burst"] = chi_squared(galaxy, fit)

    fit_instructions.pop("delayed", None)
    fit_instructions["lognormal"] = lognormal
    fit = pipes.fit(galaxy, fit_instructions, run="r2_lognormal_burst")
    fit.fit(verbose=False)
    chi_squ_vals["r2_lognormal_burst"] = chi_squared(galaxy, fit)

    # Get the functional form with the lowest chi-squared value.
    best_func = min(chi_squ_vals, key=lambda k: chi_squ_vals[k])
    # Select the fit with lowest chi-squared value and plot it.
    fit = pipes.fit(galaxy, fit_instructions, run=best_func)
    fig = fit.plot_sfh_posterior()

    print_posterior(fit)
