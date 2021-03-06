import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
import sys
from pipes_utils import *

from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

rframe_lines = [9532, 9068, 8662, 8542, 8489,7150, 6731, 6716, 6583, 6563, 6180,
                5892, 5175, 5007, 4959, 4861, 4340, 4101, 3727, 2800, 1549, 1400
                ]

rframe_line_labels = ["SIII", "SIII", "CaII", "CaII", "CaII", "TiO", "SII",
                      "SII", "NII", "H-alpha", "TiO", "NaD", "MgIb", "OIII",
                      "OIII", "H-beta", "H-gamma", "H-delta", "OII", "MgII",
                      "CIV", "SiIV"
                      ]

datafiles = [sys.argv[1]]
# datafiles = ["phil_model_02", "phil_model_03", "phil_model_04"]
# datafiles = ["phil_model_05", "phil_model_06", "phil_model_07"]
# datafiles = ["phil_model_08", "phil_model_09", "phil_model_10"]

exponential = {}
exponential["age"] = (3.5, 10.0)   # Gyr
exponential["tau"] = (0.1, 2.0)    # Gyr
exponential["massformed"] = (0.0, 15.0)   # log_10(M*/M_solar)
exponential["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

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
    "redshift"    : (0.0, 0.1),   # Obs. redshift from 0-10.
    "t_bc"        : (0.005, 0.015),
    "veldisp"     : (150.0, 200.0),
    "exponential" : exponential, # Add the exp SFH component.
    "dust"        : dust
    }

    # Do an initial fit with only an exponential compontent, over a large parameter space.
    galaxy, model_components = import_spectrum(filename)
    fit = pipes.fit(galaxy, fit_instructions, run="r1_exponential_noburst")
    fit.fit(verbose=False)

    # Now run a new fit with other functional components with age constrained to one standard deviation around the exponential age.
    age_lower_bound = np.percentile(fit_instructions["exponential"]["age"], 16)
    age_upper_bound = np.percentile(fit_instructions["exponential"]["age"], 84)

    # Create a dictionary for storying posterior sample distribution widths.
    chi_squ_vals = {"exponential" : chi_squared(galaxy, fit)}

    fit_instructions.pop("exponential", None)
    delayed["age"] = (age_lower_bound, age_upper_bound)
    fit_instructions["delayed"] = delayed
    fit = pipes.fit(galaxy, fit_instructions, run="r1_delayed_noburst")
    fit.fit(verbose=False)
    chi_squ_vals["delayed"] = chi_squared(galaxy, fit)

    fit_instructions.pop("delayed", None)
    lognormal["tmax"] = (13.5 - age_upper_bound, 13.5 - age_lower_bound)
    fit_instructions["lognormal"] = lognormal
    fit = pipes.fit(galaxy, fit_instructions, run="r1_lognormal_noburst")
    fit.fit(verbose=False)
    chi_squ_vals["lognormal"] = chi_squared(galaxy, fit)

    # Get the functional form with the lowest chi-squared value.
    best_func = min(chi_squ_vals, key=lambda k: chi_squ_vals[k])

    # Now run a fit with that functional form AND a burst in recent history.
    fit_instructions.pop("lognormal", None)
    fit_instructions[best_func] = eval(best_func)

    burst = {}
    burst["age"] = (0.1, age_lower_bound)
    burst["massformed"] = (0.0, 15.0)
    burst["metallicity"] = (0.0, 2.5)

    fit_instructions["burst"] = burst
    fit = pipes.fit(galaxy, fit_instructions, run="r1_"+best_func+"_burst")
    fit.fit(verbose=False)

    chi_squ_vals[best_func+"_burst"] = chi_squared(galaxy, fit)

    # Check if the burst component improves the fit and save the appropriate plot.
    if chi_squ_vals[best_func] >= chi_squ_vals[best_func+"_burst"]:
        fig = fit.plot_sfh_posterior(save=True, show=False)
    else:
        fit = pipes.fit(galaxy, fit_instructions, run="r1_"+best_func+"_noburst")
        fig = fit.plot_sfh_posterior(save=True, show=False)

    print_posterior(fit)
