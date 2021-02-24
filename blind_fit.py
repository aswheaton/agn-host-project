import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
from pipes_utils import *

from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

# datafiles = ["phil_model_01"]
# datafiles = ["phil_model_02", "phil_model_03", "phil_model_04"]
# datafiles = ["phil_model_05", "phil_model_06", "phil_model_07"]
datafiles = ["phil_model_08", "phil_model_09", "phil_model_10"]

exponential = {}
exponential["age"] = (3.5, 10.0)   # Gyr
exponential["tau"] = (0.1, 2.0)    # Gyr
exponential["massformed"] = (0.0, 15.0)   # log_10(M*/M_solar)
exponential["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

# dblplaw = {}
# dblplaw["tau"] = 11.72
# dblplaw["alpha"] = 290.39
# dblplaw["beta"] = 240.33
# dblplaw["massformed"] = 9.72
# dblplaw["metallicity"] = 0.85

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

fit_instructions = {
"redshift"    : (0.0, 0.1),   # Obs. redshift from 0-10.
"t_bc"        : (0.005, 0.015),
"veldisp"     : (150.0, 200.0),
"exponential" : exponential, # Add the exp SFH component.
"dust"        : dust
}

for filename in datafiles:

    print("Running initial exponential fit for {}...".format(filename)),

    # Do an initial fit with only an exponential compontent, over a large parameter space.
    galaxy, model_components = import_spectrum(filename)
    fit = pipes.fit(galaxy, fit_instructions, run="exponential_noburst")
    fit.fit(verbose=False)

    # Now run a new fit with other functional components with age constrained to one standard deviation around the exponential age.
    age_lower_bound = np.percentile(fit_instructions["exponential"]["age"], 16)
    age_upper_bound = np.percentile(fit_instructions["exponential"]["age"], 84)

    # Create a dictionary for storying posterior sample distribution widths.
    chi_squ_vals = {"exponential" : chi_squared(galaxy, fit)}

    # dblplaw["age"] = ()
    # fit_instructions["dblplaw"] = dblplaw
    # fit = pipes.fit(galaxy, fit_instructions, run="dblplaw_noburst")
    # fit.fit(verbose=True)
    # fit_instructions.pop("dblplaw", None)

    fit_instructions.pop("exponential", None)
    delayed["age"] = (age_lower_bound, age_upper_bound)
    fit_instructions["delayed"] = delayed
    fit = pipes.fit(galaxy, fit_instructions, run="delayed_noburst")
    fit.fit(verbose=False)
    chi_squ_vals["delayed"] = chi_squared(galaxy, fit)

    fit_instructions.pop("delayed", None)
    lognormal["tmax"] = (13.5 - age_upper_bound, 13.5 - age_lower_bound)
    fit_instructions["lognormal"] = lognormal
    fit = pipes.fit(galaxy, fit_instructions, run="lognormal_noburst")
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
    fit = pipes.fit(galaxy, fit_instructions, run=best_func+"_burst")
    fit.fit(verbose=False)

    chi_squ_vals[best_func+"_burst"] = chi_squared(galaxy, fit)

    # Check if the burst component improves the fit and save the appropriate plot.
    if chi_squ_vals[best_func] >= chi_squ_vals[best_func+"_burst"]:
        fig = fit.plot_sfh_posterior(save=True, show=False)
    else:
        fit = pipes.fit(galaxy, fit_instructions, run=best_func+"_noburst")
        fig = fit.plot_sfh_posterior(save=True, show=False)

print ('parameter     median     16th percentile     84th percentile')
for key in fit.posterior.samples.keys():
    print(key+": ", np.median(fit.posterior.samples[key]), np.percentile(fit.posterior.samples[key], 16), np.percentile(fit.posterior.samples[key], 84))
