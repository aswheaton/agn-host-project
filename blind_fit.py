import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
from pipes_utils import *

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# datafiles = ["phil_model_01", "phil_model_02", "phil_model_03", "phil_model_04", "phil_model_05",
#              "phil_model_06", "phil_model_07", "phil_model_08", "phil_model_09", "phil_model_10"]

datafiles = ["phil_model_01"]

exponential = {}
exponential["age"] = (3.5, 10.0)   # Gyr
exponential["tau"] = (0.0, 2.0)    # Gyr
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
delayed["tau"] = (0.0, 2.0)           # Timescale of decrease: Gyr
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

    # Do an initial fit with only an exponential compontent, over a large parameter space.
    print("Running initial exponential fit for {}...".format(filename)),
    galaxy, model_components = import_spectrum(filename)
    fit = pipes.fit(galaxy, fit_instructions, run="exp_noburst")
    fit.fit(verbose=False)
    print("Done!")

    # Now run a new fit with other functional components with age constrained to one standard deviation around the exponential age.
    age_lower_bound = np.percentile(fit_instructions["exponential"]["age"], 16)
    age_upper_bound = np.percentile(fit_instructions["exponential"]["age"], 84)

    # Create a dictionary for storying posterior sample distribution widths.
    sigma_values = {"exponential" : np.percentile(fit.posterior.samples["exponential:age"], 84) - np.percentile(fit.posterior.samples["exponential:age"], 16)}

    fit_instructions.pop("exponential", None)
    # dblplaw["age"] = ()
    # fit_instructions["dblplaw"] = dblplaw
    # fit = pipes.fit(galaxy, fit_instructions, run="dblplaw_noburst")
    # fit.fit(verbose=True)
    # sigma_values["dblplaw"] = np.percentile(fit.posterior.samples["age"], 84) - np.percentile(fit.posterior.samples["age"], 16)
    #
    print("Running delayed fit for {}...".format(filename)),
    # fit_instructions.pop("dblplaw", None)
    delayed["age"] = (age_lower_bound, age_upper_bound)
    fit_instructions["delayed"] = delayed
    fit = pipes.fit(galaxy, fit_instructions, run="delayed_noburst")
    fit.fit(verbose=False)
    sigma_values["delayed"] = np.percentile(fit.posterior.samples["delayed:age"], 84) - np.percentile(fit.posterior.samples["delayed:age"], 16)
    print("Done!")

    print("Running lognormal fit for {}...".format(filename)),
    fit_instructions.pop("delayed", None)
    lognormal["tmax"] = (13.5 - age_upper_bound, 13.5 - age_lower_bound)
    fit_instructions["lognormal"] = lognormal
    fit = pipes.fit(galaxy, fit_instructions, run="lognormal_noburst")
    fit.fit(verbose=False)
    sigma_values["lognormal"] = np.percentile(fit.posterior.samples["lognormal:tmax"], 84) - np.percentile(fit.posterior.samples["lognormal:tmax"], 16)
    print("Done!")

    # Get the functional form with the "tightest" contraints on age.
    best_func = min(sigma_values, key=lambda k: sigma_values[k])

    print("Running a {} and burst fit for {}...".format(best_func, filename)),
    # Now run a fit with that functional form AND a burst in recent history.
    fit_instructions.pop("lognormal", None)
    fit_instructions[best_func] = eval(best_func)

    burst = {}
    burst["age"] = (0.0, age_lower_bound)
    burst["massformed"] = (0.0, 15.0)
    burst["metallicity"] = (0.0, 2.5)

    fit_instructions["burst"] = burst
    fit = pipes.fit(galaxy, fit_instructions, run=best_func+"_burst")
    fit.fit(verbose=False)
    print("Done!")


print ('parameter     median     16th percentile     84th percentile')
for key in fit_instructions.keys():
    try:
        for subkey in fit_instructions[key].keys():
            i = key+":"+subkey
            print(i + ': ', np.median(fit.posterior.samples[i]), np.percentile(fit.posterior.samples[i], 16), np.percentile(fit.posterior.samples[i], 84))
    except AttributeError:
        i = key
        print(i + ': ', np.median(fit.posterior.samples[i]), np.percentile(fit.posterior.samples[i], 16), np.percentile(fit.posterior.samples[i], 84))

fig = fit.plot_spectrum_posterior(save=False, show=True)
fig = fit.plot_sfh_posterior(save=False, show=True)
