import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
from pipes_utils import *

def make_model():
    exp = {"age" : 3.0,          # Gyr
           "tau" : 0.75,         # Gyr
           "massformed" : 9.0,   # log_10(M*/M_solar)
           "metallicity" : 0.5   # Z/Z_oldsolar
           }

    dust = {"type" : "Calzetti", # Define the shape of the attenuation curve
            "Av"   : 0.2         # Extinction, in magnitudes.
            }

    model_components = {"redshift"    : 1.0, # Observed redshift.
                        "exponential" : exp,
                        "dust"        : dust
                        }

    obs_wavs = np.arange(2500., 7500., 5.)
    model = pipes.model_galaxy(model_components, spec_wavs=obs_wavs)
    model.plot()
    export_spectrum("model1", model)

def main():

    burst = {"age"         : (0.0, 3.5), # Vary age from 10 to 15 Gyr
             "metallicity" : (0.0, 2.5),  # Vary metallicity from 0 to 2.5 Solar
             "massformed"  : (0.0, 13.,0) # Vary log_10(mass formed) from 0 to 13
            }

    sfh_exponential = {"age"         : (3.5, 10.0),   # Gyr
                       "tau"         : (0.0, 2.0),    # Gyr
                       "massformed"  : (0.0, 15.0),   # log_10(M*/M_solar)
                       "metallicity" : (0.0, 2.5)     # Z/Z_oldsolar
                       }

    fit_instructions = {"burst"       : burst,
                        "exponential" : sfh_exponential, # Add the exp SFH component.
                        "redshift"    : (0.0, 0.1)      # Obs. redshift from 0-10.
                        }

    # galaxy, model_components = import_spectrum("model1")
    # galaxy.plot()
    # fit = pipes.fit(galaxy, fit_instructions)
    # fit.fit(verbose=True)

    # galaxy = pipes.galaxy("20200127_xshoot_corr", load_xshooter, photometry_exists=False)
    # galaxy.plot()
    # fit = pipes.fit(galaxy, fit_instructions)
    # fit.fit(verbose=True)

    galaxy = pipes.galaxy("host_hyz_specwerr", load_xshooter_2, photometry_exists=False)
    galaxy.plot()
    fit = pipes.fit(galaxy, fit_instructions)
    fit.fit(verbose=False)

    # galaxy = pipes.galaxy("15oi_ubv+vis+nir_med5_errmask", load_xshooter_2, photometry_exists=False)
    # galaxy.plot()
    # fit = pipes.fit(galaxy, fit_instructions)
    # fit.fit(verbose=True)

    fit.plot_spectrum_posterior()  # Shows the input and fitted spectrum/photometry
    fit.plot_sfh_posterior()       # Shows the fitted star-formation history
    fit.plot_1d_posterior()        # Shows 1d posterior probability distributions
    fit.plot_corner()              # Shows 1d and 2d posterior probability distributions

    # data = np.loadtxt("data/20200127_xhoot_med15.asci", dtype="float")
    # lambdas, fluxes = data[:,0], data[:,1]
    # plt.plot(lambdas, fluxes)
    # plt.show()
    #
    # data = np.loadtxt("data/20200127_xshoot_corr.asci", dtype="float")
    # lambdas, fluxes = data[:,0], data[:,1]
    # plt.plot(lambdas, fluxes)
    # plt.show()


main()
