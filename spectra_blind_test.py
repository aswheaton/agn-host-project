import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
from pipes_utils import *

# Make plots nice.
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

rframe_lines = [9532, 9068, 8662, 8542, 8489,7150, 6731, 6716, 6583, 6563, 6180,
                5892, 5175, 5007, 4959, 4861, 4340, 4101, 3727, 2800, 1549, 1400
                ]

rframe_line_labels = ["SIII", "SIII", "CaII", "CaII", "CaII", "TiO", "SII",
                      "SII", "NII", "H-alpha", "TiO", "NaD", "MgIb", "OIII",
                      "OIII", "H-beta", "H-gamma", "H-delta", "OII", "MgII",
                      "CIV", "SiIV"
                      ]

def main():

    exp = {"age" : 3.0,          # Gyr
           "tau" : 0.75,         # Gyr
           "massformed" : 10.0,  # log_10(M*/M_solar)
           "metallicity" : 0.5   # Z/Z_oldsolar
           }

    dblplaw = {}
    dblplaw["tau"] = 5.
    dblplaw["alpha"] = 2.5
    dblplaw["beta"] = 10.
    dblplaw["massformed"] = 10.
    dblplaw["metallicity"] = 0.5

    burst1 = {}
    burst1["age"] = 0.2
    burst1["massformed"] = 7.
    burst1["metallicity"] = 1.

    burst2 = {}
    burst2["age"] = 3.
    burst2["massformed"] = 7.5
    burst2["metallicity"] = 0.2

    dust = {"type" : "Calzetti", # Define the shape of the attenuation curve
            "Av"   : 0.2,        # Extinction, in magnitudes.
            "eta"  : 3.0         # Extra dust for young stars: multiplies Av
            }

    nebular = {"logU" : -3.0}    # log_10(ionization parameter)

    model_components = {"redshift"    : 1.0,     # Observed redshift.
                        "exponential" : exp,
                        "dblplaw"     : dblplaw,
                        "burst1"      : burst1,
                        "burst2"      : burst2,
                        "dust"        : dust,
                        "nebular"     : nebular,
                        "t_bc"        : 0.01,    # Lifetime of birth clouds (Gyr)
                        "veldisp"     : 200.0    # km/s
                        }

    # Create the model galaxy object with the defined parameters.
    model_ID = "model001"
    obs_wavs = np.arange(2500.0, 9500.0, 5.0)
    model = pipes.model_galaxy(model_components, spec_wavs=obs_wavs)
    fig, ax = model.plot(show=False)
    sfh = model.sfh.plot(show=False)
    ax[0].vlines(rframe_lines, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='r', linestyles='dashed', label=rframe_line_labels)
    # plt.tight_layout()
    plt.savefig("pipes/plots/model001_sfh.pdf")
    plt.show()

    export_spectrum(model_ID, model, 0.01)

    # burst = {"age"         : (0.0, 3.5), # Vary age from 10 to 15 Gyr
    #          "metallicity" : (0.0, 2.5),  # Vary metallicity from 0 to 2.5 Solar
    #          "massformed"  : (0.0, 13.,0) # Vary log_10(mass formed) from 0 to 13
    #          }

    # exponential = {"age"         : (3.5, 10.0),   # Gyr
    #                "tau"         : (0.0, 2.0),    # Gyr
    #                "massformed"  : (0.0, 15.0),   # log_10(M*/M_solar)
    #                "metallicity" : (0.0, 2.5)     # Z/Z_oldsolar
    #                }

    # fit_instructions = {"burst"       : burst,
    #                     "exponential" : exponential, # Add the exp SFH component.
    #                     "redshift"    : (0.0, 0.1)   # Obs. redshift from 0-10.
    #                     }

    # galaxy, model_components = import_spectrum("model1")
    # galaxy.plot()
    # fit = pipes.fit(galaxy, fit_instructions)
    # fit.fit(verbose=True)

    # galaxy = pipes.galaxy("host_hyz_specwerr", load_xshooter, photometry_exists=False)
    # fig, ax = galaxy.plot(show=False)
    # ax[0].vlines(rframe_lines, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='r', linestyles='dashed', label=rframe_line_labels)
    # plt.show()
    # fit = pipes.fit(galaxy, fit_instructions)
    # fit.fit(verbose=False)

    # fit.plot_spectrum_posterior()  # Shows the input and fitted spectrum/photometry
    # fit.plot_sfh_posterior()       # Shows the fitted star-formation history
    # fit.plot_1d_posterior()        # Shows 1d posterior probability distributions
    # fit.plot_corner()              # Shows 1d and 2d posterior probability distributions

main()
