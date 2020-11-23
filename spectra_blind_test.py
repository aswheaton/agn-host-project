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

    exp = {"age" : 2.85,          # Gyr
           "tau" : 0.5,         # Gyr
           "massformed" : 9.05,  # log_10(M*/M_solar)
           "metallicity" : 0.87  # Z/Z_oldsolar
           }

    dblplaw = {}
    dblplaw["tau"] = 2.25
    dblplaw["alpha"] = 8.5
    dblplaw["beta"] = 11.
    dblplaw["massformed"] = 9.5
    dblplaw["metallicity"] = 0.14

    burst1 = {}
    burst1["age"] = 5.0
    burst1["massformed"] = 8.05
    burst1["metallicity"] = 0.1

    burst2 = {}
    burst2["age"] = 3.8
    burst2["massformed"] = 7.73
    burst2["metallicity"] = 0.15

    burst3 = {}
    burst3["age"] = 2.15
    burst3["massformed"] = 7.25
    burst3["metallicity"] = 0.52



    delayed = {}                   # Delayed Tau model t*e^-(t/tau)
    delayed["age"] = 10.0           # Time since SF began: Gyr
    delayed["tau"] = 1.45           # Timescale of decrease: Gyr
    delayed["massformed"] = 10.8
    delayed["metallicity"] = 1.25

    lognormal1 = {}
    lognormal1["tmax"] = 10.9
    lognormal1["fwhm"] = 0.15
    lognormal1["massformed"] = 8.15
    lognormal1["metallicity"] = 0.82

    lognormal2 = {}
    lognormal2["tmax"] = 5.65
    lognormal2["fwhm"] = 0.05
    lognormal2["massformed"] = 8.14
    lognormal2["metallicity"] = 0.52

    lognormal3 = {}                       # lognormal SFH
    lognormal3["tmax"] = 8.15              # Age of Universe at peak SF: Gyr
    lognormal3["fwhm"] = 1.5
    lognormal3["massformed"] = 8.85
    lognormal3["metallicity"] = 1.25

    dust = {"type" : "Calzetti", # Define the shape of the attenuation curve
            "Av"   : 0.25,        # Extinction, in magnitudes.
            "eta"  : 1.5         # Extra dust for young stars: multiplies Av
            }

    nebular = {"logU" : -3.5}    # log_10(ionization parameter)

    model_components = {"redshift"    : 0.75,     # Observed redshift.
                        # "exponential" : exp,
                        # "dblplaw"     : dblplaw,
                        "burst1"      : burst1,
                        "burst2"      : burst2,
                        "burst3"      : burst3,
                        "lognormal1"   : lognormal1,
                        # "lognormal2"   : lognormal2,
                        # "lognormal3"   : lognormal3,
                        # "delayed"     : delayed,
                        "dust"        : dust,
                        "nebular"     : nebular,
                        "t_bc"        : 0.01,    # Lifetime of birth clouds (Gyr)
                        "veldisp"     : 323.0    # km/s
                        }

    # Create the model galaxy object with the defined parameters.
    model_ID = "model010"
    obs_wavs = np.arange(2500.0, 9500.0, 5.0)
    model = pipes.model_galaxy(model_components, spec_wavs=obs_wavs)
    fig, ax = model.plot(show=False)
    sfh = model.sfh.plot(show=False)
    ax[0].vlines(rframe_lines, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='r', linestyles='dashed', label=rframe_line_labels)
    # plt.tight_layout()
    plt.savefig("pipes/plots/" + model_ID + "_sfh.pdf")
    plt.show()

    export_spectrum(model_ID, model, 0.12)

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
