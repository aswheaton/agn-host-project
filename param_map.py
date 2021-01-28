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

lognormal = {}
lognormal["tmax"] = 3.5
lognormal["fwhm"] = 0.66
lognormal["massformed"] = 10.0
lognormal["metallicity"] = 0.01

burst = {}
burst["age"] = 10.5
burst["massformed"] = 9.0
burst["metallicity"] = 0.01

dust = {"type" : "Calzetti", # Define the shape of the attenuation curve
        "Av"   : 1.00,        # Extinction, in magnitudes.
        "eta"  : 1.00         # Extra dust for young stars: multiplies Av
        }

nebular = {"logU" : -2.87}    # log_10(ionization parameter)

model_components = {"redshift"    : 0.01,     # Observed redshift.
                    "lognormal"   : lognormal,
                    "burst"       : burst,
                    "dust"        : dust,
                    "nebular"     : nebular,
                    "t_bc"        : 0.01,    # Lifetime of birth clouds (Gyr)
                    "veldisp"     : 250.0    # km/s
                    }

lambda_min, lambda_max, lambda_step = 3500.0, 10000.0, 5.0
burst_age_min, burst_age_max, burst_age_step = 0.0, 9.5, -0.1
lambda_steps = int((lambda_max - lambda_min) / abs(lambda_step))
burst_steps = int((burst_age_max - burst_age_min) / abs(burst_age_step))

spectra = np.zeros((lambda_steps, burst_steps)) # Extra row for burst age.
# spectra[0,:] = np.arange(burst_age_max, burst_age_min, burst_age_step)

obs_wavs = np.arange(lambda_min, lambda_max, lambda_step)
model = pipes.model_galaxy(model_components, spec_wavs=obs_wavs)


for i in range(burst_steps):
    burst["age"] = np.arange(burst_age_max, burst_age_min, burst_age_step)[i]
    model.update(model_components)
    spectra[:,i] = model.spectrum[:,1]
    # fig, ax = model.plot(show=False)
    # sfh = model.sfh.plot(show=False)
    # ax[0].vlines(rframe_lines, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='r', linestyles='dashed', label=rframe_line_labels)
    # plt.tight_layout()
    # plt.savefig("pipes/plots/" + model_ID + "_sfh.pdf")
    plt.show()

np.savetxt("data/spectra.csv", spectra, delimiter=" ")
plt.imshow(spectra)
plt.show()

# flux_variance = np.zeros(len(obs_wavs))
# for i range(len(obs_wavs)):
flux_variance = np.var(spectra, axis=1)

plt.plot(obs_wavs, flux_variance)
reddened_lines = np.array(rframe_lines)[:-3] * 1.01
plt.vlines(reddened_lines, plt.ylim()[0], plt.ylim()[1], colors='r', linestyles='dashed', label=np.array(rframe_line_labels)[:-3])
plt.xlabel('Wavelength, Angstroms')
plt.ylabel('Flux Variance, Burst at 0.01 < z < 2.2')
plt.legend()
plt.savefig("figure.png")
plt.show()
