import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
from pipes_utils import *

from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

datafiles = ["phil_model_01", "phil_model_02", "phil_model_03", "phil_model_04",
             "phil_model_05", "phil_model_06", "phil_model_07", "phil_model_08",
             "phil_model_09", "phil_model_10"]

for filename in datafiles:

    galaxy, model_components = import_spectrum(filename)
    # Create the model galaxy object with the defined parameters.
    model_ID = filename
    obs_wavs = np.arange(1000.0, 10000.0, 5.0)
    model = pipes.model_galaxy(model_components, spec_wavs=obs_wavs)

    sfh = model.sfh.plot(show=False)
    plt.tight_layout()
    plt.savefig("pipes/plots/prior/" + model_ID + "_sfh.pdf")
