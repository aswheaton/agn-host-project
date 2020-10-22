import bagpipes as pipes
import numpy as np

exp = {}                          # Tau model star formation history component
exp["age"] = 3.                   # Gyr
exp["tau"] = 0.75                 # Gyr
exp["massformed"] = 9.            # log_10(M*/M_solar)
exp["metallicity"] = 0.5          # Z/Z_oldsolar

dust = {}                         # Dust component
dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
dust["Av"] = 0.2                  # magnitudes

model_components = {}                   # The model components dictionary
model_components["redshift"] = 1.0      # Observed redshift
model_components["exponential"] = exp
model_components["dust"] = dust

def bin_spec(spectrum, n_bin):
    """
    Args:
        spectrum (2Darray): two or three column array of wavelengths, fluxes,
        and (optionally) flux errors.
        n_bin (int): number of bins into which the fluxes are sorted.
    Returns:
        binspec (2Darray): bins up two or three column spectral data by a specified factor.
    """

    n_bin = int(n_bin)
    n_bins = len(spectrum)/n_bin
    binspec = np.zeros((n_bins, spectrum.shape[1]))

    for i in range(binspec.shape[0]):
        spec_slice = spectrum[i*n_bin:(i+1)*n_bin, :]
        binspec[i, 0] = np.mean(spec_slice[:, 0])
        binspec[i, 1] = np.mean(spec_slice[:, 1])

        if spectrum.shape[1] == 3:
            binspec[i,2] = (1./float(n_bin)
                            *np.sqrt(np.sum(spec_slice[:, 2]**2)))
    return(binspec)

def load_xshooter_spec(ID):
    data = np.loadtxt("data/20200127_xshoot_corr.asci", dtype="float")
    lambdas, fluxes = data[:,0], data[:,1]
