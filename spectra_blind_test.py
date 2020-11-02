import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt

def export_spectrum(filename, model):
    """
    Args: receives a Bagpipes model_galaxy object and writes its spectrum to a
        text file with a header containing the "secret" parameters used to
        generate that spectrum.
    Returns: None
    """
    # Generate the "secret" model components file header.
    header = ""
    for key, value in model.model_comp.items():
        try: # The value is another dictionary.
            for subkey, subvalue in value.items():
                header += subkey + ":" + str(subvalue) + "\n"
        except AttributeError: # The value is not another dictionary.
            header += key + ":" + str(value) + "\n"
    # Strip the last newline from the header.
    header = header[:-1]

    # Extract the wavelength and flux value of full spectrum from the object.
    spectrum = model.spectrum
    # Write the data with "secret" model components header to a file.
    if filename[-4:] != ".csv": filename += ".csv"
    np.savetxt("data/"+filename, spectrum, delimiter=", ", header=header)

def export_sfh(filename, model):
    pass

def load_data(filename):
    """
    Data import function for the __init__() method of bagpipes.galaxy.
    """
    # Load the data and extract into wavelength and flux arrays.
    spectrum = np.loadtxt("data/" + filename, delimiter=", ")
    spectrum_with_errs = np.zeros((spectrum.shape[0], 3))
    spectrum_with_errs[:,0:2] = np.copy(spectrum)
    spectrum_with_errs[:,2] = np.copy(np.sqrt(spectrum[:,1]))
    return(spectrum_with_errs)

def load_xshooter(ID):
    """
    Data import function for the __init__ method of bagpipes.galaxy. Loads in
    wavelengths and fluxes from the  XSHOOTER .asci files.
    """
    spectrum = np.loadtxt("data/20200127_xshoot_corr.asci", dtype="float")
    spectrum_with_errs = np.zeros((spectrum.shape[0], 3))
    spectrum_with_errs[:,0:2] = np.copy(spectrum)
    spectrum_with_errs[:,2] = np.copy(np.sqrt(spectrum[:,1]))
    return(spectrum_with_errs)

def import_spectrum(filename):
    """
    Args: receives a filename for a file containing an array of wavelengths and
        corresponding fluxes.
    Returns: Fits a new Bagpipes galaxy object to that data using certain
        assumptions about dust and nebulae emmission, and returns the object.
    """
    # Create empty model components dictionary.
    model_components = {"redshift" : 0.0,
                        "exponential" : {},
                        "dust" : {}
                        }

    if filename[-4:] != ".csv": filename += ".csv"
    # Reassemble the model components dictionary from the header.
    file = open("data/"+filename)
    for line in file.readlines():
        if line[0] == "#":  # The line is a header line, unpack and sort.
            line = line[2:] # Trim header character & whitespace from the line.
            key, value = line.split(":")
            if key == "redshift":
                model_components[key] = float(value)
            if key in ["age", "tau", "tau", "massformed", "metallicity"]:
                model_components["exponential"][key] = float(value)
            if key in ["type", "Av"]:
                try:
                    model_components["dust"][key] = float(value)
                except ValueError:
                    model_components["dust"][key] = value[:-1]

        else: # The line is a data line, exit the loop.
            break

    # Create new galaxy object from the spectrum.
    galaxy = pipes.galaxy(filename, load_data, photometry_exists=False)
    # Return the newly created object and the reconstructed dictionary.
    return(galaxy, model_components)

def import_sfh(filename):
    pass

def main():
    """
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
    """
    # burst = {}
    # burst["age"] = (0., 15.)            # Vary age from 0 to 15 Gyr
    # burst["metallicity"] = (0., 2.5)    # Vary metallicity from 0 to 2.5 Solar
    # burst["massformed"] = (0., 13.)     # Vary log_10(mass formed) from 0 to 13

    sfh_exponential = {"age"         : (0.0, 15.0),   # Gyr
                       "tau"         : (0.0, 2.0),    # Gyr
                       "massformed"  : (0.0, 15.0),   # log_10(M*/M_solar)
                       "metallicity" : (0.0, 2.5)     # Z/Z_oldsolar
                       }

    fit_instructions = {"exponential" : sfh_exponential, # Add the exp SFH component.
                        "redshift"    : (0.0, 10.0)      # Obs. redshift from 0-10.
                        }

    # galaxy, model_components = import_spectrum("model1")
    # galaxy.plot()
    # fit = pipes.fit(galaxy, fit_instructions)
    # fit.fit(verbose=True)

    galaxy = pipes.galaxy("20200127_xshoot_corr", load_xshooter, photometry_exists=False)
    galaxy.plot()
    fit = pipes.fit(galaxy, fit_instructions)
    fit.fit(verbose=True)

    # fit.plot_spectrum_posterior()  # Shows the input and fitted spectrum/photometry
    # fit.plot_sfh_posterior()       # Shows the fitted star-formation history
    # fit.plot_1d_posterior()        # Shows 1d posterior probability distributions
    # fit.plot_corner()              # Shows 1d and 2d posterior probability distributions

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
