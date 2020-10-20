import bagpipes as pipes
import numpy as np

def export_spectrum(filename, model):
    """
    Args: receives a Bagpipes model object and writes its spectrum to a text
        file with a header containing the "secret" parameters used to generate
        that spectrum.
    Returns: None
    """
    # Generate the "secret" parameters file header.
    header = ""
    for key in model.components.items():
        try: # The value is another dictionary.
            for subkey in model.components.items().items():
                header.append(key+":"+model.components.items()[subkey])
                header.append("\n")
        except AttributeError: # The value is not another dictionary.
            header.append(key+":"+model.components[key])
            header.append("\n")
    # Extract the wavelength and flux value of full spectrum from the object.
    lamdas = model.spectrum_full
    fluxes = model.fluxes_full
    # Write the data with parameters header to a file.
    if filename[-4:] != ".csv": filename.append(".csv")
    array = np.vstack(lambdas, fluxes)
    np.savetxt("data/"+filename, array, delimiter=", ", header=header)

def import_spectrum(filename):
    """
    Args: receives a filename for a file containing an array of wavelengths and
        corresponding fluxes.
    Returns: Fits a new Bagpipes galaxy object to that data using certain
        assumptions about dust and nebulae emmission, and returns the object.
    """
    # Load the data and extract into wavelength and flux arrays.


    # Reassemble the model components dictionary.
    file = load(filename)
    for line in file.readlines():
        if line[0] = "#": # The line is a header line, unpack and sort.
            pass
        else: # The line is a data line.
            break
    return(lambdas, fluxes, model_components)



# exp = {}                          # Tau model star formation history
# exp["age"] = 3.                   # Gyr
# exp["tau"] = 0.75                 # Gyr
# exp["massformed"] = 9.            # log_10(M*/M_solar)
# exp["metallicity"] = 0.5          # Z/Z_oldsolar
#
# dust = {}                         # Dust component
# dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
# dust["Av"] = 0.2                  # magnitudes
#
# model_components = {}                   # The model components dictionary\n",
# model_components["redshift"] = 1.0      # Observed redshift  \n",
# model_components["exponential"] = exp
# model_components["dust"] = dust
