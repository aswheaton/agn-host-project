import bagpipes as pipes
import numpy as np

def export_spectrum(filename, model):
    """
    Args: receives a Bagpipes model_galaxy object and writes its spectrum to a
        text file with a header containing the "secret" parameters used to
        generate that spectrum.
    Returns: None
    """
    # Generate the "secret" parameters file header.
    header = ""
    for key, value in model.model_comp.items():
        try: # The value is another dictionary.
            for subkey, subvalue in model.model_comp[key].items():
                print(model.model_comp[key][subkey])
                header += key + ":" + str(model.model_comp[key][subkey])
                header += "\n"
        except AttributeError: # The value is not another dictionary.
            print(model.model_comp[key])
            header += key + ":" + str(model.model_comp[key])
            header += "\n"
    # Extract the wavelength and flux value of full spectrum from the object.
    lamdas = model.spectrum_full
    fluxes = model.fluxes_full
    # Write the data with parameters header to a file.
    if filename[-4:] != ".csv": filename += ".csv"
    array = np.vstack(lambdas, fluxes)
    np.savetxt("data/"+filename, array, delimiter=", ", header=header)

def export_sfh(filename, model):
    pass

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
        if line[0] == "#": # The line is a header line, unpack and sort.
            pass
        else: # The line is a data line.
            break
    return(lambdas, fluxes, model_components)

def import_sfh(filename):

def main():

    exp = {}                          # Tau model star formation history
    exp["age"] = 3.                   # Gyr
    exp["tau"] = 0.75                 # Gyr
    exp["massformed"] = 9.            # log_10(M*/M_solar)
    exp["metallicity"] = 0.5          # Z/Z_oldsolar

    dust = {}                         # Dust component
    dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
    dust["Av"] = 0.2                  # magnitudes

    model_components = {}                   # The model components dictionary\n",
    model_components["redshift"] = 1.0      # Observed redshift  \n",
    model_components["exponential"] = exp
    model_components["dust"] = dust

    goodss_filt_list = np.loadtxt("examples/filters/goodss_filt_list.txt", dtype="str")
    model = pipes.model_galaxy(model_components, filt_list=goodss_filt_list)

    export_spectrum("data/model1.csv", model)

main()
