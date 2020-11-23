import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt

def redden(redshift, spectrum):
    spectrum[:,0] *= (1.0 + redshift)
    return(spectrum)

def deredden(redshift, spectrum):
    spectrum[:,0] /= (1.0 + redshift)
    return(spectrum)

def noisify(values, percent_error):
    """
    Receives an array of values and desired percent error. For each value,
    draws a number from a Gaussian distribution around that value, with
    a width equal to twice the median of the data times the percent error.
    Returns an array of these errors.
    """
    sigma = percent_error * np.median(values)
    errors = np.abs(np.random.normal(loc=0.0, scale=sigma, size=(values.shape[0],1)))
    return(errors)

def export_spectrum(filename, model, percent_error):
    """
    Receives a Bagpipes model_galaxy object and writes its spectrum to a text
    file with a header containing the "secret" parameters used to generate that
    spectrum.
    Args:
        filename (str): filename to write (can end in .csv or no file extension)
        model (model_galaxy): a Bagpipes model_galaxy object
        percent_error (float): the percent error, as a decimal, of the median
        flux value around which to draw noise values for the spectrum.
    Returns: None
    """
    # Generate the "secret" model components file header.
    header = ""
    for key, value in model.model_comp.items():
        try: # The value is another dictionary.
            for subkey, subvalue in value.items():
                header += key + ":" + subkey + ":" + str(subvalue) + "\n"
        except AttributeError: # The value is not another dictionary.
            header += key + ":" + str(value) + "\n"
    # Strip the last newline from the header.
    header = header[:-1]

    # Extract the wavelength and flux value of full spectrum from the object.
    spectrum = model.spectrum
    # Generate noisy errors for the spectrum.
    errors = noisify(spectrum[:,1], percent_error)
    # Combine the wavelengths, fluxes, and flux errors into one array.
    spectrum_with_errs = np.append(spectrum, errors, axis=1)
    # Write the data with "secret" model components header to a file.
    if filename[-4:] != ".csv": filename += ".csv"
    np.savetxt("data/"+filename, spectrum_with_errs, delimiter=" ", header=header)

def load_data(filename):
    """
    Data import function for the __init__() method of bagpipes.galaxy.
    """
    # Load the data and extract into wavelength and flux arrays.
    spectrum_with_errs = np.loadtxt("data/" + filename, delimiter=" ")
    return(spectrum_with_errs)

def load_xshooter(ID):
    """
    Data import function for the __init__ method of bagpipes.galaxy. Loads in
    wavelengths and fluxes from the  XSHOOTER .asci files.
    """
    spectrum_with_errs = np.loadtxt("data/"+ID+".txt", dtype="float")
    spectrum_with_errs = deredden(0.046, spectrum_with_errs)
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

def export_sfh(filename, model):
    pass

def import_sfh(filename):
    pass

def print_fit_params(fit):
    labels = ['dust:Av', 'dblplaw:tau', 'dblplaw:alpha', 'dblplaw:beta', 'dblplaw:massformed', 'dblplaw:metallicity', 'redshift']
    print ('parameter     median     16th percentile     84th percentile')
    for i in labels:
        print ( i + ': ', np.median( fit.posterior.samples[i]), np.percentile(fit.posterior.samples[i], 16), np.percentile(fit.posterior.samples[i], 84) )
