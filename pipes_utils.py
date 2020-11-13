import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt

def noisify_1(value, percent_error):
    """
    Receives a value and desired percent error on that value. Draws a number
    from a Gaussian distribution around that value with width equal to twice
    the value times the desired percent error. Returns the error.
    """
    sigma = percent_error * value
    error = np.random.normal(loc=value, scale=sigma, size=None)
    return(error)

def noisify_2(values, percent_error):
    """
    Receives an array of values and desired percent error. For each value,
    draws a number from a Gaussian distribution around that value, with
    a width equal to twice the median of the data times the percent error.
    Returns an array of these errors.
    """
    sigma = percent_error * np.median(values)
    errors = np.random.normal(loc=values, scale=sigma, size=values.shape)
    return(errors)

def export_spectrum(filename, model, percent_error):
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
    # Generate noisy errors for the spectrum.
    errors = noisify_1(spectrum[:,1], percent_error)
    # Combine the wavelengths, fluxes, and flux errors into one array.
    spectrum_with_errs = np.concatenate(spectrum, errors, axis=1)
    # Write the data with "secret" model components header to a file.
    if filename[-4:] != ".csv": filename += ".csv"
    np.savetxt("data/"+filename, spectrum_with_errs, delimiter=", ", header=header)

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
    print("This is legacy xshooter!")
    return(spectrum_with_errs)

def load_xshooter_2(ID):
    """
    Data import function for the __init__ method of bagpipes.galaxy. Loads in
    wavelengths and fluxes from the  XSHOOTER .asci files.
    """
    spectrum_with_errs = np.loadtxt("data/"+ID+".txt", dtype="float")
    # spectrum_with_errs = np.zeros((spectrum.shape[0], 3))
    # spectrum_with_errs[:,0:2] = np.copy(spectrum)
    # spectrum_with_errs[:,2] = np.copy(np.sqrt(spectrum[:,1]))
    print("This is the new xshooter!")
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

def print_fit_params(fit):
    labels = ['dust:Av', 'dblplaw:tau', 'dblplaw:alpha', 'dblplaw:beta', 'dblplaw:massformed', 'dblplaw:metallicity', 'redshift']
    print ('parameter     median     16th percentile     84th percentile')
    for i in labels:
        print ( i + ': ', np.median( fit.posterior.samples[i]), np.percentile(fit.posterior.samples[i], 16), np.percentile(fit.posterior.samples[i], 84) )
