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
    spectrum_with_errs = np.loadtxt("data/mods/" + filename + ".csv", delimiter=" ")
    return(spectrum_with_errs)

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

# def load_xshooter_spec(ID):
#     data = np.loadtxt("data/20200127_xshoot_corr.asci", dtype="float")
#     lambdas, fluxes = data[:,0], data[:,1]

def load_xshooter(ID):
    """
    Data import function for the __init__ method of bagpipes.galaxy. Loads in
    wavelengths, fluxes, and flux errors from the  XSHOOTER datafiles.
    """
    spectrum_with_errs = np.loadtxt("data/tdes/"+ID+".txt", dtype="float")
    return(spectrum_with_errs)

def import_spectrum(filename):
    """
    Args: receives a filename for a file containing an array of wavelengths and
        corresponding fluxes.
    Returns: Fits a new Bagpipes galaxy object to that data using certain
        assumptions about dust and nebulae emmission, and returns the object.
    """
    # Create empty model components dictionary.
    model_components = {}
    # Reassemble the model components dictionary from the header.
    file = open("data/mods/" + filename + ".csv")
    for line in file.readlines():
        if line[0] == "#":  # The line is a header line, unpack and sort.
            line = line[2:-1] # Trim header character & newline from the line.
            values = line.split(":")
            # Create top-level key-value pair.
            if len(values) == 2:
                try:
                    model_components[values[0]] = float(values[1])
                except ValueError:
                    model_components[values[0]] = values[1]
            # Add key-value pair to a sub-dictionary if it exists, otherwise create the sub-dicitonary.
            if len(values) == 3:
                if values[0] in list(model_components.keys()):
                    model_components[values[0]][values[1]] = float(values[2])
                else:
                    try:
                        model_components[values[0]] = {values[1] : float(values[2])}
                    except ValueError:
                        model_components[values[0]] = {values[1] : values[2]}


            # if components[0] == "redshift":
            #     model_components[components[0]] = float(components[1])
            # if key in ["age", "tau", "massformed", "metallicity"]:
            #     model_components["exponential"][key] = float(value)
            # if key in ["type", "Av"]:
            #     try:
            #         model_components["dust"][key] = float(value)
            #     except ValueError:
            #         model_components["dust"][key] = value[:-1]

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

def chi_squared(galaxy, fit):
    fit.posterior.get_advanced_quantities()
    # Calculate the median posterior spectrum.
    spec = fit.posterior.samples["spectrum"]
    # TODO: these keys don't seem to exist for fit.posterior.samples
    # spec /= fit.posterior.samples["calib"]
    # spec += fit.posterior.samples["noise"]
    posterior_spectrum = np.percentile(spec, 50, axis=0)
    # Find the wavelength interval matching the a priori spectrum.
    posterior_wavs = fit.galaxy.spectrum[:,0]
    # prior_wavs = galaxy.spectrum[:,0]
    # for i in range(len(posterior_wavs)):
    #     if np.all(prior_wavs == posterior_wavs[i:i+len(prior_wavs)]):
    #         indices = np.arange(i, i+len(prior_wavs), 1, dtype=int)
    #         break
    # Calculate chi_squared for flux over the a priori wavelength range.
    chi_squared = np.sum((galaxy.spectrum[:,1]-posterior_spectrum[:])**2 / galaxy.spectrum[:,2]**2)
    return(chi_squared)

def print_posterior(fit):
    print ('parameter     median     16th percentile     84th percentile')
    for key in fit.posterior.samples.keys():
        print(key+": ", np.median(fit.posterior.samples[key]), np.percentile(fit.posterior.samples[key], 16), np.percentile(fit.posterior.samples[key], 84))



def plot_corner(fit, names=[], show=False, save=True, bins=25, type="fit_params"):
    """ Make a corner plot of the fitted parameters. """

    import corner

    tex_on = True

    update_rcParams()

    if names == []:
        names = fit.fitted_model.params
        samples = np.copy(fit.posterior.samples2d)
    else:
        for name in names:
            index = fit.fitted_model.params.index(name)
            column = np.array([fit.posterior.samples2d[:,index]]).T
            try:
                samples = np.concatenate((samples, column), axis=1)
            except UnboundLocalError:
                samples = column

    # Set up axis labels
    if tex_on:
        labels = fix_param_names(names)
    else:
        labels = fit.fitted_model.params

    # Log any parameters with log_10 priors to make them easier to see
    for name in names:
        i = fit.fitted_model.params.index(name) # Index among all samples.
        j = names.index(name)                   # Index among plot samples.
        if fit.fitted_model.pdfs[i] == "log_10":
            samples[:, j] = np.log10(samples[:, j])
            if tex_on:
                labels[j] = "$\\mathrm{log_{10}}(" + labels[j][1:-1] + ")$"
            else:
                labels[j] = "log_10(" + labels[j] + ")"

    # Make the corner plot
    fig = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
                        show_titles=True, title_kwargs={"fontsize": 13},
                        smooth=1., smooth1d=1., bins=bins)

    # Save the corner plot to file
    if save:
        plotpath = ("pipes/plots/" + fit.run + "/" + fit.galaxy.ID + "_corner.jpg")
        plt.savefig(plotpath, bbox_inches="tight")
        plt.close(fig)

    # Alternatively show the corner plot
    if show:
        plt.show()
        plt.close(fig)

    return fig

def fix_param_names(fit_params):

    latex_names = {"redshift": "z",
    "metallicity": "Z",
    "massformed": "\\mathrm{log_{10}(M",
    "mass": "\\mathrm{log_{10}(M_*",
    "stellar_mass": "\\mathrm{log_{10}(M_*",
    "tau": "\\tau",
    "alpha": "\\alpha",
    "beta": "\\beta",
    "age": "\\mathrm{Age}",
    "age_min": "\\mathrm{Min\\ Age}",
    "age_max": "\\mathrm{Max\\ Age}",
    "Av": "{A_V}",
    "n": "n",
    "veldisp": "\\sigma_{vel}",
    "0": "\\mathrm{N}0",
    "1": "\\mathrm{N}1",
    "2": "\\mathrm{N}2",
    "3": "\\mathrm{N}3",
    "4": "\\mathrm{N}4",
    "5": "\\mathrm{N}5",
    "6": "\\mathrm{N}6",
    "7": "\\mathrm{N}7",
    "8": "\\mathrm{N}8",
    "9": "\\mathrm{N}9",
    "10": "\\mathrm{N}10",
    "sfr": "\\mathrm{SFR}",
    "mass_weighted_age": "\\mathrm{Age_{MW}}",
    "tform": "\\mathrm{t_{form}}",
    "tquench": "\\mathrm{t_{quench}}",
    "ssfr": "\\mathrm{log_{10}(sSFR",
    "sig_exp": "\\Delta",
    "prob": "P",
    "mu": "\\mu",
    "sigma": "\\sigma",
    "tau_q": "\\tau_\\mathrm{quench}",
    "length": "l",
    "norm": "n",
    "scaling": "s",
    "t_bc": "t_{BC}",
    "B": "B",
    "delta": "\delta",
    "fwhm": "\\mathrm{FWHM}"}

    latex_units = {"metallicity": "Z_{\\odot}",
    "massformed": "M_{\\odot})}",
    "mass": "M_{\\odot})}",
    "stellar_mass": "M_{\\odot})}",
    "tau": "\\mathrm{Gyr}",
    "age": "\\mathrm{Gyr}",
    "age_min": "\\mathrm{Gyr}",
    "age_max": "\\mathrm{Gyr}",
    "Av": "\\mathrm{mag}",
    "veldisp": "\\mathrm{km s^{-1}}",
    "sfr": "\\mathrm{M_\\odot\\ yr}^{-1}",
    "ssfr": "\\mathrm{yr}^{-1})}",
    "mass_weighted_age": "\\mathrm{Gyr}",
    "tform": "\\mathrm{Gyr}",
    "tau_q": "\\mathrm{Gyr}",
    "tquench": "\\mathrm{Gyr}",
    "t_bc": "\\mathrm{Gyr}",
    "fwhm": "\\mathrm{Gyr}"}

    latex_comps = {"dblplaw": "dpl",
    "exponential": "exp",
    "exponential2": "burst",
    "constant": "const",
    "delayed": "del",
    "calibration": "calib",
    "nebular": "neb",
    "lognormal": "",
    "iyer2019": "GP"}

    new_params = []

    if not isinstance(fit_params, list):
        fit_params = [fit_params]

    for fit_param in fit_params:
        split = fit_param.split(":")

        if len(split) == 1:
            comp = None
            param = split[0]
        if len(split) == 2:
            comp = split[0]
            param = split[1]

        if param in list(latex_names):
            new_param = latex_names[param]

            if comp is not None:
                if comp in list(latex_comps):
                    new_param += "_\\mathrm{" + latex_comps[comp] + "}"
                else:
                    new_param += "_\\mathrm{" + comp + "}"

            if param in list(latex_units):
                new_param = new_param + "/" + latex_units[param]

            new_param = "$" + new_param + "$"

        else:
            new_param = fit_param

        new_params.append(new_param)

    if len(new_params) == 1:
        new_params = new_params[0]

    return new_params



def update_rcParams():

    import matplotlib as mpl

    tex_on = True

    mpl.rcParams["lines.linewidth"] = 2.
    mpl.rcParams["axes.linewidth"] = 1.5
    mpl.rcParams["axes.labelsize"] = 18.
    mpl.rcParams["xtick.top"] = True
    mpl.rcParams["xtick.labelsize"] = 14
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.right"] = True
    mpl.rcParams["ytick.labelsize"] = 14
    mpl.rcParams["ytick.direction"] = "in"

    if tex_on:
        mpl.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
        mpl.rc('text', usetex=True)
        mpl.rcParams["text.usetex"] = True

    else:
        mpl.rcParams["text.usetex"] = False
