import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
from pipes_utils import *

# oldfiles = ["phil_model_01", "phil_model_02","phil_model_03", "phil_model_04","phil_model_05", "phil_model_06","phil_model_07", "phil_model_08","phil_model_09", "phil_model_10"]
datafiles = ["phil_model_1", "phil_model_2", "phil_model_3", "phil_model_4", "phil_model_5", "phil_model_6", "phil_model_7", "phil_model_8", "phil_model_9", "phil_model_10"]

# redshifts = [0.0206, 0.0484, 0.06, 0.026211, 0.022, 0.0512, 0.01513, 0.018]
# datafiles = ["ASASSN14li","ASASSN15oi","AT2018fyk","AT2019ahk","AT2019azh","AT2019dsg", "AT2019qiz", "iPTF16fnl"]

# runs = ["r1_delayed_noburst","r1_exponential_noburst","r1_lognormal_noburst"]
# runs = ["r1_delayed_burst","r1_exponential_burst","r1_lognormal_burst"]
# runs = ["r2_dblplaw_burst_wide","r2_delayed_burst_wide", "r2_exponential_burst_wide", "r2_lognormal_burst_wide"]
# runs = ["r3_dblplaw_burst_wide_zconst", "r3_exponential_burst_wide_zconst"]
runs = ["r4_exponential_burst", "r4_delayed_burst", "r4_lognormal_burst", "r4_dblplaw_burst"]

if True:
    exponential1 = {}
    exponential1["age"] = (7.5, 12.5)   # Gyr
    exponential1["tau"] = (0.5, 2.0)    # Gyr
    exponential1["massformed"] = (5.0, 12.5)   # log_10(M*/M_solar)
    exponential1["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

    exponential2 = {}
    exponential2["age"] = (0.0, 3.5)   # Gyr, lifetime of F type stars
    exponential2["tau"] = (0.1, 2.0)    # Gyr
    exponential2["massformed"] = (0.0, 12.5)   # log_10(M*/M_solar)
    exponential2["metallicity"] = (0.0, 2.5)     # Z/Z_oldsolar

    dblplaw = {}
    dblplaw["tau"] = (0.0, 4.5) # Do not let the peak occur beyond ~6.5 Gyr
    dblplaw["alpha"] = (5.0, 10.0) # Formation must begin and end reasonably fast
    dblplaw["beta"] = (5.0, 10.0) # Formation must begin and end reasonably fast
    dblplaw["massformed"] = (5.0, 12.5)
    dblplaw["metallicity"] = (0.0, 2.5)

    delayed = {}                   # Delayed Tau model t*e^-(t/tau)
    delayed["age"] = (7.5, 12.5)           # Time since SF began: Gyr
    delayed["tau"] = (0.1, 2.0)           # Timescale of decrease: Gyr
    delayed["massformed"] = (5.0, 12.5)
    delayed["metallicity"] = (0.0, 2.5)

    lognormal = {}
    lognormal["tmax"] = (1.0, 6.0)
    lognormal["fwhm"] = (0.1, 3.0)
    lognormal["massformed"] = (5.0, 12.5)
    lognormal["metallicity"] = (0.0, 2.5)

    dust = {}
    dust["type"] = "Calzetti"
    dust["Av"] = (0.0, 2.0)

    nebular = {}
    nebular["logU"] = (-4.0,-2.0) # Prior dictionaries.

def fix_h5_enumeration():
    for oldname, newname in zip(oldfiles, datafiles):
        for run in runs:
            try:
                oldpath = "pipes/posterior/" + run + "/" + oldname + ".h5"
                newpath = "pipes/posterior/" + run + "/" + newname + ".h5"
                os.rename(oldpath, newpath)
            except FileNotFoundError:
                pass

            if os.path.isfile(newpath):
                galaxy, model_components = import_spectrum(newname)
                fit_instructions = {}
                fit = pipes.fit(galaxy, fit_instructions, run=run)
                fit.fit(verbose=False)
                fit.plot_sfh_posterior()
            else:
                pass

def fix_enumeration():
    oldfiles = ["model001","model002","model003","model004","model005","model006","model007","model008","model009","model010","model011","model012","model013","model014","model015","model016","model017","model018","model019","model020"]
    newfiles = ["alex_model_1","alex_model_2","alex_model_3","alex_model_4","alex_model_5","alex_model_6","alex_model_7","alex_model_8","alex_model_9","alex_model_10","alex_model_11","alex_model_12","alex_model_13","alex_model_14","alex_model_15","alex_model_16","alex_model_17","alex_model_18","alex_model_19","alex_model_20"]

    runs = ["r0_priors"]

    for oldname, newname in zip(oldfiles, newfiles):
        for run in runs:
            try:
                oldpath = "pipes/plots/" + run + "/" + oldname + "_sfh.pdf"
                newpath = "pipes/plots/" + run + "/" + newname + "_sfh.pdf"
                os.rename(oldpath, newpath)
            except FileNotFoundError:
                pass

def corner():
    # from matplotlib import rcParams
    # rcParams.update({'figure.autolayout': True})
    galaxy, model_components = import_spectrum("phil_model_4")
    fit_instructions = {}
    fit = pipes.fit(galaxy, fit_instructions, run="r4_exponential_burst")
    fit.fit(verbose=False)
    plot_corner(fit, names=["exponential2:massformed", "exponential1:tau"], save=True)

def get_best_fit():

    for filename in datafiles:

        fit_instructions = {}
        chi_squ_vals = {}

        galaxy, model_components = import_spectrum(filename)
        # Calculate redshift constraints.
        z_low, z_high = model_components["redshift"] - 0.001,  model_components["redshift"] + 0.001

        fit = pipes.fit(galaxy, fit_instructions, run="r4_exponential_burst")
        fit.fit(verbose=False)
        chi_squ_vals = {"r4_exponential_burst" : chi_squared(galaxy, fit)}

        fit = pipes.fit(galaxy, fit_instructions, run="r4_dblplaw_burst")
        fit.fit(verbose=True)
        chi_squ_vals["r4_dblplaw_burst"] = chi_squared(galaxy, fit)

        fit = pipes.fit(galaxy, fit_instructions, run="r4_delayed_burst")
        fit.fit(verbose=False)
        chi_squ_vals["r4_delayed_burst"] = chi_squared(galaxy, fit)

        fit = pipes.fit(galaxy, fit_instructions, run="r4_lognormal_burst")
        fit.fit(verbose=False)
        chi_squ_vals["r4_lognormal_burst"] = chi_squared(galaxy, fit)

        # Get the functional form with the lowest chi-squared value.
        best_func = min(chi_squ_vals, key=lambda k: chi_squ_vals[k])
        # Select the fit with lowest chi-squared value and plot it.
        print(filename)
        print(best_func)
        print(chi_squ_vals[best_func])

def fixed_aspect_ratio(ratio, ax):
    '''
    Set a fixed aspect ratio on matplotlib plots
    regardless of axis units
    '''
    xvals,yvals = ax.get_xlim(), ax.get_ylim()
    xrange = xvals[1]-xvals[0]
    yrange = yvals[1]-yvals[0]
    ax.set_aspect(ratio*(xrange/yrange), adjustable='box')

def tde_histograms():

    update_rcParams()
    nrow, ncol = 2, 4
    row, col = 0, 0
    fig, axes = plt.subplots(nrow, ncol, sharey=True, sharex=True)

    # fig, axes = plt.subplots(nrow, ncol,
    #             gridspec_kw=dict(wspace=0.0, hspace=0.0,
    #             top=1. - 0.5 / (nrow + 1), bottom=0.5 / (nrow + 1),
    #             left=0.5 / (ncol + 1), right=1 - 0.5 / (ncol + 1)),
    #             figsize=(6.8, 3.4), sharey='row', sharex='col',
    #             )
    for file in datafiles:
        spec_with_errs = np.loadtxt("data/tdes/"+file+".txt")
        sn_ratio = abs(spec_with_errs[:,1] / spec_with_errs[:,2])
        axes[row,col].hist(sn_ratio, bins=[0,5,10,15,20,25,30,35,40,45,50,55,60])
        axes[row,col].set_box_aspect(1)
        axes[row,col].set_xticks([5,15,30,45])
        # axes[row,col].set_yticks([2000,4000,6000,8000,10000])
        axes[row,col].set_xlim(0,60)
        axes[row,col].set_ylim(0,25000)
        axes[row,col].annotate(file, (20,20000))

        col += 1
        if col == 4:
            col = 0
            row += 1

    # axes[1, 0].set_xticks([0,10,20,30,40,50,60])
    # fig.set_figwidth(6.75)
    # fig.set_size_inches(6.75, )


    fig.subplots_adjust(wspace=0, hspace=0)
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("$S/N$")
    plt.tight_layout()
    # plt.savefig("report/img/tde_sn_hist.pdf", dpi=300)
    plt.show()

def plot_posterior_quantities(filename):

    fit_instructions = {}
    chi_squ_vals = {}

    galaxy, model_components = import_spectrum(filename)
    # Calculate redshift constraints.
    z_low, z_high = model_components["redshift"] - 0.001,  model_components["redshift"] + 0.001

    fit = pipes.fit(galaxy, fit_instructions, run="r4_exponential_burst")
    fit.fit(verbose=False)
    chi_squ_vals = {"r4_exponential_burst" : chi_squared(galaxy, fit)}

    fit = pipes.fit(galaxy, fit_instructions, run="r4_dblplaw_burst")
    fit.fit(verbose=True)
    chi_squ_vals["r4_dblplaw_burst"] = chi_squared(galaxy, fit)

    fit = pipes.fit(galaxy, fit_instructions, run="r4_delayed_burst")
    fit.fit(verbose=False)
    chi_squ_vals["r4_delayed_burst"] = chi_squared(galaxy, fit)

    fit = pipes.fit(galaxy, fit_instructions, run="r4_lognormal_burst")
    fit.fit(verbose=False)
    chi_squ_vals["r4_lognormal_burst"] = chi_squared(galaxy, fit)

    # Get the functional form with the lowest chi-squared value.
    best_func = min(chi_squ_vals, key=lambda k: chi_squ_vals[k])
    fit = pipes.fit(galaxy, fit_instructions, run=best_func)
    # Select the fit with lowest chi-squared value and plot it.
    print(filename)
    print(best_func)
    for key in chi_squ_vals.keys():
        print(key, chi_squ_vals[key])

    if best_func == "r4_exponential_burst":
        labels = ["exponential2:age", "exponential2:massformed", "exponential1:age", "exponential1:massformed"]
    if best_func == "r4_delayed_burst":
        labels = ["exponential2:age", "exponential2:massformed", "delayed:age", "delayed:massformed"]
    if best_func == "r4_lognormal_burst":
        labels = ["exponential2:age", "exponential2:massformed", "lognormal:tmax", "lognormal:massformed"]
    if best_func == "r4_dblplaw_burst":
        labels = ["exponential2:age", "exponential2:massformed", "dblplaw:tau", "dblplaw:massformed"]

    post_quantities = dict(zip(labels, [fit.posterior.samples[l] for l in labels]))

    tex_on = True
    update_rcParams()

    fig, axes = plt.subplots(1, 4)

    for i in range(4):
        hist1d(post_quantities[labels[i]], axes[i], smooth=True, label=labels[i])
        fixed_aspect_ratio(1,axes[i])

    fig.set_figwidth(6.75)
    plt.savefig("report/img/"+filename+"_posterior.pdf", bbox_inches = 'tight', pad_inches = 0)

for filename in datafiles:
    plot_posterior_quantities(filename)
