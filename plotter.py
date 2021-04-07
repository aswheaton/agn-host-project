import os
import bagpipes as pipes
from pipes_utils import *

# oldfiles = ["phil_model_01", "phil_model_02",
#             "phil_model_03", "phil_model_04",
#             "phil_model_05", "phil_model_06",
#             "phil_model_07", "phil_model_08",
#             "phil_model_09", "phil_model_10"
#             ]
#
# datafiles = ["phil_model_1", "phil_model_2",
#              "phil_model_3", "phil_model_4",
#              "phil_model_5", "phil_model_6",
#              "phil_model_7", "phil_model_8",
#              "phil_model_9", "phil_model_10"
#              ]

# runs = ["r1_delayed_noburst","r1_exponential_noburst","r1_lognormal_noburst"]
# runs = ["r1_delayed_burst","r1_exponential_burst","r1_lognormal_burst"]
# runs = ["r2_dblplaw_burst_wide","r2_delayed_burst_wide", "r2_exponential_burst_wide", "r2_lognormal_burst_wide"]
# runs = ["r3_dblplaw_burst_wide_zconst", "r3_exponential_burst_wide_zconst"]

# for oldname, newname in zip(oldfiles, datafiles):
#     for run in runs:
#         try:
#             oldpath = "pipes/posterior/" + run + "/" + oldname + ".h5"
#             newpath = "pipes/posterior/" + run + "/" + newname + ".h5"
#             os.rename(oldpath, newpath)
#         except FileNotFoundError:
#             pass
#
#         if os.path.isfile(newpath):
#             galaxy, model_components = import_spectrum(newname)
#             fit_instructions = {}
#             fit = pipes.fit(galaxy, fit_instructions, run=run)
#             fit.fit(verbose=False)
#             fit.plot_sfh_posterior()
#         else:
#             pass

# oldfiles = ["model001",
#             "model002",
#             "model003",
#             "model004",
#             "model005",
#             "model006",
#             "model007",
#             "model008",
#             "model009",
#             "model010",
#             "model011",
#             "model012",
#             "model013",
#             "model014",
#             "model015",
#             "model016",
#             "model017",
#             "model018",
#             "model019",
#             "model020"
#             ]
#
# newfiles = ["alex_model_1",
#             "alex_model_2",
#             "alex_model_3",
#             "alex_model_4",
#             "alex_model_5",
#             "alex_model_6",
#             "alex_model_7",
#             "alex_model_8",
#             "alex_model_9",
#             "alex_model_10",
#             "alex_model_11",
#             "alex_model_12",
#             "alex_model_13",
#             "alex_model_14",
#             "alex_model_15",
#             "alex_model_16",
#             "alex_model_17",
#             "alex_model_18",
#             "alex_model_19",
#             "alex_model_20"
#             ]
#
# runs = ["r0_priors"]
#
# for oldname, newname in zip(oldfiles, newfiles):
#     for run in runs:
#         try:
#             oldpath = "pipes/plots/" + run + "/" + oldname + "_sfh.pdf"
#             newpath = "pipes/plots/" + run + "/" + newname + "_sfh.pdf"
#             os.rename(oldpath, newpath)
#         except FileNotFoundError:
#             pass


galaxy, model_components = import_spectrum("phil_model_4")
fit_instructions = {}
fit = pipes.fit(galaxy, fit_instructions, run="r4_exponential_burst")
fit.fit(verbose=False)
plot_corner(fit, names=["exponential2:massformed", "exponential1:tau"], save=True)
