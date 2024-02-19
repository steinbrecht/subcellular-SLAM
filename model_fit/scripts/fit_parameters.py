# the parameters for the fit procedure are defined here

# --- PACKAGES ----------------------
from pathlib import Path
import numpy as np
rng = np.random.default_rng()


# --- fixed parameters for fit --------------------------
# standard fitting methods: 'leastsq' (LM), 'least_squares' (LM with reflecive bounds)
# fitting method (others: 'nelder', 'dual_annealing', 'ampgo', 'emcee')
method = 'least_squares'
# number of times the fit is run
samples = 200

# delay until labeling is incorparated and transciption starts (exp. observation)
overall_delay = 5

# make lower and upper bounds for parameter variation in fit
lower_bound, upper_bound = np.log(1e-04), np.log(2)
lower_bound_decay = np.log(1e-06)
lb_erate, ub_erate = np.log(0.1), np.log(10)

# choose lower and upper bounds for random starting values
lower_start, upper_start = np.log(1e-03), np.log(1.1)
lower_start_decay = np.log(1e-05)
ls_erate, us_erate = np.log(.5), np.log(7.)

# offset for boundary cost, so that value does not go to infinity
offset = 0.05

# -- function that calculates random starting values ---
def random_parameter(lower_start, upper_start):
    x = lower_start + rng.random() * (upper_start - lower_start)
    return x