import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

df = pd.read_csv("time_mother_length_wild_type.csv")
df2 = pd.read_csv('time_flour_wild_type.csv')
observed_times = df["time"]
events = [] 

# find cell division events and create a dataframe for these
cell_data = df.drop(columns=["time"])
for column in cell_data.columns:
    mother_length = cell_data[column]
    birth_size = mother_length[0]

    for i in range(1, len(mother_length)):
        # a cell division is a drop of over 10 
        if mother_length[i] - mother_length[i-1] < -10:
            division_time = observed_times[i]
            added_size = mother_length[i - 1] - birth_size

            events.append({
                "cell": int(column),
                "event_time": division_time,
                "added_size": added_size
            })

            birth_size = mother_length[i] 

cell_divisions = pd.DataFrame(events)

# define CIF
def conditional_intensity(t, event_times, added_sizes, alpha, beta, rho, gamma, delta):
    added = 0
    clock = 0
    clock_func = lambda t: np.cos(0.1955 * t + 2.5460) 
    N=0 # count for the number of events before t 
    for i in range(0, len(event_times)):
        if event_times[i] < t:
            added = added_sizes[i]
            clock = clock_func(event_times[i])
            N += 1

    exponent = alpha + beta* ( t - rho * N + gamma*added + delta*clock)
    return np.exp(exponent)

def log_likelihood(params):
    # find log likelihood 
    alpha, beta, rho, gamma, delta = params
    logL = 0

    # summation of the log of CIF term
    for cell in range(1, 64):
        cell_events = cell_divisions[cell_divisions["cell"] == cell]
        T = cell_events["event_time"].max()
        event_times = cell_events["event_time"].values
        added_sizes = cell_events["added_size"].values

        log_intensities = np.log([
           conditional_intensity(t, event_times, added_sizes, alpha, beta, rho, gamma, delta)
           for t in event_times
        ])
        logL += np.sum(log_intensities)

        # use numpy.trapezoidal to find the integral term 
        t_grid = np.linspace(0, T, 500)
        intensity_vals = [conditional_intensity(t, event_times, added_sizes, alpha, beta, rho, gamma, delta) for t in t_grid]
        integral = np.trapz(intensity_vals, t_grid)
        logL -= integral

    # return the negative to minimise the negative function (maximise the log likelihood)
    return -logL

# positivity constraints on parameters else CIF becomes an overall decreasing function
bounds = [(None, None), 
          (0.00001, None),   
          (0.00001, None),  
          (0.00001, None),
          (0.00001, None)] 

initial_params = [0.1, 0.1, 0.1, 0.1, 0.1]  # initial parameter guesses for alpha, beta, rho, gamma, delta 

# optimise the function to obtain our parameter estimates 
parameters = scipy.optimize.minimize(log_likelihood, initial_params, method='L-BFGS-B', bounds=bounds) 
print(parameters.x)
