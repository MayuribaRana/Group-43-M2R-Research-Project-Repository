import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("time_mother_length_wild_type.csv")
df2 = pd.read_csv('time_flour_wild_type.csv')
observed_times = df["time"]
birth_lengths = []
events = [] 

cell_data = df.drop(columns=["time"])
for column in cell_data.columns:
    mother_length = cell_data[column]
    birth_size = mother_length[0]
    birth_lengths.append(birth_size)

    for i in range(1, len(mother_length)):
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

def conditional_intensity(t, event_times, added_sizes, alpha, beta, rho, gamma, delta):
    summation_added = 0
    summation_clock = 0
    
    clock_func = lambda t: np.cos(0.1955 * t + 2.5460) 
    N=0
    for i in range(0, len(event_times)):
        if event_times[i] < t:
            summation_added = added_sizes[i]
            summation_clock = clock_func(event_times[i])
            N += 1

    exponent = alpha + beta* ( t - rho * N + gamma*summation_added + delta*summation_clock)
    return np.exp(exponent)

def hazard(t, events, added):
    return conditional_intensity(t, events, added, -2.58024047e+00,  5.92850811e-02,  1.10708964e+01,  1.00000000e-05, 8.47101935e-01)

def thinning(T, cell_max, L_birth):
    t = 0
    S = [0]
    added_sizes = [L_birth]
    while t < T:
        m = 5
        l = 20
        s = np.random.exponential(m * cell_max)
        u = np.random.uniform(0, 1)
        added_sizes.append(L_birth * np.exp(0.0628*(t - S[-1])))
        h = hazard(t+s, S, added_sizes)
        if s > l:
            t = t+l
        elif t + s > T or u > h/m:
            t = t+s
        else:
            t = t+s
            S.append(t)
    return S

def cell_info(cell):
    cell_events = cell_divisions[cell_divisions["cell"] == cell]
    L_birth = birth_lengths[cell - 1]
    event_times = cell_events["event_time"].values
    added_sizes = cell_events["added_size"].values
    T = cell_events["event_time"].max()
    lambda_max = conditional_intensity(T, event_times, added_sizes, -2.58024047e+00,  5.92850811e-02,  1.10708964e+01,  1.00000000e-05, 8.47101935e-01)
    
    return L_birth, lambda_max


inter_divisions = []
for i in range(1, 65):
    birth, lambda_max = cell_info(i)
    division_timings  = thinning(500, lambda_max, birth)
    for j in range(20, len(division_timings)):
        inter_divisions.append(division_timings[j] - division_timings[j-1])

m, s = np.mean(inter_divisions), np.std(inter_divisions)
sns.histplot(inter_divisions, binwidth = 0.5, kde= True, color= 'black', edgecolor = 'red')
plt.text(0.95, 0.95, f"Mean = {m:.4f}\nStd = {s:.4f}", 
         horizontalalignment = "right", 
         verticalalignment = "top", 
         transform = plt.gca().transAxes, 
         fontsize = 11)
plt.legend([],[], frameon=False)
plt.show()