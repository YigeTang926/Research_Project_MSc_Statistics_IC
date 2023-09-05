import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('seaborn')

# Get the path
args = {}
args['seed'] = 111
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['simulation_data'] = os.path.join(args['cwd'], 'data', 'derived', 'TBATS_GARCH_simulation')
args['simulation_figure'] = os.path.join(args['cwd'], 'figures', 'TBATS_GARCH_simulation')

# Number of observeds
n = 1000

# Trend
trend = 0.05 * np.arange(n)

# Seasonal
periods = [7, 30.42, 100]
amptitudes = [3, 3, 3]
seasonal_1 = amptitudes[0] * np.cos(2 * np.pi * np.arange(n) / periods[0])
seasonal_2 = amptitudes[1] * np.sin(2 * np.pi * np.arange(n) / periods[1])
seasonal_3 = amptitudes[2] * np.sin(2 * np.pi * np.arange(n) / periods[2] + np.pi / 4)
seasonal = seasonal_1 + seasonal_2 + seasonal_3

# Simulate ARMA(1,1)+GARCH(1,1) process
phi = np.array([0.8]) # ar coefficients
theta = np.array([0.5]) # ma coefficients
omega = 1 # garch intercept
alpha = np.array([0.5]) # arch coefficients
beta = np.array([0.4]) # garch coefficients

def simulate_arma_garch(n, phi, theta, omega, alpha, beta, seed):
    p = len(phi)
    q = len(theta)
    v = len(alpha)
    u = len(beta)
    # Simulate GARCH process
    np.random.seed(seed)
    z = np.random.normal(size=n)
    e = np.zeros(n)
    s = np.zeros(n)
    for i in range(max(u,v), n):
        s[i] = np.sqrt(omega + np.sum(alpha * np.flip(e[(i-v):i])**2) + np.sum(beta * np.flip(s[(i-u):i])**2))
        e[i] = s[i] * z[i]
    # Simulate ARMA process
    x = np.zeros(n)
    for i in range(max(p,q), n):
        x[i] = np.sum(phi * np.flip(x[(i-p):i])) + np.sum(theta * np.flip(e[(i-q):i]))
    return x, e, s

error, epsilon, sigma = simulate_arma_garch(n, phi, theta, omega, alpha, beta, seed=args['seed'])

# Process with trend, seasonality, and ARMA + GARCH errors
observed = trend + seasonal + error

# Save simulated data
df = pd.DataFrame({'observed': observed, 'trend': trend, 'seasonal_1': seasonal_1, 'seasonal_2': seasonal_2, 'seasonal_3': seasonal_3, 'error': error, 'epsilon': epsilon, 'sigma': sigma})
df.to_csv(os.path.join(args['simulation_data'],'TG_one_simulation_data.csv'), index=False)

# Plot simulated data
fig, ax = plt.subplots(6, 1, figsize=(12, 9))

ax[0].plot(observed, color='black', linewidth=0.5)
ax[0].set_ylabel('observed')
ax[0].set_xticks([])
ax[0].set_yticks([0,25,50])

ax[1].plot(trend, color='red', linewidth=0.5)
ax[1].set_ylabel('trend')
ax[1].set_xticks([])
ax[1].set_yticks([0,25,50])

ax[2].plot(seasonal_1, color='blue', linewidth=0.5)
ax[2].set_ylabel('seasonal1')
ax[2].set_xticks([])
ax[2].set_yticks([-amptitudes[0],0,amptitudes[0]])

ax[3].plot(seasonal_2, color='blue', linewidth=0.5)
ax[3].set_ylabel('seasonal2')
ax[3].set_xticks([])
ax[3].set_yticks([-amptitudes[1],0,amptitudes[1]])

ax[4].plot(seasonal_3, color='blue', linewidth=0.5)
ax[4].set_ylabel('seasonal3')
ax[4].set_xticks([])
ax[4].set_yticks([-amptitudes[2],0,amptitudes[2]])

ax[5].plot(error, color='green', linewidth=0.5)
ax[5].set_ylabel('error')
ax[5].set_xlabel('Time')

plt.tight_layout()
plt.savefig(os.path.join(args['simulation_figure'],'TG_one_simulation_data_plot.pdf'))