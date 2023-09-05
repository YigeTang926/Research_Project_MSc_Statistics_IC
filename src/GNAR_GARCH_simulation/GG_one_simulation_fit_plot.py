import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['simulation_data'] = os.path.join(args['cwd'], 'data', 'derived', 'GNAR_GARCH_simulation')
args['simulation_figure'] = os.path.join(args['cwd'], 'figures', 'GNAR_GARCH_simulation')

# Read the data
observations = pd.read_csv(os.path.join(args['simulation_data'], 'GG_one_simulation_data.csv'),index_col=0)
fitted_values1 = pd.read_csv(os.path.join(args['simulation_data'], 'GG_one_simulation_fittedvalues_gnar.csv'),index_col=0)
fitted_values2 = pd.read_csv(os.path.join(args['simulation_data'], 'GG_one_simulation_fittedvalues_gnar_garch.csv'),index_col=0)
residuals1 = pd.read_csv(os.path.join(args['simulation_data'], 'GG_one_simulation_residuals_gnar.csv'),index_col=0)
conditional_std_node1 = pd.read_csv(os.path.join(args['simulation_data'], 'GG_one_simulation_conditional_std_gnar_garch_node1.csv'),index_col=0)

# Plot the data
fig, ax = plt.subplots(2, 1, figsize=(12, 8))

ax[0].plot(range(1,1001), observations.iloc[:,0], color='black', label='Observations', linewidth=0.5)
ax[0].plot(range(3,1001), fitted_values1.iloc[:,0], color='red', label='Fitted values', linewidth=0.5)
ax[0].fill_between(range(3,1001), fitted_values1.iloc[:,0]+1.96*np.std(residuals1.values), 
                   fitted_values1.iloc[:,0]-1.96*np.std(residuals1.iloc[:,0]), color='blue', 
                   alpha=0.15, label='95% interval')
ax[0].set_ylabel('')
ax[0].set_xticks([])
ax[0].set_ylim([-8,8])
ax[0].set_title('GNAR(2,{2,1})', fontsize=15)
ax[0].legend(loc = 'lower right', fontsize=12)

ax[1].plot(range(1,1001), observations.iloc[:,0], color='black', label='Observations', linewidth=0.5)
ax[1].plot(range(3,1001), fitted_values2.iloc[:,0], color='red', label='Fitted values', linewidth=0.5)
ax[1].fill_between(range(3,1001), fitted_values1.iloc[:,0]+1.96*conditional_std_node1.iloc[:,0], 
                   fitted_values1.iloc[:,0]-1.96*conditional_std_node1.iloc[:,0], color='blue', 
                   alpha=0.15, label='95% interval')
ax[1].set_ylabel('')
ax[1].set_ylim([-8,8])
ax[1].set_title('GNAR(2,{2,1}) + GARCH(1,1)', fontsize=15)
ax[1].legend(loc = 'lower right', fontsize=12)

plt.tight_layout()
plt.savefig(os.path.join(args['simulation_figure'],'GG_one_simulation_fit_plot.pdf'))