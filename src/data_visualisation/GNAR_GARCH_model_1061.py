import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['simulation_data'] = os.path.join(args['cwd'], 'data', 'derived', 'czech', 'gnar_garch')
args['simulation_figure'] = os.path.join(args['cwd'], 'figures', 'czech')

# Read the data
dataframe = pd.read_csv(os.path.join(args['simulation_data'], 'GG_results_1061.csv'),index_col=0)
T = dataframe.shape[0]

# Plot the data
fig, ax = plt.subplots(2, 1, figsize=(12, 8))

ax[0].plot(range(T),dataframe.iloc[:,0], color='black', label='Observations', linewidth=0.5)
ax[0].plot(range(T),dataframe.iloc[:,1], color='red', label='Fitted values', linewidth=0.5)
ax[0].fill_between(range(T),dataframe.iloc[:,1]+1.96*np.sqrt(dataframe.iloc[:,2]), 
                   dataframe.iloc[:,1]-1.96*np.sqrt(dataframe.iloc[:,2]), color='blue', 
                   alpha=0.15, label='95% interval')
ax[0].set_ylabel('')
ax[0].set_title('Wind: GNAR(2,{2,1}) + GARCH(1,1)', fontsize=15)
ax[0].legend(loc = 'lower right', fontsize=12)

ax[1].plot(range(T),dataframe.iloc[:,3], color='black', label='Observations', linewidth=0.5)
ax[1].plot(range(T),dataframe.iloc[:,4], color='red', label='Fitted values', linewidth=0.5)
ax[1].fill_between(range(T),dataframe.iloc[:,4]+1.96*np.sqrt(dataframe.iloc[:,5]), 
                   dataframe.iloc[:,4]-1.96*np.sqrt(dataframe.iloc[:,5]), color='blue', 
                   alpha=0.15, label='95% interval')
ax[1].set_xticks([0,182,366,547,731,912,1096])
ax[1].set_xticklabels(['2012.1','2012.7','2013.1','2013.7','2014.1','2014.7','2015.1'])
ax[1].set_ylabel('')
ax[1].set_title('Solar: GNAR(2,{2,1}) + GARCH(1,1)', fontsize=15)
ax[1].legend(loc = 'lower right', fontsize=12)

plt.tight_layout()
plt.savefig(os.path.join(args['simulation_figure'],'GNAR_GARCH_model_1061.pdf'))