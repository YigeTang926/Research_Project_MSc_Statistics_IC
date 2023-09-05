import os
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['simulation_data'] = os.path.join(args['cwd'], 'data', 'derived', 'TBATS_GARCH_simulation')
args['simulation_figure'] = os.path.join(args['cwd'], 'figures', 'TBATS_GARCH_simulation')

# Read the data
components1 = pd.read_csv(os.path.join(args['simulation_data'], 'TG_one_simulation_fit_components1.csv'),index_col=0)
components2 = pd.read_csv(os.path.join(args['simulation_data'], 'TG_one_simulation_fit_components2.csv'),index_col=0)
components3 = pd.read_csv(os.path.join(args['simulation_data'], 'TG_one_simulation_fit_components3.csv'),index_col=0)
components4 = pd.read_csv(os.path.join(args['simulation_data'], 'TG_one_simulation_fit_components4.csv'),index_col=0)

# Plot the data
fig, ax = plt.subplots(5, 1, figsize=(12, 7.5))

ax[0].plot(components1['level'], color='black', label='ARMA(0,0)+GARCH(0,0)', linewidth=0.5)
ax[0].plot(components2['level'], color='blue', label='ARMA(1,1)+GARCH(0,0)', linewidth=0.5)
ax[0].plot(components3['level'], color='green', label='ARMA(0,0)+GARCH(1,1)', linewidth=0.5)
ax[0].plot(components4['level'], color='red', label='ARMA(1,1)+GARCH(1,1)', linewidth=0.5)
ax[0].set_ylabel('level')
ax[0].set_xticks([])
ax[0].legend(loc = 'upper left', fontsize=8)
ax[0].set_yticks([0,25,50])

ax[1].plot(components1['slope'], color='black', linewidth=0.5)
ax[1].plot(components2['slope'], color='blue', linewidth=0.5)
ax[1].plot(components3['slope'], color='green', linewidth=0.5)
ax[1].plot(components4['slope'], color='red', linewidth=0.5)
ax[1].set_ylabel('slope')
ax[1].set_xticks([])
ax[1].set_yticks([-0.55,-0.25,0.05])

ax[2].plot(components1['seasonal1'], color='black', linewidth=0.5)
ax[2].plot(components2['seasonal1'], color='blue', linewidth=0.5)
ax[2].plot(components3['seasonal1'], color='green', linewidth=0.5)
ax[2].plot(components4['seasonal1'], color='red', linewidth=0.5)
ax[2].set_ylabel('seasonal1')
ax[2].set_xticks([])
ax[2].set_yticks([-3,0,3])

ax[3].plot(components1['seasonal2'], color='black', linewidth=0.5)
ax[3].plot(components2['seasonal2'], color='blue', linewidth=0.5)
ax[3].plot(components3['seasonal2'], color='green', linewidth=0.5)
ax[3].plot(components4['seasonal2'], color='red', linewidth=0.5)
ax[3].set_ylabel('seasonal2')
ax[3].set_xticks([])
ax[3].set_yticks([-3,0,3])

ax[4].plot(components1['seasonal3'], color='black', linewidth=0.5)
ax[4].plot(components2['seasonal3'], color='blue', linewidth=0.5)
ax[4].plot(components3['seasonal3'], color='green', linewidth=0.5)
ax[4].plot(components4['seasonal3'], color='red', linewidth=0.5)
ax[4].set_ylabel('seasonal3')
ax[4].set_xlabel('Time')
ax[4].set_yticks([-3,0,3])

plt.tight_layout()
plt.savefig(os.path.join(args['simulation_figure'],'TG_one_simulation_fit_plot.pdf'))