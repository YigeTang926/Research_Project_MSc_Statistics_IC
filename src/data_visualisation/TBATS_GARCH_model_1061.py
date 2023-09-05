import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['derived_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'czech')
args['figures_dir'] = os.path.join(args['cwd'], 'figures', 'czech')

import sys
sys.path.append(os.path.join(args['cwd'], 'src', 'data_analysis_czech'))
from find_periods import periodogram, direct, find_peaks

# Read the data
wind_1061 = pd.read_csv(os.path.join(args['derived_dir'], 'wind_time_series.csv'),index_col=0)['1061']
solar_1061 = pd.read_csv(os.path.join(args['derived_dir'], 'solar_time_series.csv'),index_col=0)['1061']

# Center the time series
wind_center = wind_1061 - wind_1061.values.mean(axis=0)
solar_center = solar_1061 - solar_1061.values.mean(axis=0)
# Compute the direct spectral estimator
f_wind, sp_wind = direct(wind_center.values, 0.5, shift=1)
f_solar, sp_solar = direct(solar_center.values, 0.5, shift=1)
# Get the periods for each node
peak_T_wind, peak_f_wind, peak_sp_wind = find_peaks(f_wind[f_wind>0],sp_wind[np.where(f_wind>0)],4,5e5,20)
peak_T_solar, peak_f_solar, peak_sp_solar = find_peaks(f_solar[f_solar>0],sp_solar[np.where(f_solar>0)],4,5e5,20)

# Figure 1: Direct spectral estimator using a 50% cosine taper
fig, axes = plt.subplots(1,2,figsize=(30,12))
axes[0].plot(f_wind,sp_wind,color='b')
axes[0].set_title('Wind',fontsize=40)
axes[0].set_xlabel('f',fontsize=35)
axes[0].set_ylabel(r'$\widehat{S}^{(d)}(f)$',fontsize=35)
axes[0].set_xticks(np.arange(-0.5,0.6,0.5))
axes[0].yaxis.set_major_locator(plt.MaxNLocator(5)) 
axes[0].tick_params(axis='both', which='both', labelsize=30)
for T, f, sp in zip(peak_T_wind, peak_f_wind, peak_sp_wind):
    axes[0].text(f, sp, f'{T}', ha='center', va='bottom', fontsize=18)
axes[1].plot(f_solar,sp_solar,color='b')
axes[1].set_title('Solar',fontsize=40)
axes[1].set_xlabel('f',fontsize=35)
axes[1].set_ylabel(r'$\widehat{S}^{(d)}(f)$',fontsize=35)
axes[1].set_xticks(np.arange(-0.5,0.6,0.5))
axes[1].yaxis.set_major_locator(plt.MaxNLocator(5)) 
axes[1].tick_params(axis='both', which='both', labelsize=30)
for T, f, sp in zip(peak_T_solar, peak_f_solar, peak_sp_solar):
    axes[1].text(f, sp, f'{T}', ha='center', va='bottom', fontsize=18)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'node_1061_spectra.pdf'))


# Plot the seasonal decomposition
wind_components_1061 = pd.read_csv(os.path.join(args['derived_dir'], 'seasonal_decomposition', 'wind_components.csv'),index_col=0).filter(like='1061', axis=1)
solar_components_1061 = pd.read_csv(os.path.join(args['derived_dir'], 'seasonal_decomposition', 'solar_components.csv'),index_col=0).filter(like='1061', axis=1)
# transform the index to datetime
daily_time_format = '%Y-%m-%d'
wind_components_1061.index = pd.to_datetime(wind_1061.index,format=daily_time_format)
solar_components_1061.index = pd.to_datetime(solar_1061.index,format=daily_time_format)
# Calculate the remainder
wind_remainder_1061 = wind_components_1061.iloc[:,0] - wind_components_1061.iloc[:,1:].sum(axis=1)
solar_remainder_1061 = solar_components_1061.iloc[:,0] - solar_components_1061.iloc[:,1:].sum(axis=1)

fig, ax = plt.subplots(wind_components_1061.shape[1]+1, 1, figsize=(12, 15))
for i in range(wind_components_1061.shape[1]):
    wind_components_1061.iloc[:,i].plot(ax=ax[i], color='black', linewidth=0.5)
    ax[i].set_ylabel(wind_components_1061.columns[i].split('.')[1])
    ax[i].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax[i].set_xlabel('')
wind_remainder_1061.plot(ax=ax[i+1], color='black', linewidth=0.5)
ax[i+1].set_ylabel('remainder')
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'],'wind_decomposition_1061.pdf'))

fig, ax = plt.subplots(wind_components_1061.shape[1]+1, 1, figsize=(12, 15))
for i in range(solar_components_1061.shape[1]):
    solar_components_1061.iloc[:,i].plot(ax=ax[i], color='black', linewidth=0.5)
    ax[i].set_ylabel(solar_components_1061.columns[i].split('.')[1])
    ax[i].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax[i].set_xlabel('')
solar_remainder_1061.plot(ax=ax[i+1], color='black', linewidth=0.5)
ax[i+1].set_ylabel('remainder')
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'],'solar_decomposition_1061.pdf'))