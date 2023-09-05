import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['derived_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'production_time_series')
args['figures_dir'] = os.path.join(args['cwd'], 'figures', 'production_time_series')

# Read the data
average_hourly = pd.read_csv(os.path.join(args['derived_dir'], 'average_production_hourly.csv'),index_col=0)
average_daily = pd.read_csv(os.path.join(args['derived_dir'], 'average_production_daily.csv'),index_col=0)
average_wind_by_country = pd.read_csv(os.path.join(args['derived_dir'], 'wind_average_by_country.csv'),index_col=0)
average_solar_by_country = pd.read_csv(os.path.join(args['derived_dir'], 'solar_average_by_country.csv'),index_col=0)

# transform the index to datetime
daily_time_format = '%Y-%m-%d'
average_daily.index = pd.to_datetime(average_daily.index,format=daily_time_format)

# Get the rolling mean
average_daily_rollmean_wind = average_daily['wind'].rolling(window=30).mean()
average_daily_rollmean_solar = average_daily['solar'].rolling(window=30).mean()
# Get the rolling standard deviation
average_daily_rollstd_wind = average_daily['wind'].rolling(window=30).std()
average_daily_rollstd_solar = average_daily['solar'].rolling(window=30).std()

# Set the color palette
sns_blue = '#1f77b4'
sns_orange = '#ff7f0e'

# Figure 1.1: Plot the wind production for each day
fig, ax = plt.subplots(nrows=4, ncols=3, figsize=[15, 15], sharey=True)
ax = ax.flatten()
MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct','Nov','Dec']
for ix, month in enumerate(MONTHS):
    # Plot individual ts
    daily_ts = []
    for _, ts in (
        average_hourly[['wind', 'day_of_month', 'month']]
        .query(f'month == {ix+1}')
        .groupby('day_of_month')
    ):
        if len(ts) == 72:
            daily_ts.append(ts.reset_index()['wind'])
            ts.reset_index()['wind'].plot(
                alpha=0.3, ax=ax[ix], color=sns_blue, label='_no_legend_', fontsize=20
            )
    # Plot the mean ts
    mean_ts = pd.concat(daily_ts, axis=1).mean(axis=1)
    mean_ts.plot(
        ax=ax[ix], color='blue', label='mean', legend=True, linewidth=2
    )
    ax[ix].legend(loc='upper left', fontsize=25, frameon=False)
    # Set the xticks and labels
    ax[ix].set_xticks(np.arange(0, len(mean_ts) + 1, 8),labels=[0,8,16,24,8,16,24,8,16,24],fontsize=20)
    props = dict(boxstyle='round', facecolor='white', alpha=0.1)
    ax[ix].text(0.12, 0.1, '2012', transform=ax[ix].transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax[ix].text(0.43, 0.1, '2013', transform=ax[ix].transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax[ix].text(0.74, 0.1, '2014', transform=ax[ix].transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax[ix].set_title(month, fontsize=30)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'wind_daily_patterns_for_each_month.pdf'))

# Figure 1.2: Plot the solar production for each day
fig, ax = plt.subplots(nrows=4, ncols=3, figsize=[15, 15], sharey=True)
ax = ax.flatten()
MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct','Nov','Dec']
for ix, month in enumerate(MONTHS):
    # Plot individual ts
    daily_ts = []
    for _, ts in (
        average_hourly[['solar', 'day_of_month', 'month']]
        .query(f'month == {ix+1}')
        .groupby('day_of_month')
    ):
        if len(ts) == 72:
            daily_ts.append(ts.reset_index()['solar'])
            ts.reset_index()['solar'].plot(
                alpha=0.3, ax=ax[ix], color=sns_orange, label='_no_legend_', fontsize=20
            )
    # Plot the mean ts
    mean_ts = pd.concat(daily_ts, axis=1).mean(axis=1)
    mean_ts.plot(
        ax=ax[ix], color='red', label='mean', legend=True, linewidth=2
    )
    ax[ix].legend(loc='upper left', fontsize=25, frameon=False)
    # Set the xticks and labels
    ax[ix].set_xticks(np.arange(0, len(mean_ts) + 1, 8),labels=[0,8,16,24,8,16,24,8,16,24],fontsize=20)
    props = dict(boxstyle='round', facecolor='white', alpha=0.1)
    ax[ix].text(0.12, 0.1, '2012', transform=ax[ix].transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax[ix].text(0.43, 0.1, '2013', transform=ax[ix].transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax[ix].text(0.74, 0.1, '2014', transform=ax[ix].transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax[ix].set_title(month, fontsize=30)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'solar_daily_patterns_for_each_month.pdf'))

## Figure 2: Plot the daily wind/solar production time series average over all nodes
fig, axes = plt.subplots(2,1,figsize=(10, 8))
average_daily['wind'].plot(ax=axes[0],label='Original Time Series (daily)',color=sns_blue,fontsize=12)
average_daily_rollmean_wind.plot(ax=axes[0],color='black',label='Rolling Mean (30 days)')
average_daily_rollstd_wind.plot(ax=axes[0],color='orange',label='Rolling Standard Deviation (30 days)')
axes[0].set_ylabel('MWh',fontsize=15)
axes[0].set_xlabel('')
axes[0].set_title('Wind',fontsize=15)
axes[0].legend(fontsize=10,loc='upper right')
average_daily['solar'].plot(ax=axes[1],label='Original Time Series (daily)',color=sns_orange,fontsize=12)
average_daily_rollmean_solar.plot(ax=axes[1],color='black',label='Rolling Mean (30 days)')
average_daily_rollstd_solar.plot(ax=axes[1],color='orange',label='Rolling Standard Deviation (30 days)')
axes[1].set_ylabel('MWh',fontsize=15)
axes[1].set_xlabel('')
axes[1].set_title('Solar',fontsize=15)
axes[1].legend(fontsize=10,loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'average_production_daily.pdf'))

# Figure 3: Plot the daily wind/solar production time series average over all nodes for each country
fig, ax = plt.subplots(5, 5, figsize=(20, 20))
for i in range(5):
    for j in range(5):
        average_wind_by_country.iloc[:,i*5+j].plot(ax=ax[i,j], label='wind', fontsize=12, color=sns_blue)
        average_solar_by_country.iloc[:,i*5+j].plot(ax=ax[i,j], label='solar', fontsize=12, color=sns_orange)
        ax[i,j].set_title(average_wind_by_country.columns[i*5+j], fontsize=24)
        ax[i,j].set_xlabel('')
        ax[i,j].set_xticks([182, 547, 912])
        ax[i,j].set_xticklabels(average_wind_by_country.index[[182, 547, 912]])
        ax[i,j].legend(fontsize=15, loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'average_production_daily_by_country.pdf'))

## Figure 4: Seasonal boxplots of the daily wind/solar production average over all nodes
seasons = ['Spring','Summer','Autumn','Winter']
daily_time_locations = []
daily_time_names = []
for i in range(2012,2015):
    for j in range(3,13,3):
        time = str(i) + '-' + '{:02d}'.format(j) + '-01'
        location = average_daily.index.get_loc(time)
        daily_time_locations.append(location)
        daily_time_names.append(str(i)+'-'+'{:02d}'.format(j))

daily_wind_seasons = []
daily_solar_seasons = []
for i in range(len(daily_time_locations)-1):
    daily_wind_seasons.append(average_daily['wind'].values[daily_time_locations[i]:daily_time_locations[i+1]])
    daily_solar_seasons.append(average_daily['solar'].values[daily_time_locations[i]:daily_time_locations[i+1]])

fig, axes = plt.subplots(2,1,figsize=(12, 10))
labels = ['2012 '+i for i in seasons]+['2013 '+i for i in seasons]+['2014 '+i for i in seasons[:-1]]
bp1 = axes[0].boxplot(daily_wind_seasons,labels=labels,
                      medianprops={'color':'red','linewidth':'1.5'},
                      meanline=True,
                      showmeans=True,
                      meanprops={'color':'blue','ls':'--','linewidth':'1.5'})
axes[0].legend(handles=[bp1['medians'][0],bp1['means'][0]], labels=['median','mean'],fontsize=12)
axes[0].set_ylabel('MWh',fontsize=12)
axes[0].set_title('Wind',fontsize=15)
bp2 = axes[1].boxplot(daily_solar_seasons,labels=labels,
                      medianprops={'color':'red','linewidth':'1.5'},
                      meanline=True,
                      showmeans=True,
                      meanprops={'color':'blue','ls':'--','linewidth':'1.5'})
axes[1].legend(handles=[bp1['medians'][0],bp1['means'][0]], labels=['median','mean'],fontsize=12)
axes[1].set_ylabel('MWh',fontsize=12)
axes[1].set_title('Solar',fontsize=15)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'seasonal_boxplots.pdf'))