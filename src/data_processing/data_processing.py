import os
import pandas as pd

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['raw_dir'] = os.path.join(args['cwd'], 'data', 'raw')
args['derived_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'production_time_series')
args['czech_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'czech')

# Read the data
network_nodes = pd.read_csv(os.path.join(args['raw_dir'], 'network_nodes.csv'), usecols=['ID','country'])
wind_signal = pd.read_csv(os.path.join(args['raw_dir'], 'wind_signal_COSMO.csv'),index_col=0)
solar_signal = pd.read_csv(os.path.join(args['raw_dir'], 'solar_signal_COSMO.csv'),index_col=0)
wind_layouts = pd.read_csv(os.path.join(args['raw_dir'], 'wind_layouts_COSMO.csv'),index_col=0)
solar_layouts = pd.read_csv(os.path.join(args['raw_dir'], 'solar_layouts_COSMO.csv'),index_col=0)

# fill the NaN values with 0 for solar_signals
solar_signal = solar_signal.fillna(0)

# Multiply the signals with the layouts to get the production time series
wind_production_hourly = wind_signal.multiply(wind_layouts['Proportional'].values.flatten(),axis=1)
solar_production_hourly = solar_signal.multiply(solar_layouts['Proportional'].values.flatten(),axis=1)
# Transform the colnames to integer
wind_production_hourly.columns = wind_production_hourly.columns.astype(int)
solar_production_hourly.columns = solar_production_hourly.columns.astype(int)
# Output the hourly production time series for wind and solar
wind_production_hourly.to_csv(os.path.join(args['derived_dir'], 'wind_production_hourly.csv'))
solar_production_hourly.to_csv(os.path.join(args['derived_dir'], 'solar_production_hourly.csv'))

# Transform the index to datetime
time_format = '%Y-%m-%d %H:%M:%S'
wind_production_hourly.index = pd.to_datetime(wind_production_hourly.index,format=time_format)
solar_production_hourly.index = pd.to_datetime(solar_production_hourly.index,format=time_format)
# Resample to daily
wind_production_daily = wind_production_hourly.resample('D').sum()
solar_production_daily = solar_production_hourly.resample('D').sum()
# Transform the colnames to integer
wind_production_daily.columns = wind_production_daily.columns.astype(int)
solar_production_daily.columns = solar_production_daily.columns.astype(int)
# Output the daily production time series for wind and solar
wind_production_daily.to_csv(os.path.join(args['derived_dir'], 'wind_production_daily.csv'))
solar_production_daily.to_csv(os.path.join(args['derived_dir'], 'solar_production_daily.csv'))

# Average the hourly production over all nodes
average_production_hourly = pd.DataFrame()
average_production_hourly['wind'] = wind_production_hourly.sum(axis=1)/len(wind_production_hourly.columns)
average_production_hourly['solar'] = solar_production_hourly.sum(axis=1)/len(wind_production_hourly.columns)
# Compute date time variables
average_production_hourly['year'] = average_production_hourly.index.year
average_production_hourly['week'] = average_production_hourly.index.isocalendar().week
average_production_hourly['month'] = average_production_hourly.index.month
average_production_hourly['day_of_month'] = average_production_hourly.index.day
# Output the average hourly production time series for wind and solar
average_production_hourly.to_csv(os.path.join(args['derived_dir'], 'average_production_hourly.csv'))

# Average the daily production over all nodes
average_production_daily = pd.DataFrame()
average_production_daily['wind'] = wind_production_daily.sum(axis=1)/len(wind_production_daily.columns)
average_production_daily['solar'] = solar_production_daily.sum(axis=1)/len(wind_production_daily.columns)
# Compute date time variables
average_production_daily['year'] = average_production_daily.index.year
average_production_daily['week'] = average_production_daily.index.isocalendar().week
average_production_daily['month'] = average_production_daily.index.month
average_production_daily['day_of_month'] = average_production_daily.index.day
# Output the average daily production time series for wind and solar
average_production_daily.to_csv(os.path.join(args['derived_dir'], 'average_production_daily.csv'))

# Average the time series over the nodes with the same country
# Melt the dataframe to long format
melted_wind = wind_production_daily.reset_index().melt(id_vars='Time', var_name='ID', value_name='production')
melted_solar = solar_production_daily.reset_index().melt(id_vars='Time', var_name='ID', value_name='production')
# Merge the two dataframes on 'ID'
merged_wind = pd.merge(melted_wind, network_nodes, on='ID')
merged_solar = pd.merge(melted_solar, network_nodes, on='ID')
# Group by 'country' and 'Time', and then take the mean
average_wind_by_country = merged_wind.groupby(['country', 'Time']).mean().reset_index()
average_solar_by_country = merged_solar.groupby(['country', 'Time']).mean().reset_index()
# Pivot the dataframe to get the desired format
average_wind_by_country = average_wind_by_country.pivot(index='Time', columns='country', values='production')
average_solar_by_country = average_solar_by_country.pivot(index='Time', columns='country', values='production')
# Output the average daily production time series for wind and solar
average_wind_by_country.to_csv(os.path.join(args['derived_dir'], 'wind_average_by_country.csv'))
average_solar_by_country.to_csv(os.path.join(args['derived_dir'], 'solar_average_by_country.csv'))

# Get the ID for the nodes in czech
czech_nodes = network_nodes[network_nodes['country']=='CZE']['ID'].values
# Get the production time series for the nodes in czech
czech_wind_production_daily = wind_production_daily.loc[:,czech_nodes]
czech_solar_production_daily = solar_production_daily[czech_nodes]
# Output the daily production time series for wind and solar
czech_wind_production_daily.to_csv(os.path.join(args['czech_dir'], 'wind_time_series.csv'))
czech_solar_production_daily.to_csv(os.path.join(args['czech_dir'], 'solar_time_series.csv'))