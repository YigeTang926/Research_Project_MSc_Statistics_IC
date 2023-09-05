import os
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['raw_dir'] = os.path.join(args['cwd'], 'data', 'raw')
args['figures_dir'] = os.path.join(args['cwd'], 'figures', 'layouts_with_network')

# Read the data of wind/solar layout
wind_layouts = pd.read_csv(os.path.join(args['raw_dir'], 'wind_layouts_COSMO.csv'),index_col=0)
solar_layouts = pd.read_csv(os.path.join(args['raw_dir'], 'solar_layouts_COSMO.csv'),index_col=0)

# Read the network
network_nodes = pd.read_csv(os.path.join(args['raw_dir'], 'network_nodes.csv'))
network_edges = pd.read_csv(os.path.join(args['raw_dir'], 'network_edges.csv'))

# Create a dictionary of node locations with latitude and longitude
nodes_locations = {}
for i in range(len(network_nodes)):
    ID = network_nodes['ID'][i]
    lon = network_nodes['longitude'][i]
    lat = network_nodes['latitude'][i]
    nodes_locations[ID] = (lon,lat)

# Initialize a Graph object and add edges to it
G = nx.Graph()
for i in range(len(network_nodes)):
    G.add_node(network_nodes['ID'][i])
for i in range(len(network_edges)):
    G.add_edge(network_edges['fromNode'][i],network_edges['toNode'][i])

# Create a color palette for the plots
palette = sns.color_palette('tab20')+sns.color_palette('viridis',5)
# color the nodes with countries
countries = network_nodes['country'].unique()
countries = {country:palette[i] for i,country in enumerate(countries)}
node_colors = [countries[network_nodes['country'][i]] for i in range(len(network_nodes))]

## Figure 1.1: Plot the wind layouts on a map
fig = plt.figure(figsize=(12, 7.5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-10, 30, 35, 58])
# Add map features (coastlines, land, oceans, etc.)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAKES, alpha=0.5)
# Get the layouts
layouts = dict(wind_layouts['Proportional'])
# Draw the nodes with sizes proportional to the layouts
nx.draw_networkx_nodes(
    G,
    pos=nodes_locations,
    nodelist=list(layouts.keys()),
    node_size=[x/30 for x in layouts.values()],
    node_color=node_colors,
    edgecolors='k',
    linewidths=0.2,
    ax=ax,
)
nx.draw_networkx_edges(
    G,
    pos=nodes_locations,
    edge_color='k',
    width=0.2,
    ax=ax,
)
plt.title('Wind Layouts in Europe',fontsize=18)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'wind_layouts.pdf'))

## Figure 1.2: Plot the solar layouts on a map
fig = plt.figure(figsize=(12, 7.5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-10, 30, 35, 58])
# Add map features (coastlines, land, oceans, etc.)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAKES, alpha=0.5)
# Get the layouts
layouts = dict(solar_layouts['Proportional'])
# Draw the nodes with sizes proportional to the layouts
nx.draw_networkx_nodes(
    G,
    pos=nodes_locations,
    nodelist=list(layouts.keys()),
    node_size=[x/30 for x in layouts.values()],
    node_color=node_colors,
    edgecolors='k',
    linewidths=0.2,
    ax=ax,
)
nx.draw_networkx_edges(
    G,
    pos=nodes_locations,
    edge_color='k',
    width=0.2,
    ax=ax,
)
plt.title('Solar Layouts in Europe',fontsize=18)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'solar_layouts.pdf'))