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
args['derived_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'czech')
args['figures_dir'] = os.path.join(args['cwd'], 'figures', 'czech')

# Read the data
wind_czech = pd.read_csv(os.path.join(args['derived_dir'], 'wind_time_series.csv'),index_col=0)
solar_czech = pd.read_csv(os.path.join(args['derived_dir'], 'solar_time_series.csv'),index_col=0)

# Read the network
network_nodes = pd.read_csv(os.path.join(args['raw_dir'], 'network_nodes.csv'))
network_edges = pd.read_csv(os.path.join(args['raw_dir'], 'network_edges.csv'))
# get the ID of nework nodes of country=CZE
selected_nodes = network_nodes.loc[network_nodes['country'] == 'CZE']['ID'].values

# Create a dictionary of node locations with latitude and longitude
nodes_locations = {}
for id in selected_nodes:
    node_row = network_nodes.loc[network_nodes['ID'] == id]
    lon = node_row['longitude'].values[0]
    lat = node_row['latitude'].values[0]
    nodes_locations[id] = (lon,lat)

# Initialize a Graph object and add edges to it
G = nx.DiGraph()
for id in selected_nodes:
    G.add_node(id)
for id in range(len(network_edges)):
    if network_edges['fromNode'][id] in selected_nodes and network_edges['toNode'][id] in selected_nodes:
        G.add_edge(network_edges['fromNode'][id],network_edges['toNode'][id])

# Figure 1: Plot the network
fig = plt.figure(figsize=(8,4.4))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([12, 19, 48.5, 51.2])
ax.set_aspect(aspect=1.4) # Set the aspect ratio of the plot
# Add map features (coastlines, land, oceans, etc.)
ax.add_feature(cfeature.BORDERS)
# Draw latitude and longitude gridlines
ax.gridlines(draw_labels=True, linewidth=1, color='white', zorder=0)
# Draw the nodes and edges
nx.draw_networkx_nodes(
    G,
    pos=nodes_locations,
    node_size=240,
    node_color='white',
    edgecolors='black',
    linewidths=0.5,
    ax=ax,
)
nx.draw_networkx_edges(
    G,
    pos=nodes_locations,
    edge_color='black',
    width=0.5,
    arrowsize=8,
    connectionstyle='arc3,rad=0.2',
    ax=ax,
)
# Add labels to the nodes
adjusted_nodes_locations = {node: (lon, lat) for node, (lon, lat) in nodes_locations.items()}
nx.draw_networkx_labels(G, pos=adjusted_nodes_locations, font_size=7, font_color='blue', ax=ax)
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'nodes_czech.pdf'))

# Figure 2: Plot the correlation matrix as a heatmap
wind_czech.columns = ['W'+str(id) for id in wind_czech.columns]
solar_czech.columns = ['S'+str(id) for id in solar_czech.columns]
wind_solar_czech = pd.concat([wind_czech,solar_czech],axis=1)
correlation_matrix = wind_solar_czech.corr()

plt.figure(figsize=(7.5, 6))
sns.heatmap(correlation_matrix, cmap='viridis')
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'correlation_matrix_czech.pdf'))