import pandas as pd
import numpy as np
import networkx as nx
from pyvis.network import Network

#load and clean data
normazlize_linkage_df = pd.read_csv('normalize_linkage.csv')
normazlize_linkage_df.set_index('Window', inplace=True)

#extracting normalize linkage values, and turn into a list
normazlize_linkage_values = normazlize_linkage_df.values.flatten()
normazlize_linkage_values = normazlize_linkage_values[~np.isnan(normazlize_linkage_values)] # remove NaN values

q3 = np.percentile(normazlize_linkage_values, 75)
print(f'Q3 value: {q3:.4f}') # calculate Q3 value

adjacency_matrix = normazlize_linkage_df.copy()
adjacency_matrix = adjacency_matrix.applymap(lambda x: np.nan if np.isnan(x) else (1 if x > q3 else 0)) # set values greater than q3 to 1, and others to 0
#adjacency_matrix = adjacency_matrix.fillna(np.nan).astype('Int64') # convert to Int64 and keep NaN values

num_windows = adjacency_matrix.shape[0]
degree_centrality = []

for index, row in adjacency_matrix.iterrows():
    degree = np.nansum(row) # count number of 1s in each row ignoring NaN values
    centrality = degree / (num_windows - 1) # calculate degree centrality
    degree_centrality.append((index, centrality)) # store as tuple (window name, centrality)

degree_centrality.sort(key=lambda x: x[1]) # sort by centrality in ascending order

min_centrality = degree_centrality[0]
print(f'Window {min_centrality[0]} has Minimum Degree Centrality: {min_centrality[1]:.4f}')
max_centrality = degree_centrality[-1]
print(f'Window {max_centrality[0]} has Maximum Degree Centrality: {max_centrality[1]:.4f}')
avg_centrality = np.mean([x[1] for x in degree_centrality])
print(f'Average Degree Centrality: {avg_centrality:.4f}')
print('--------------------------------------------------------------------------------------')

for window, centrality in degree_centrality:
    print(f'Window {window} has Degree Centrality: {centrality:.4f}')

#draw graph
G = nx.Graph() # create an empty undirected graph

for node, centrality in degree_centrality:
    G.add_node(node, size=centrality) # add node with centrality for size referecen

for i,row_label in enumerate(adjacency_matrix.index):
    for j,col_label in enumerate(adjacency_matrix.columns):
        if adjacency_matrix.iloc[i,j] == 1:
            G.add_edge(row_label, col_label) # add edge between nodes if there is a connection

# Create a Pyvis network
net = Network(
    #visual size
    height="800px",
    width="100%",
    # background color
    bgcolor="#222222",
    # label color
    font_color="#F4F4F4",
    # allow search by window name
    select_menu=True,
    #filter_menu=True,
    cdn_resources='in_line'  # Allow offline usage
)

#Set the physics layout of the network
# net.barnes_hut(
#     gravity=40000,
#     central_gravity=0.4,
#     spring_length=250,
#     spring_strength=0.001,
#     damping=0.09,
#     overlap=1
# )

# Convert the NetworkX graph to a Pyvis graph
net.from_nx(G)

# Scale node sizes based on centrality
degree_centrality_dict = dict(degree_centrality)
for node in net.nodes:
    centrality = degree_centrality_dict[node['id']]
    node['size'] = 5 + 50 * centrality  # scale size based on centrality
    node['color'] = {
        'background': '#2B7CE9', # node color
        'border': '#2B7CE9',
        'highlight': {
            'background': '#FFA500', #node color when highlighted
            'border': '#FFA500'
        }
    }
    node['borderWidth'] = 2
    node['borderWidthSelected'] = 3

# Configure edges
for edge in net.edges:
    edge['color'] = {'color': '#FFFFFF', 'highlight': '#FFA500'}
    edge['width'] = 0.5
    edge['smooth'] = {'type': 'continuous'} #curve edges

# Additional options json configuration for rendering rule
net.set_options("""
{
  "nodes": {
    "scaling": {
      "min": 10,
      "max": 30
    }
  },
  "edges": {
    "color": {
      "inherit": true
    },
    "smooth": {
      "enabled": true,
      "type": "continuous"
    }
  },
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -80000,
      "centralGravity": 0.3,
      "springLength": 250,
      "springConstant": 0.001,
      "damping": 0.09,
      "avoidOverlap": 0.1
    },
    "minVelocity": 0.75,
    "solver": "barnesHut"
  },
  "interaction": {
    "hover": true,
    "tooltipDelay": 200,
    "hideEdgesOnDrag": false,
    "multiselect": true,
    "navigationButtons": true,
    "keyboard": true
  }
}
""")

# Save the network as an HTML file
net.write_html("interactive_network.html", notebook=False)