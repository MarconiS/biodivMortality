import geopandas as gpd
import networkx as nx
import torch
from torch_geometric.data import Data
from torch_geometric.nn import MessagePassing, GATConv
from shapely.geometry import Point
from itertools import combinations
from rtree import index
import networkx as nx
from rtree import index
from geopandas import sjoin
from shapely.geometry import Point
import pandas as pd

gdf = gpd.read_file('data/final_dataset.shp')

# Define a spatial proximity threshold (in meters)
proximity_threshold = 100  # Adjust as needed

# Get the centroids
centroids = gdf.centroid

gdf['Easting'] = centroids.x
gdf['Northing'] = centroids.y

# Create a Rtree index
idx = index.Index()

# Create an empty graph
G = nx.Graph()

# Create nodes for each point in the DataFrame
for count, row in gdf.iterrows():
    node = count  
    lat = row['Easting'] 
    lon = row['Northing'] 
    G.add_node(node, pos=(lon, lat))
    idx.insert(node, (lon, lat, lon, lat))

# Defining proximity_threshold
proximity_threshold = 100

# Create edges based on spatial proximity
for node1 in G.nodes(data=True):
    potential_edges = list(idx.nearest((node1[1]['pos']), 
                num_results=None, objects=False, 
                distance_less_equal=proximity_threshold))

    # Remove self-loop
    potential_edges.remove(node1[0])

    for node2_index in potential_edges:
        node2 = G.nodes(data=True)[node2_index]

        pos1 = node1[1]['pos']
        pos2 = node2['pos']

        # Calculate the distance between points
        distance = Point(pos1).distance(Point(pos2))

        # If the distance is within the proximity threshold, create an edge
        if distance <= proximity_threshold:
            G.add_edge(node1[0], node2_index)