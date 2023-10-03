from scipy.spatial import KDTree
from shapely.ops import transform
from functools import partial
import pyproj
from math import sqrt
import networkx as nx
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, Point
import pickle
import rasterio
import rasterstats as rs

def create_grid(minx, miny, maxx, maxy, xsize, ysize):
    cols_list = np.arange(minx, maxx, xsize)
    rows_list = np.arange(miny, maxy, ysize)
    polygons = []
    for x in cols_list:
        for y in rows_list:
            polygons.append(Polygon([(x,y), (x+xsize, y), (x+xsize, y-ysize), (x, y-ysize)]))
    return gpd.GeoDataFrame({'geometry':polygons})

def calc_distance(point1, point2):
    proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'), pyproj.Proj(init='epsg:3857'))
    point1_transformed = transform(proj, point1)
    point2_transformed = transform(proj, point2)
    return point1_transformed.distance(point2_transformed)

def edge_strength(dist, threshold):
    if dist > threshold:
        dist = threshold
    return 1 - (dist ** 2) / (threshold ** 2)

# Define a function to extract the 95th percentile height of each tree
def extract_tree_height(row, chm):
    # Get the geometry of the crown polygon
    geom = row.geometry
    # change nan from chm to 0
    chm_values = chm.read(1)
    chm_values[np.isnan(chm_values)] = 0
    # Extract the values of the CHM raster within the crown polygon
    chm_values = rs.zonal_stats(geom, chm_values, nodata=0, stats='percentile_98', affine=chm.transform)
    return chm_values[0]['percentile_98']

#create a function that finds areas in the grid that have no dead trees, and creates a network centered on a randomly selected tree in that area
def get_networks_with_no_dead(site, threshold_value, threshold_sample):
    gx_df = gpd.read_file(f'data/{site}h.gpkg')  # insert file path to your shapefile
    bounds = gx_df.total_bounds
    grid_gdf = create_grid(bounds[0], bounds[1], bounds[2], bounds[3], threshold_sample, threshold_sample)
    joined = gpd.sjoin(gx_df, grid_gdf, how='left', op='within')
    # select only grids where there are no dead trees
    picked_entries = joined.groupby('index_right').apply(lambda x: x.sample(1) if x.dead_label.sum() == 0 else None)
    # pick only the left index and convert to list
    gap_trees = picked_entries.index.droplevel(0)
    #calculate centroid of each polygon
    gx_df['geometry'] = gx_df['geometry'].centroid
    # from geopandas extract the coordinates of the points
    gx_df['x'] = gx_df['geometry'].apply(lambda p: p.x)
    gx_df['y'] = gx_df['geometry'].apply(lambda p: p.y)
    # Convert coordinates to a list and build spatial index
    points_list = gx_df[['x', 'y']].values.tolist()
    gx_df_spindex = KDTree(points_list)
    networks = {}
    for row in gx_df.itertuples():
        # if row index is in the list of dead trees
        if row.Index in gap_trees:
            point = Point(row.x, row.y)
            indices_within_threshold = gx_df_spindex.query_ball_point((row.x, row.y), threshold_value)
            if len(indices_within_threshold) > 0:
                network = nx.DiGraph()  # Create a new network for each DEAD tree
                for index in indices_within_threshold:
                    neighbor_point = Point(*points_list[index])
                    #calculate distance between point and dead_point
                    dist = point.distance(neighbor_point)
                    strength = edge_strength(dist, threshold_value)
                    network.add_edge(row.Index, index, tag = gx_df.sci_name[index], dist=dist, strength=strength, is_dead = gx_df.dead_label[index],
                                    crown_area = gx_df.crown_area[index], height = gx_df.height[index], east = gx_df.x[index], north = gx_df.y[index])           
                networks[row.Index] = network
    #save networks to file
    for idx, (name, network) in enumerate(networks.items()):
        with open(f'networks_no_dead/{site}/{name}.gpickle', 'wb') as f:
            pickle.dump(network, f, pickle.HIGHEST_PROTOCOL)



def create_networks(site, threshold_value, threshold_sample):
    # Create a new directed graph G
    gx_df = gpd.read_file(f'data/{site}h.gpkg')  # insert file path to your shapefile
    #select only dead trees
    deads = gx_df[gx_df['dead_label'] == 1]
    bounds = deads.total_bounds
    grid_gdf = create_grid(bounds[0], bounds[1], bounds[2], bounds[3], threshold_sample, threshold_sample)
    joined = gpd.sjoin(deads, grid_gdf, how='left', op='within')
    picked_entries = joined.groupby('index_right').apply(lambda x: x.sample(1))
    # pick only the left index and convert to list
    gap_trees = picked_entries.index.droplevel(0)
    #calculate centroid of each polygon
    gx_df['geometry'] = gx_df['geometry'].centroid
    # from geopandas extract the coordinates of the points
    gx_df['x'] = gx_df['geometry'].apply(lambda p: p.x)
    gx_df['y'] = gx_df['geometry'].apply(lambda p: p.y)
    # Convert coordinates to a list and build spatial index
    points_list = gx_df[['x', 'y']].values.tolist()
    gx_df_spindex = KDTree(points_list)
    networks = {}
    for row in gx_df.itertuples():
        # if row index is in the list of dead trees
        if row.Index in gap_trees:
            point = Point(row.x, row.y)
            indices_within_threshold = gx_df_spindex.query_ball_point((row.x, row.y), threshold_value)
            if len(indices_within_threshold) > 0:
                network = nx.DiGraph()  # Create a new network for each DEAD tree
                for index in indices_within_threshold:
                    neighbor_point = Point(*points_list[index])
                    #calculate distance between point and dead_point
                    dist = point.distance(neighbor_point)
                    strength = edge_strength(dist, threshold_value)
                    network.add_edge(row.Index, index, tag = gx_df.sci_name[index], dist=dist, strength=strength, is_dead = gx_df.dead_label[index],
                                    crown_area = gx_df.crown_area[index], height = gx_df.height[index], east = gx_df.x[index], north = gx_df.y[index])
                networks[row.Index] = network
    #save networks to file
    for idx, (name, network) in enumerate(networks.items()):
        with open(f'networks/{site}/{name}.gpickle', 'wb') as f:
            pickle.dump(network, f, pickle.HIGHEST_PROTOCOL)


site = ["SERC", "OSBS", "HARV"]
threshold_value = 50
threshold_sample = 100
for st in site:
    create_networks(st, threshold_value, threshold_sample)

for st in site:
    get_networks_with_no_dead(st, threshold_value, threshold_sample)
