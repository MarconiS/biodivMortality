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

def create_networks(site, threshold_value, threshold_sample):
    # Create a new directed graph G
    gx_df = gpd.read_file(f'data/{site}.gpkg')  # insert file path to your shapefile

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
                                    crown_area = gx_df.crown_area[index], east = gx_df.x[index], north = gx_df.y[index])
                
                networks[row.Index] = network

    #save networks to file
    for idx, (name, network) in enumerate(networks.items()):
        with open(f'networks/{site}/{name}.gpickle', 'wb') as f:
            pickle.dump(network, f, pickle.HIGHEST_PROTOCOL)


site = "SERC"
threshold_value = 50
threshold_sample = 25
create_networks(site, threshold_value, threshold_sample)