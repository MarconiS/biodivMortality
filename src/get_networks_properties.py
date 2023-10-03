# load network saved in pickle data file

import pickle
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt
import pandas as pd
import math
import geopandas as gpd
import ee
from ee import batch
import rasterio
def calculate_hill_biodiversity(G, q = 2):
    # When q = 1, the Hill number is Shannon diversity. When q = 2, the Hill number is Simpson diversity.

    # Extract tags from edges
    edge_tags = [d.get('tag') for u,v,d in G.edges(data=True)]
    # Count frequency of each tag
    tag_freq = dict()
    for tag in edge_tags:
        if tag in tag_freq:
            tag_freq[tag] += 1
        else:
            tag_freq[tag] = 1
    num_tags = len(tag_freq)
    total_individuals = sum(tag_freq.values())
    # Calculate Hull index
    # D = (SUM p_i^q)^1/(1-q)
    # S is the number of species, p_i is the proportion of species i, and q is the Hill order
    if q == 1:
        # Handle Shannon diversity calculation (when q = 1)
        shannon_entropy = -sum([(freq/total_individuals) * math.log(freq/total_individuals) for freq in tag_freq.values()])
        hill_index = math.exp(shannon_entropy)
    else:
        # Calculate Hull index for other q values
        hill_index = sum([(freq/total_individuals)**q for tag, freq in tag_freq.items()])**(1/(1-q))    #
    # hill_diversity([(freq)**q for tag, freq in tag_freq.items()], q = q)
    return hill_index, num_tags, tag_freq


def count_dead_tags(directory):
    gap_features = []
    for filename in os.listdir(directory):
        if filename.endswith(".gpickle"):
            pt = os.path.join(directory, filename)
            with open(pt, 'rb') as f:
                network = pickle.load(f)

            #get the sum of the strength of edges with is_dead == True
            dead_count = sum([edge[2].get('strength') for edge in network.edges(data=True) if edge[2].get('is_dead') == 1])

            # get the sum of the crown_area * strength of edges for which is_dead == False
            surface_alive = sum([edge[2].get('crown_area') * edge[2].get('strength') for edge in network.edges(data=True) if edge[2].get('is_dead') == 0])

            # get the sum of the crown_area * strength of edges for which is_dead == True
            surface_dead = sum([edge[2].get('crown_area') * edge[2].get('strength') for edge in network.edges(data=True) if edge[2].get('is_dead') == 1])

            # get the relative hill diversity after dropping all edges with is_dead == True
            G = network.copy()
            G.remove_edges_from([(u,v) for u,v,d in G.edges(data=True) if d.get('is_dead') == 1])
            # if G is empty, set hill_div to 0
            if len(G.edges(data=True)) == 0:
                hill_div = 0
                species_count = 0
                sp_frequency = {}
            else:           
                hill_div_2, species_count, sp_frequency = calculate_hill_biodiversity(G, q = 2)
                hill_div_1, _, _ = calculate_hill_biodiversity(G, q = 1)
                hill_div_0, _, _ = calculate_hill_biodiversity(G, q = 0)


            # measure of network length
            network_length = sum([edge[2].get('strength') for edge in network.edges(data=True)])

            # get a subset of networks with top 2 dominant species. Calculate the area of dead of dominant species, dominant alive_species and  alive rare species
            # get the top 2 dominant species
            dominant_species = sorted(sp_frequency.items(), key=lambda x: x[1], reverse=True)[:2]
            dominant_species = [sp[0] for sp in dominant_species]
            #Calculate the area of is_dead == 1 && dominant_species, is_dead == 0 && dominant_species
            dominant_dead = (sum([edge[2].get('crown_area') * edge[2].get('strength') for edge in network.edges(data=True) 
                                  if edge[2].get('is_dead') == 1 and edge[2].get('tag') in dominant_species]))
            dominant_alive = (sum([edge[2].get('crown_area') * edge[2].get('strength') for edge in network.edges(data=True) 
                                   if edge[2].get('is_dead') == 0 and edge[2].get('tag') in dominant_species]))
            
            #calculate the ratio of dominant dead over all dead
            if dead_count > 0:
                dominant_dead_ratio = dominant_dead / dead_count
            else:
                dominant_dead_ratio = 0

            # to identify new trees, pick trees wth a height in the lowest 10% of the height distribution
            # get the height of the trees
            heights = [edge[2].get('height') for edge in network.edges(data=True)]
            # get the 10th percentile
            height_threshold = np.percentile(heights, 10)
            # get the number of trees with height < height_threshold, that are alive and belong to the domiant species
            new_trees = sum([1 for edge in network.edges(data=True) if edge[2].get('height') < height_threshold and edge[2].get('is_dead') == 0 and edge[2].get('tag') in dominant_species])
            tot_new_trees = sum([1 for edge in network.edges(data=True) if edge[2].get('height') < height_threshold and edge[2].get('is_dead') == 0 ])
            # calculate the ratio of new trees over all trees
            if tot_new_trees > 0:
                new_trees_ratio = new_trees / tot_new_trees
            else:
                new_trees_ratio = 0
            #Calculate the area of is_dead == 0 && rare_species
            rare_alive = sum([edge[2].get('crown_area') * edge[2].get('strength') for edge in network.edges(data=True) if edge[2].get('is_dead') == 0 and edge[2].get('tag') not in dominant_species])
            rare_dead = sum([edge[2].get('crown_area') * edge[2].get('strength') for edge in network.edges(data=True) if edge[2].get('is_dead') == 1 and edge[2].get('tag') not in dominant_species])
            # create a pandas series with the measures and the filename
            idrow = (pd.Series([filename, dead_count, surface_dead, surface_alive, species_count, hill_div_2, hill_div_1, hill_div_0, 
                                sp_frequency, network_length, dominant_dead, dominant_alive, rare_alive, 
                                rare_dead, dominant_dead_ratio, new_trees_ratio, tot_new_trees], 
                              index=['focal_id', 'dead_count', 'surface_dead', 'surface_alive', 'species_count', 
                                     'hill_simpson', 'hill_shannon', 'hill_abundance', 'sp_frequency', 'network_length', 
                                     'dominant_dead', 'dominant_alive', 
                                     'rare_alive', 'rare_dead', 'dominant_dead_ratio', 'new_trees_ratio', 'tot_new_trees'])) 
                
            # append idrow series       
            gap_features.append(idrow)

    return gap_features


site = ["SERC", "OSBS", "HARV"]
for st in site:
    directory = f"networks_no_dead/{st}"
    gap_features = count_dead_tags(directory)
    gx_df = gpd.read_file(f'data/{st}.gpkg')  # insert file path to your shapefile

    # create a dataframe from the list of series
    df = pd.DataFrame(gap_features)
    # remove .gpickle form the focal_id
    df['focal_id'] = df['focal_id'].str.replace('.gpickle', '')

    # select only the trees with index in gap_features focal_id
    gx_df_ = gx_df[gx_df.index.isin([int(focal_id) for focal_id in df['focal_id']])]

    # add focal_id = to the index of gx_df 
    gx_df_['focal_id'] = gx_df_.index
    df.reset_index(inplace=True)
    gx_df_.reset_index(inplace=True)

    # add gx_df geometry and sci_name to df based on focal_id using pd.concat
    df['focal_id'] = df['focal_id'].astype('int64')

    df_ = pd.merge(df, gx_df_[['focal_id', 'sci_name', 'geometry']], on='focal_id')

    # turn  df_ to geodataframe and save as geopackage
    df_ = gpd.GeoDataFrame(df_, geometry='geometry')

    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling
    from rasterio.mask import mask
    import geopandas as gpd
    from shapely.geometry import box
    #load nlcd data with rasterio and assign value to each point in df_ based on spatial join 
    raster =  rasterio.open(f'data/{st}_NLCD2021.tif') 

    # function to get raster value at a point
    def get_value(row, raster):
        x = row.geometry.centroid.x
        y = row.geometry.centroid.y
        row['raster_value'] = next(raster.sample([(x, y)]))[0]
        return row

    # apply the function to each point
    points = df_.apply(get_value, raster=raster, axis=1)
    points.to_file(f'data/{st}_no_dead_gap_features.gpkg', driver='GPKG')

site = ["SERC", "OSBS", "HARV"]
for st in site:
    directory = f"networks/{st}"
    gap_features = count_dead_tags(directory)
    gx_df = gpd.read_file(f'data/{st}.gpkg')  # insert file path to your shapefile

    # create a dataframe from the list of series
    df = pd.DataFrame(gap_features)
    # remove .gpickle form the focal_id
    df['focal_id'] = df['focal_id'].str.replace('.gpickle', '')

    # select only the trees with index in gap_features focal_id
    gx_df_ = gx_df[gx_df.index.isin([int(focal_id) for focal_id in df['focal_id']])]

    # add focal_id = to the index of gx_df 
    gx_df_['focal_id'] = gx_df_.index
    df.reset_index(inplace=True)
    gx_df_.reset_index(inplace=True)

    # add gx_df geometry and sci_name to df based on focal_id using pd.concat
    df['focal_id'] = df['focal_id'].astype('int64')

    df_ = pd.merge(df, gx_df_[['focal_id', 'sci_name', 'geometry']], on='focal_id')

    # turn  df_ to geodataframe and save as geopackage
    df_ = gpd.GeoDataFrame(df_, geometry='geometry')

    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling
    from rasterio.mask import mask
    import geopandas as gpd
    from shapely.geometry import box
    #load nlcd data with rasterio and assign value to each point in df_ based on spatial join 
    raster =  rasterio.open(f'data/{st}_NLCD2021.tif') 

    # function to get raster value at a point
    def get_value(row, raster):
        x = row.geometry.centroid.x
        y = row.geometry.centroid.y
        row['raster_value'] = next(raster.sample([(x, y)]))[0]
        return row

    # apply the function to each point
    points = df_.apply(get_value, raster=raster, axis=1)
    points.to_file(f'data/{st}_gap_features.gpkg', driver='GPKG')