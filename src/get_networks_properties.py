# load network saved in pickle data file

import pickle
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt
import pandas as pd
import math
import geopandas as gpd

def calculate_hull_biodiversity(G):
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

    # Calculate Hull index
    hull_index = num_tags - (sum([(-freq/num_tags)*math.log(freq/num_tags) for tag, freq in tag_freq.items()]))
    
    return hull_index, num_tags, tag_freq



def count_dead_tags(directory):
    gap_features = []
    for filename in os.listdir(directory):
        if filename.endswith(".gpickle"):
            dead_count = 0
            pt = os.path.join(directory, filename)
            with open(pt, 'rb') as f:
                network = pickle.load(f)

            #get the number of edges with is_dead == True
            dead_count = len([edge for edge in network.edges(data=True) if edge[2].get('is_dead') == 0])

            # get the sum of the crown_area of edges for which is_dead == False
            surface_alive = sum([edge[2].get('crown_area') for edge in network.edges(data=True) if edge[2].get('is_dead') == 0])

            # get the sum of the crown_area of edges for which is_dead == True
            surface_dead = sum([edge[2].get('crown_area') for edge in network.edges(data=True) if edge[2].get('is_dead') == 1])

            # get the relative hill diversity
            hill_div, species_count, sp_frequency = calculate_hull_biodiversity(network)

            # create a pandas series with the dead count and the filename
            idrow = pd.Series([filename, dead_count, surface_dead, surface_alive, species_count, hill_div, sp_frequency], 
                              index=['focal_id', 'dead_count', 'surface_dead', 'surface_alive', 'species_count', 'hill_div', 'sp_frequency']) 
            # append idrow series       
            gap_features.append(idrow)
        
    return gap_features



site = "SERC"
directory = f"networks/{site}"
gap_features = count_dead_tags(directory)
gx_df = gpd.read_file(f'data/{site}.gpkg')  # insert file path to your shapefile


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
df_.to_file(f'data/{site}_gap_features.gpkg', driver='GPKG')
