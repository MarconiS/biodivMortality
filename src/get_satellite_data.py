import openeo
import geopandas as gpd
import numpy as np
site = "HARV"
# get the area of interest from the AOP footprint
aoi = gpd.read_file('data/AOP_flightBoxes/AOP_flightboxesAllSites.shp')
aoi = aoi[aoi['siteID'] == site]

#extract extent of aoi
extent = aoi.total_bounds

# Connect to openEO back-end.
connection = openeo.connect("openeo.vito.be").authenticate_oidc()
connection.list_collection_ids()
# Load data cube from TERRASCOPE_S2_NDVI_V2 collection
year = 2017

#loop through the range of years
for year in range(2000, 2023):
    try:
        cube = connection.load_collection(
            "LANDSAT8-9_L2",
            spatial_extent={"west": extent[0], "south": extent[1], "east": extent[2], "north": extent[3]},
            temporal_extent=[f"{year}-05-01", f"{year}-06-30"],
            bands=['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B10', 'dataMask'],
        )
        cube.download(f"satellite/{year}_{site}.tiff")
    except:
        print(f"no data for {year}")
        continue

import ee
import openeo
import geopandas as gpd
import numpy as np
from ee import batch

import openeo
import geopandas as gpd
import numpy as np
site = ["SERC", "OSBS", "HARV"]
for st in site:
    # get the area of interest from the AOP footprint
    aoi = gpd.read_file('data/AOP_flightBoxes/AOP_flightboxesAllSites.shp')
    aoi = aoi[aoi['siteID'] == st]
    #extract extent of aoi
    extent = aoi.total_bounds
    ee.Initialize()
    for year in range(2003, 2023):
        # Create an image collection.
        collection = ee.ImageCollection('USDA/NAIP/DOQQ').filterDate(f'{year}-06-01', f'{year}-7-31')
        # Filter spatially using a geometry object.
        collection = collection.filterBounds(ee.Geometry.Rectangle(tuple(extent)))
        if collection.size().getInfo() > 0:
            image = collection.median()
            taskConfig = {
                'image': image,
                'description': f'summer_NAIP_{year}',
                'region': ee.Geometry.Rectangle(tuple(extent)).getInfo()["coordinates"],
                'folder' : 'NAIP_forestGeo',
                'fileFormat': 'GeoTIFF',
                'maxPixels': 1E10,
            }
            task = batch.Export.image.toDrive(**taskConfig)
            task.start()
        else:
            print(f"No data for {year}")