library(terra)
library(sf)
library(exactextractr)

# Pull municipal boundaries for Brazil
itcs <- sf::read_sf('data/SERC.gpkg')

# Pull gridded precipitation data
chm <- terra::rast('data/SERC.tif')

# Calculate vector of mean December precipitation amount for each municipality
itcs$height <- exact_extract(chm[[1]], itcs, 'quantile', quantiles = 0.98)

sf::write_sf(itcs, 'data/SERCh.gpkg')
