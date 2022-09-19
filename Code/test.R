library(stars)
library(here)
library(glue)

EPSG = 3163
nodat = -9999
proj.s <- "EPSG:4326"
proj.t <- paste("EPSG:", EPSG, sep = "")
ISO_country_code = "NCL"
# border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
#                       layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
# border <- st_transform(border[1], crs = EPSG)

# env <- split(read_stars(here("output", "environNC.tif")))
forest <- read_stars(here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km%.tif"))
clim <- st_as_stars(read_stars(here("output","current_chelsaNC.tif")))
system(glue('gdal_merge.py  -o {here("output", "all_currentNC.tif")} -of GTiff -co "COMPRESS=LZW" -co "PREDICTOR=2"  -ot Int16 -a_nodata {nodat} \\
           {here("output", "current_chelsaNC.tif")} {here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km%.tif")} '))
plot(read_stars(here("output", "all_currentNC.tif")))

plot(clim[,ceiling(1706/2):ceiling(1706 / 2 + 477),(570 - 353 - 1):570])
# all_env <- st_crop(read_stars(here("output", "all_currentNC.tif")), border)
# write_stars(all_env, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
#             type = "Int16", dsn = here("output", "all_currentNC.tif"))
