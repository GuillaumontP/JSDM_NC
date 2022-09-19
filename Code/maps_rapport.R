library(here)
library(tidyverse)
library(OpenStreetMap)
library(ggplot2)
Extent <- readLines(here("output/extent_short_latlong.txt"))
lat1 <- -23; lat2 <- -19.5 ; lon1 <- 162.5; lon2 <- 169

coord_site <- read.csv2("data_raw/NCpippn/coord_site.csv", sep = ",")[3:4]

# other 'type' options are "osm", "maptoolkit-topo", "bing", "stamen-toner",
# "stamen-watercolor", "esri", "esri-topo", "nps", "apple-iphoto", "skobbler";
# play around with 'zoom' to see what happens; 10 seems just right to me
sa_map <- openmap(c(lat2, lon1), c(lat1, lon2),
                  type = "bing", zoom = 10, mergeTiles = TRUE)

# reproject onto WGS84
sa_map2 <- openproj(sa_map)

# use instead of 'ggplot()'
sa_map2_plt <- OpenStreetMap::autoplot.OpenStreetMap(sa_map2) + 
  geom_point(data = coord_site, aes(x =  as.numeric(coord_site[,1]), y = as.numeric(coord_site[,2])),
             colour = "red", size = 0.5)
sa_map2_plt
