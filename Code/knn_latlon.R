library(glue)
library(here)
library(stars)
library(rgdal)
library(ggplot2)
library(stringr)
library(rgdal)

forest <- read_sf( here("output", "moist_forest.shp"))
EPSG <- 3163
proj.t <- paste0("EPSG:", EPSG) 
coord <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")[, 3:4]
data <- st_crop(read_stars(here("output", "dataNC.tif"))[,,,1], forest)
coord_pixel_latlon <- coordinates(spTransform(as_Spatial(forest), CRS("+proj=longlat +datum=WGS84")))
coord_pixel_stars <- which(!is.na(data[[1]]), arr.ind = TRUE)[, 1:2]
distance <- split(data)
knn <- 5


for (site in 1:dim(coord)[1]) # distance between pixels & each site ~16min
{
  for (pixel in 1:dim(coord_pixel_latlon)[1]) 
  {
    distance[[1]][coord_pixel_stars[pixel, 1], 
                  coord_pixel_stars[pixel, 2]] <- sqrt((coord_pixel_latlon[pixel, 1] - as.numeric(coord[site, 1]))^2 + 
                                                         (coord_pixel_latlon[pixel, 2] - as.numeric(coord[site, 2]))^2)
  }
  if (site == 1)
  {
    all_distance <- distance
  }else{
    all_distance <- c(all_distance, distance)
  }
}
write_stars(merge(all_distance), dsn = here("output", "distance_pixel_sites.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))
all_distance <- read_stars(here("output", "distance_pixel_sites.tif"))

alpha <- split(data)
alpha[[1]] <- 0
alpha <- W1 <- W2 <- st_normalize(st_crop(alpha, forest))
names(alpha) <- "alpha"
names(W1) <- "W1"
names(W2) <- "W2"
load(here("output", "RData", "jSDM_binom_pro.RData"))
coord[,3] <- colMeans(jSDM_binom_pro$mcmc.alpha)
coord[,4] <- colMeans(jSDM_binom_pro$mcmc.latent$lv_1)
coord[,5] <- colMeans(jSDM_binom_pro$mcmc.latent$lv_2)
names(coord) <- c("longitude", "latitude", "alpha")

for (site in 1:dim(coord_pixel_latlon)[1]) {
  dis_pixel <- all_distance[[1]][coord_pixel_stars[site, 1], coord_pixel_stars[site, 2],]
  knn_site_nb <- which(dis_pixel >= sort(dis_pixel, decreasing = TRUE)[knn]) # site number
  knn_dis <- dis_pixel[knn_site_nb] # distance to each site
  
  alpha_value <- sum(coord[knn_site_nb, 3] * (1 - knn_dis / sum(knn_dis)) / (length(knn_site_nb) - 1))
  W1_value <- sum(coord[knn_site_nb, 4] * (1 - knn_dis / sum(knn_dis)) / (length(knn_site_nb) - 1))
  W2_value <- sum(coord[knn_site_nb, 5] * (1 - knn_dis / sum(knn_dis)) / (length(knn_site_nb) - 1))
  
  alpha[[1]][coord_pixel_stars[site, 1], coord_pixel_stars[site, 2]] <- alpha_value
  W1[[1]][coord_pixel_stars[site, 1], coord_pixel_stars[site, 2]] <- W1_value
  W2[[1]][coord_pixel_stars[site, 1], coord_pixel_stars[site, 2]] <- W2_value
    
}
alpha[[1]] <- alpha[[1]] - mean(alpha[[1]], na.rm = TRUE)
W1[[1]] <- scale(W1[[1]])
W2[[1]] <- scale(W2[[1]])
write_stars(alpha, here("output", "alpha_knn.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))
write_stars(W1, here("output", "W1_knn.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))
write_stars(W2, here("output", "W2_knn.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))


