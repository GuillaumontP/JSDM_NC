library(here)
library(stars)
library(ggplot2)
library(terra)
library(glue)
set.seed(1234)
EPSG <- 3163
proj.t <- paste0("EPSG:", EPSG) 

cell_site <- split(read_stars(here("output", "jSDM_data_final.tif")))[1,,]
load(here("output", "one_LV", "RData", "jSDM_binom_pro.RData"))
for (x in 1:st_dimensions(cell_site)$x$to) {
  for (y in 1:st_dimensions(cell_site)$y$to) {
    cell_site[[1]][x,y] <- st_dimensions(cell_site)$x$to * x + y 
  }
}
write_stars(cell_site, dsn = here("output", "cell_site.tif"),  options = c("COMPRESS=LZW","PREDICTOR=2"))

latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
value_cell_site <- rep(0, dim(coord)[1])
cell <- matrix(0, ncol = 5, nrow = dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  value_cell_site[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "cell_site.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
  cell[i, 1:2] <- c(value_cell_site[i] %/% st_dimensions(cell_site)$x$to, value_cell_site[i] %% st_dimensions(cell_site)$x$to)
}
colnames(cell) <- c("x", "y", "alpha", "W_1", "dist")
cell[,3] <- colMeans(jSDM_binom_pro$mcmc.alpha)
cell[,4] <- colMeans(jSDM_binom_pro$mcmc.latent$lv_1)

alpha_matrix <- matrix(0, nrow = st_dimensions(cell_site)$y$to,
                       ncol = st_dimensions(cell_site)$x$to)
W_1_matrix <- matrix(0, nrow = st_dimensions(cell_site)$y$to,
                     ncol = st_dimensions(cell_site)$x$to)
knn <- 5
for (x in 1:st_dimensions(cell_site)$x$to) {
  for (y in 1:st_dimensions(cell_site)$y$to) {
    cell[,5] <- sqrt((cell[,1] - x)^2 + (cell[,2] - y)^2) # distance from inventory sites
    value_alpha <- 0 
    value_W_1 <- 0 
    value_W_2 <- 0 
    for (i in unique(sort(cell[,5])[1:knn])) {
      # value times distance from site
      value_alpha <- value_alpha + sum(cell[cell[,5] == i, 3]) * i
      value_W_1 <- value_W_1 + sum(cell[cell[,5] == i, 4]) * i
    }
    # total sum of distance 
    weight <- sum(cell[,5][cell[,5] <= sort(cell[,5])[knn]])
    alpha_matrix[y, x] <- value_alpha / weight
    W_1_matrix[y, x] <- value_W_1 / weight
  }
}
for (i  in 1:length(cell[,1])) {
  alpha_matrix[cell[i, 2], cell[i, 1]] <- cell[i,3]
  W_1_matrix[cell[i, 2], cell[i, 1]] <- cell[i,4]
}
alpha_stars <- cell_site
W1_stars <- cell_site
alpha_stars[[1]] <- t(alpha_matrix)
W1_stars[[1]] <- t(W_1_matrix)
forest <- read_sf(here("output", "moist_forest.shp"))
GT <- read_sf(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))
write_stars(st_crop(st_crop(alpha_stars, GT), forest), here("output", "one_LV", "alpha_knn.tif"))
write_stars(st_crop(st_crop(W1_stars, GT), forest), here("output", "one_LV", "W1_knn.tif"))

