library(here)
library(stars)
library(rgdal)
library(rgrass7)
library(glue)
library(sf)
library(ggplot2)
library(terra)

proj.t <- "EPSG:3163"
load(here("output", "RData", "jSDM_binom_pro.RData"))

##==============
##
## Alpha interp vs alpha model inventory site
##
##==============

params_sites <- read.csv(here("data_raw", "NCpippn", "var_site.csv"), sep = ",")
coord <- read.csv(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
params_sites$latitude <- coord$latitude
params_sites$longitude <- coord$longitude
alpha_site_model <- colMeans(jSDM_binom_pro$mcmc.alpha)
alpha_site_interp <- rep(0, length(alpha_site_model))
for (i in 1:length(alpha_site_interp)) 
{
  alpha_site_interp[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "alphas.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE)) 
}
png(filename = here("output", "plot", "alpha_interp_vs_model.png"))
plot(alpha_site_interp, alpha_site_model, main = "Alpha model vs alpha interpolated by RST")
abline(a = 0, b = 1, col = "red")
dev.off()

##=============
##
## Alpha interpolated forest 
##
##=============

alpha <- read_stars(here("output", "alphas.tif"))
border_forest <- read_sf(here("output", "moist_forest.shp"))
alpha_crop <- st_crop(alpha, border_forest)
write_stars(alpha_crop, here("output", "alpha_forest.tif"))
alpha_interp <- alpha_crop[[1]][!is.na(alpha_crop[[1]])]
summary(alpha_interp) # interpolation higher than alpha model
summary(alpha_site_model)
colors <- c("Alpha interpolated" = "red", "Alpha site model" = "blue")
ggplot() + 
  geom_point(aes(y = sort(alpha_interp), x = 1:length(alpha_interp), color = "Alpha interpolated")) +
  geom_point(aes(y = sort(alpha_site_model), 
                 x = seq(1,length(alpha_interp), length.out = length(alpha_site_model)), color = "Alpha site model")) +
  labs(x = "", y = "alpha values", color = "") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = colors)
ggsave(here("output", "plot", "alpha_site_vs_alpha_interp_NC.png"))

##=============
##
## W1 interpolated vs W1 model inventory site
##
##=============

W1_site_model <- colMeans(jSDM_binom_pro$mcmc.latent$lv_1)
W1_site_interp <- rep(0, length(W1_site_model))
for (i in 1:length(alpha_site_interp)) 
{
  W1_site_interp[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W1.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE)) 
}
png(filename = here("output", "plot", "W1_interp_vs_model.png"))
plot(W1_site_interp, W1_site_model, main = "W1 model vs W1 interpolated by RST")
abline(a = 0, b = 1, col = "red")
dev.off()

##=============
##
## W1 interpolated forest 
##
##=============

W1 <- read_stars(here("output", "lv_W1.tif"))
border_forest <- read_sf(here("output", "moist_forest.shp"))
W1_crop <- st_crop(W1, border_forest)
write_stars(W1_crop, here("output", "lv_W1_forest.tif"))
W1_interp <- W1_crop[[1]][!is.na(W1_crop[[1]])]
summary(W1_interp) 
summary(W1_site_model)
colors <- c("W1 interpolated" = "red", "W1 site model" = "blue")
ggplot() + 
  geom_point(aes(y = sort(W1_interp), x = 1:length(W1_interp), color = "W1 interpolated")) +
  geom_point(aes(y = sort(W1_site_model), 
                 x = seq(1,length(W1_interp), length.out = length(W1_site_model)), color = "W1 site model")) +
  labs(x = "", y = "W1 values", color = "") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = colors)
ggsave(here("output", "plot", "W1_site_vs_W1_interp_NC.png"))

##=============
##
## W2 interpolated vs W2 model inventory site
##
##=============

W2_site_model <- colMeans(jSDM_binom_pro$mcmc.latent$lv_1)
W2_site_interp <- rep(0, length(W2_site_model))
for (i in 1:length(alpha_site_interp)) 
{
  W2_site_interp[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W2.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE)) 
}
png(filename = here("output", "plot", "W2_interp_vs_model.png"))
plot(W2_site_interp, W2_site_model, main = "W2 model vs W2 interpolated by RST")
abline(a = 0, b = 1, col = "red")
dev.off()

##=============
##
## W2 interpolated forest 
##
##=============

W2 <- read_stars(here("output", "lv_W2.tif"))
border_forest <- read_sf(here("output", "moist_forest.shp"))
W2_crop <- st_crop(W2, border_forest)
write_stars(W2_crop, here("output", "lv_W2_forest.tif"))
W2_interp <- W2_crop[[1]][!is.na(W2_crop[[1]])]
summary(W2_interp) 
summary(W2_site_model)
colors <- c("W2 interpolated" = "red", "W2 site model" = "blue")
ggplot() + 
  geom_point(aes(y = sort(W2_interp), x = 1:length(W2_interp), color = "W2 interpolated")) +
  geom_point(aes(y = sort(W2_site_model), 
                 x = seq(1,length(W2_interp), length.out = length(W2_site_model)), color = "W2 site model")) +
  labs(x = "", y = "W2 values", color = "") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = colors)
ggsave(here("output", "plot", "W2_site_vs_W2_interp_NC.png"))

##============
##
## Pixels with high diversity alpha
##
##============
theta_sum <- read_stars(here("output", "theta_forest_sum.tif"))
theta_sum <- st_crop(theta_sum, border_forest)
coord_pixels_rich <- which(theta_sum[[1]] > 600, arr.ind = TRUE)
values_pixels <- matrix(0, ncol = 3, nrow = length(coord_pixels_rich[,1]))
for (i in 1:length(coord_pixels_rich[,1])) 
{
  values_pixels[i, 1] <- alpha_crop[[1]][coord_pixels_rich[i,1], coord_pixels_rich[i,2]]
  values_pixels[i, 2] <- W1_crop[[1]][coord_pixels_rich[i,1], coord_pixels_rich[i,2]]
  values_pixels[i, 3] <- W2_crop[[1]][coord_pixels_rich[i,1], coord_pixels_rich[i,2]]
}

# line red : median of pixel with high specific richness
# line black : median of forest pixel

plot(sort(alpha_interp), col = "darkblue", main = "alpha interpolated")
abline(b = 0, a = median(values_pixels[,1]), col = "red")
abline(b = 0, a = median(alpha_interp), col = "black")

plot(sort(W1_interp), col = "darkblue", main = "W1 interpolated")
abline(b = 0, a = median(values_pixels[,2]), col = "red")
abline(b = 0, a = median(W1_interp), col = "black")

plot(sort(W2_interp), col = "darkblue", main = "W2 interpolated")
abline(b = 0, a = median(values_pixels[,3]), col = "red")
abline(b = 0, a = median(W2_interp), col = "black")

##=============
##
## plot interpolation with & without forest
##
##=============
plot(rast(here("output", "alphas.tif")))
plot(rast(here("output", "alpha_forest.tif")))
plot(rast(here("output", "lv_W1.tif")))
plot(rast(here("output", "lv_W1_forest.tif")))
plot(rast(here("output", "lv_W2.tif")))
plot(rast(here("output", "lv_W2_forest.tif")))

##============
##
## Sum of alpha & W_i 
## forest & high diversity alpha
##============

factor_sum <- alpha + W1 + W2
write_stars(factor_sum, here("output", "sum_alpha_W_i.tif"))
plot(rast(here("output", "sum_alpha_W_i.tif")))
site_sum <- sort(colMeans(jSDM_binom_pro$mcmc.alpha) + colMeans(jSDM_binom_pro$mcmc.latent$lv_1) + 
  colMeans(jSDM_binom_pro$mcmc.latent$lv_2))
factor_sum <- st_crop(factor_sum, border_forest)
sum_pixel <- sort(factor_sum[[1]][which(!is.na(factor_sum[[1]]))])
colors <- c("pixels sum" = "red", "site sum" = "blue")
ggplot() + 
  geom_point(aes(y = sum_pixel, x = 1:length(sum_pixel), color = "pixels sum")) +
  geom_point(aes(y = site_sum, 
                 x = seq(1,length(sum_pixel), length.out = length(site_sum)), color = "site sum")) +
  labs(x = "", y = "sum alpha W_i", color = "") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = colors)
ggsave(here("output", "plot", "sum_NC_vs_high_diversity_alpha.png"))

##=============
##
## Interpolation with differents tension and smooth parameters
##
##============
dir.create(here("output", "RST"))
dir.create(here("output", "RST", "theta"))
dir.create(here("output", "RST", "plot"))
dir.create(here("output", "RST", "W_2"))
dir.create(here("output", "RST", "W_2", "plot"))
load(here("output", "RData", "jSDM_binom_pro.RData"))
border <- read_stars(here("output", "border_island.tif"))
params_sites <- read.csv(here("data_raw", "NCpippn", "var_site.csv"), sep = ",")
coord <- read.csv(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
params_sites$latitude <- coord$latitude
params_sites$longitude <- coord$longitude

#====
# Init GRASS
#====

setwd(here("output"))
Sys.setenv(LD_LIBRARY_PATH = paste("/usr/lib/grass80/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
# use a georeferenced raster
system(glue('grass -c {here("output", "jSDM_data_final.tif")} grassdata/plot'))
# connect to grass database
initGRASS(gisBase = "/usr/lib/grass80", 
          gisDbase = "grassdata", home = tempdir(), 
          location = "plot", mapset = "PERMANENT",
          override = TRUE)
#====
# Import W_2 GRASS
#====

# Import sites parameters
params_sites <- read.csv(here("data_raw", "NCpippn", "var_site.csv"), sep = ",")
coord <- read.csv(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
params_sites$latitude <- coord$latitude
params_sites$longitude <- coord$longitude
longlat <- SpatialPoints(params_sites[, c("longitude", "latitude")])
proj4string(longlat) <- CRS("+proj=longlat +ellps=GRS80 +units=m")
# lat-long to UTM58S projection 
xy = spTransform(longlat, CRS("EPSG:3163"))

W2_sp <- terra::vect(x = data.frame(W2 = colMeans(jSDM_binom_pro$mcmc.latent$lv_2),
                                    x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                     crs = xy@proj4string@projargs)
rgrass7::write_VECT(W2_sp, "W2")

#====
# RST W_2
#====
proj.t <- "EPSG:3163"
all_dis <- matrix(0, ncol = length(seq(0,2,0.2)), nrow = 50)
for (smooth in seq(0,2,0.2)) {
  for (tension in 1:50) {
    # Re-sample with RST
    system(glue("v.surf.rst --overwrite --verbose -t tension={tension} input=W2 zcolumn=W2 \\
       smooth={smooth} elevation=W2_rst "))
    # Export
    system(glue('r.out.gdal --overwrite input=W2_rst \\
             output={here("output", "RST", "W_2", paste0("RST_t", tension, "_s", smooth, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
    # Representation
    W2_in <- st_crop(read_stars(here("output", "RST", "W_2", paste0("RST_t", tension, "_s", smooth, ".tif"))),
                     border)
    write_stars(W2_in, dsn = here("output", "RST", "W_2", paste0("RST_t", tension, "_s", smooth, ".tif")), 
                options = c("COMPRESS=LZW", "PREDICTOR=2"))
    
    W2_site_model <- colMeans(jSDM_binom_pro$mcmc.latent$lv_1)
    W2_site_interp <- rep(0, length(W2_site_model))
    for (i in 1:length(W2_site_interp)) 
    {
      W2_site_interp[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} \\
                          -valonly {here("output", "RST", "W_2", paste0("RST_t", tension, "_s", smooth, ".tif"))} \\
                          -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE)) 
    }
    png(here("output", "RST", "W_2", "plot", paste0("t_", tension, "_s_", smooth, ".png")))
    plot(W2_site_interp, W2_site_model, xlab = "W2 interpolated by RST", ylab = "W2 estimated by JSDM", 
         main = "Second latent axe")
    abline(a = 0, b = 1, col = "red")
    dev.off()
    dis <- 0
    for (i in 1:length(W2_site_interp)) {
      dis <- dis + sqrt(W2_site_interp[i]^2 + W2_site_model[i]^2 - 2 * W2_site_model[i] * W2_site_interp[i])
    }
    all_dis[tension, smooth * 5 + 1] <- dis / sqrt(2)
  }
}
Sys.time()
# for (i in list.files())
  
  
best_rst <- which(all_dis == min(all_dis), arr.ind = TRUE)
predfun <- function(scaled_clim_var, params_species, rst_alpha, rst_W1, rst_W2, species.range, Xbeta_aire){
  # Get lambda values for each species
  lambda_1 <- as.matrix(params_species[, "lambda_1"])
  lambda_2 <- as.matrix(params_species[, "lambda_2"])
  beta <- as.matrix(params_species[,2:13])
  
  ## Xbeta_1
  np <- st_dimensions(merge(scaled_clim_var))$attributes$to
  Xbeta_1 <- st_as_stars(matrix(beta[1,1][[1]], ncol = ncol(rst_alpha), nrow = nrow(rst_alpha)))
  
  st_crs(Xbeta_1) <- st_crs(rst_alpha)
  st_dimensions(Xbeta_1)$X1$delta <- 1000
  st_dimensions(Xbeta_1)$X2$delta <- -1000
  st_set_bbox(Xbeta_1, st_bbox(rst_alpha))
  st_dimensions(Xbeta_1)[["X1"]]$offset <- st_dimensions(rst_alpha)[['x']]$offset
  st_dimensions(Xbeta_1)[["X2"]]$offset <- st_dimensions(rst_alpha)[['y']]$offset
  for (p in 1:np) {
    Xbeta_1 <- Xbeta_1 + scaled_clim_var[[p]] * beta[1, p + 1] 
  }
  
  ## Wlambda_1
  Wlambda_1 <- rst_W1*lambda_1[1] + rst_W2*lambda_2[1]
  ## probit_theta_1
  probit_theta_1 <- Xbeta_1 + Wlambda_1 + rst_alpha + Xbeta_aire
  probit_theta <- probit_theta_1
  remove(list = c("probit_theta_1","Wlambda_1"))
  
  ## Other species
  for (j in (species.range[1] + 1):species.range[2]) {
    
    ## Xbeta_j
    Xbeta_j <- Xbeta_1
    Xbeta_j[[1]] <- rep(beta[j,1][[1]], ncell(Xbeta_j))
    for (p in 1:(np - 1)) {
      Xbeta_j <- Xbeta_j + scaled_clim_var[[p]] * beta[j, p + 1] 
    }
    
    ## Wlambda_j
    Wlambda_j <- rst_W1 * lambda_1[j] + rst_W2 * lambda_2[j] 
    
    ## probit_theta_j
    probit_theta_j <- Xbeta_j + Wlambda_j + rst_alpha
    probit_theta <- c(probit_theta, probit_theta_j)
    remove(list=c("probit_theta_j", "Xbeta_j", "Wlambda_j"))
  }
  names(probit_theta) <- params_species$species[species.range[1]:species.range[2]]
  return(probit_theta)
}

scaled_clim_var <- split(read_stars(here("output", "jSDM_data_final.tif")))
scaled_clim_var[[1]][is.na(scaled_clim_var[[1]])] <- 0
names(scaled_clim_var) <- c("ultramafic", "bio1", "bio4", "bio15", "cwd", "prec", "bio1^2", "bio4^2", "bio15^2", "cwd^2", "prec^2")
params_species <- read.csv2(here("output", "params_species.csv"), sep = ",", dec = ".")
rst_alpha <- read_stars(here("output", "alphas.tif"))
rst_W1 <- read_stars(here("output", "lv_W1.tif"))
rst_W2 <- read_stars(here("output", "RST", "W_2", paste0("RST_t", best_rst[1], "_s", (best_rst[2] - 1) / 5, ".tif")))
np <- st_dimensions(merge(scaled_clim_var))$attributes$to
nsp <- length(colnames(jSDM_binom_pro$theta_latent))
npart <- 30
first.species <- seq(1, nsp, by = floor(nsp / npart) + 1)
latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
NC_PIPPN <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep = ",")
data_site <- matrix(0, nrow = dim(coord)[1], ncol = 12)
for (j in 1:dim(coord)[1]) {
  data_site[j, 12] <- log(as.numeric(NC_PIPPN$AREA[NC_PIPPN$plot_name == latlong$plot_name[j]][1]))
}
scale_aire <- scale(data_site[,12])
data_site[,12] <- scale_aire
mu_aire <- attr(scale_aire, "scaled:center")
sd_aire <- attr(scale_aire, "scaled:scale")
Xbeta_aire <- (log(1000) - mu_aire) / sd_aire

for (n in 1:npart){
  probit_theta <- predfun(scaled_clim_var, params_species, rst_alpha, rst_W1, rst_W2,
                          species.range = c(first.species[n], min(nsp, first.species[n] + floor(nsp / npart))),
                          Xbeta_aire = Xbeta_aire)
  write_stars(st_crop(merge(probit_theta), border), options = c("COMPRESS=LZW", "PREDICTOR=2"),
              dsn = here("output", "RST", "theta", paste0("RST_probit_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  
  # Compute probabilities of presence theta
  theta <- merge(probit_theta)
  remove(probit_theta)
  for (j in 1:st_dimensions(theta)$attributes$to) {
    theta[[1]][,,j] <- pnorm(theta[[1]][,,j])
  }
  write_stars(theta, options = c("COMPRESS=LZW", "PREDICTOR=2"), 
              dsn = here("output", "RST", "theta", paste0("RST_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  write_stars(merge(st_crop(split(theta), border_forest)), options = c("COMPRESS=LZW", "PREDICTOR=2"), 
              dsn = here("output", "RST", "theta", paste0("RST_theta_forest_", str_pad(n, width = 2, pad = "0"), ".tif")))
  remove(theta)
}

##==========
##
## Maps diversity alpha
##
##==========

list_theta <- list.files(here("output", "RST", "theta"), pattern = "RST_theta_forest", full.names = TRUE)
theta_sum <- sum(rast(list_theta[1]))
for (i in list_theta[2:npart])
{
  theta_sum <- theta_sum + sum(terra::rast(i))
}
terra::writeRaster(theta_sum, here("output", "RST", "theta_forest_sum.tif"), overwrite = TRUE)
theta_sum <- read_stars(here("output", "RST", "theta_forest_sum.tif"))

ggplot() + 
  geom_stars(data = theta_sum) +
  ggtitle("Estimated current species richness") +
  scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
  coord_fixed() +
  theme_bw() +
  labs(fill = "Number of species") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "RST", "plot", "estimated_species_richness_forest.png"))  
