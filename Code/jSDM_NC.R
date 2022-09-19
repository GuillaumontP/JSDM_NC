library(here)
library(sp)
library(jSDM)
library(stars)
library(glue)
library(rgdal)
library(readr)
library("rnaturalearth") # plotting maps
library("rnaturalearthdata")
library("rnaturalearthhires")
library(ggplot2)
library(viridis)
library(rgrass7)
library(terra)
library(stringr)
library(ade4)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)
library(terrainr)
library(ggtern)
library(matrixStats)
library(EMCluster)
library(Rtsne) # alternative to PCA
library(rgl) # plot 3d
library(cowplot) # multiple ggplot ie par(mfrow)

##================
##
## Init data
##
##================
set.seed(1234)
EPSG <- 3163
nodat <- -9999
ISO_country_code <- "NCL"
proj.t <- paste0("EPSG:", EPSG) 
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
# Ab <- read.csv2(here("data_raw", "NCpippn", "Abondance.csv"), sep = ",")
# Ab$X <- NULL
dir.create(here("output", "plot"))
dir.create(here("output", "RData"))

# Keep only Grande Terre ie main island
border <- read_sf(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))[1]
border <- st_transform(border, EPSG)
write_sf(border, dsn = here("data_raw", "Grande_Terre", "Grande_Terre.shp"))

# only compute on forest but not on mangrove (altitude > 10)
elevation <- st_crop(read_stars(here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif")), border)
elevation[[1]][elevation[[1]] < 10] <- NA
forest <- split(read_stars(here("output", "environNC.tif")))[1,,]
forest[forest < 50] <- NA
# crop forest on elevation more than 10m
forest <- st_crop(st_crop(forest, elevation), border)
forest <- st_as_sf(forest)

# crop dry forest
# https://doi.org/10.1371/journal.pone.0252063
current <- split(read_stars(here("output", "current_chelsaNC.tif")))
bio12 <- st_crop(st_crop(current[84,,], elevation), forest) # bio12
prec <- merge(st_crop(st_crop(current[37:48,,], elevation), forest))
bio12[[1]][bio12[[1]] < 1500] <- 1 # dry
bio12[[1]][bio12[[1]] >= 1500] <- 0
prec[[1]][prec[[1]] < 100] <- 1 # dry
prec[[1]][prec[[1]] >= 100] <- 0 
for (i in 2:12) {
  prec[[1]][,,1] <- prec[[1]][,,1] + prec[[1]][,,i]
}
dry_forest <- bio12 * prec[,,,1]
dry_forest[[1]][dry_forest[[1]] > 5] <- NA # remove dry forest
forest <- st_as_sf(dry_forest) # new mask for moist forest 
write_sf(forest, here("output", "moist_forest.shp"))

# Merge climate and environ in one tif file

write_stars(c(st_normalize(st_crop(read_stars(here("output", "environ_allNC.tif")), forest)), read_stars(here("output", "current_chelsaNC.tif")), 
              along = "band"), dsn = here("output", "dataNC.tif"), 
            options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)

system(glue('gdal_translate -b 15 -b 88 -b 91 -b 99 -b 102 -b 107 \\
            {here("output","dataNC.tif")} {here("output", "jSDM_data_model_original.tif")}'))

jSDM_data_model_square <- read_stars(here("output", "jSDM_data_model_original.tif"))
jSDM_data_model_original <- read_stars(here("output", "jSDM_data_model_original.tif"))
for (i in 2:6) {
  jSDM_data_model_square[[1]][,,i] <- jSDM_data_model_square[[1]][,,i] ^ 2  
}

dataF <- split(c(jSDM_data_model_original, jSDM_data_model_square[,,,2:6], along = "band"))
data <- split(read_stars(here("output", "jSDM_data_model_original.tif")))
names(dataF) <- c("ultramafic", "bio1", "bio4", "bio12", "bio15", "cwd", "bio1^2", "bio4^2", "bio12^2", "bio15^2", "cwd^2")
dataF <- merge(dataF)
dataF <- split(st_apply(dataF, "attributes", scale))
dataF <- c(data[1,,], dataF[2:11,,])
dataF <- st_crop(dataF, forest)
write_stars(merge(dataF), dsn = here("output", "jSDM_data_final.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)
unlink(here("output", "jSDM_data_model_original.tif"))
rm(jSDM_data_model_square, jSDM_data_model_original)

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
save(mu_aire, sd_aire, file = here("output", "RData", "mu_sd_aire.RData"))
for (i in 1:dim(coord)[1]) 
{
  data_site[i, 1:11] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "jSDM_data_final.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
data_site <- data.frame(data_site)
colnames(data_site) <- c("ultramafic", "bio1", "bio4", "bio12", "bio15", "cwd", "bio1^2", "bio4^2", "bio12^2", "bio15^2", "cwd^2", "log(aire)")
data_site$ultramafic <- as.numeric(data_site$ultramafic == 1)
write.csv(data_site, here("output", "data_site.csv"))

# bio1, bio4, bio12, bio15, cwd, peridotites           
var_jSDM <- data.matrix(data_site)
var_jSDM[,1] <- as.numeric(var_jSDM[,1] == 1)
write.csv(var_jSDM, here("data_raw", "NCpippn", "var_site.csv"))

nb_species <- colSums(Ab)
# names(nb_species) <- NULL
nb_min_species <- sort(nb_species, decreasing = TRUE)[301]
Ab <- Ab[,colSums(Ab) >= nb_min_species]
Ab$X<- NULL

##================
##
## jSDM Poisson Log
## about 11h10 to run
##================

# Sys.time()
# jSDM_pois_log <- jSDM_poisson_log(
#   burnin = 5000,
#   mcmc = 10000,
#   thin = 10,
#   count_data = data.matrix(Ab),
#   site_data = var_jSDM,
#   site_formula = ~. ,
#   n_latent = 2,
#   site_effect = "random",
#   V_lambda = 1,
#   V_beta = 1,
#   shape = 0.1,
#   rate = 0.1,
#   seed = 1234,
#   verbose = 1
# )
# Sys.time()
# save(jSDM_pois_log, file = here("output", "RData", "jSDM_pois_log.RData"))
# load(here("output", "RData", "jSDM_pois_log.RData"))
# 
# ## beta_j of the five first species
# top_species <- which(colSums(Ab) >= sort(colSums(Ab), decreasing = TRUE)[5])
# colSums(PA)[top_species] 
# names(top_species) <- NULL
# par(mfrow = c(3, 2))
# np <- nrow(jSDM_pois_log$model_spec$beta_start)
# for (i in top_species) {
#   for (p in 1:np) {
#     coda::traceplot(coda::as.mcmc(jSDM_pois_log$mcmc.sp[[i]][, p]))
#     coda::densplot(coda::as.mcmc(jSDM_pois_log$mcmc.sp[[i]][, p]), 
#                    main = paste(colnames(jSDM_pois_log$mcmc.sp[[i]])[p], 
#                                 ", species : ",  names(top_species[top_species == i])))
#   }
# }
# 
# ## lambda_j of the first five species
# n_latent <- jSDM_pois_log$model_spec$n_latent
# par(mfrow = c(2, 2))
# for (j in top_species) {
#   for (l in 1:n_latent) {
#     coda::traceplot(jSDM_pois_log$mcmc.sp[[j]][, np + l])
#     coda::densplot(jSDM_pois_log$mcmc.sp[[j]][, np + l], 
#                    main = paste(colnames(jSDM_pois_log$mcmc.sp[[j]])
#                                 [np + l], ", species : ",names(top_species[top_species == j])))
#   }
# }
# 
# ## Latent variables W_i for the first two sites
# 
# par(mfrow = c(2, 2))
# for (l in 1:n_latent) {
#   for (i in 1:2) {
#     coda::traceplot(jSDM_pois_log$mcmc.latent[[paste0("lv_", l)]][, i],
#                     main = paste0("Latent variable W_", l, ", site ", i))
#     coda::densplot(jSDM_pois_log$mcmc.latent[[paste0("lv_", l)]][, i],
#                    main = paste0("Latent variable W_", l, ", site ", i))
#   }
# }
# 
# ## alpha_i of the first two sites
# plot(coda::as.mcmc(jSDM_pois_log$mcmc.alpha[, 1:2]))
# 
# ## V_alpha
# par(mfrow = c(2, 2))
# coda::traceplot(jSDM_pois_log$mcmc.V_alpha)
# coda::densplot(jSDM_pois_log$mcmc.V_alpha)
# ## Deviance
# coda::traceplot(jSDM_pois_log$mcmc.Deviance)
# coda::densplot(jSDM_pois_log$mcmc.Deviance)
# 
# 
# ## probit_theta
# par (mfrow = c(2, 1))
# hist(jSDM_pois_log$log_theta_latent, main = "Predicted log theta", xlab = "predicted log theta")
# hist(jSDM_pois_log$theta_latent, main = "Predicted theta", xlab = "predicted theta", breaks = 20)

##================
##
## jSDM Binomial Probit
## about 1h to run
##================

Sys.time()
jSDM_binom_pro <- jSDM_binomial_probit(
  burnin = 5000,
  mcmc = 10000,
  thin = 10,
  presence_data = data.matrix(PA),
  site_formula = ~.,
  site_data = var_jSDM,
  n_latent = 2,
  site_effect = "random",
  V_lambda = 1,
  V_beta = c(rep(1, 12), 0.001),
  mu_beta = c(rep(0, 12), 0.25),
  shape = 0.1,
  rate = 0.1,
  seed = 1234,
  verbose = 1
)
Sys.time()
save(jSDM_binom_pro, file = here("output","RData", "jSDM_binom_pro.RData"))
load(here("output", "RData", "jSDM_binom_pro.RData"))

top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[5])
np <- nrow(jSDM_binom_pro$model_spec$beta_start)

## beta_j of the top five species

for (j in top_species[1]) {
  for (p in 1:np) {
    png(here("output", "plot", paste0("beta_jSDM_", p, ".png")))
    par(mfrow = c(1, 2))
    coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
    coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
    mtext(paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p], ", species : ", names(top_species[top_species == j])),
          side = 3, line = - 2, outer = TRUE)
    dev.off()
  }
}


## lambda_j of the top five species
n_latent <- jSDM_binom_pro$model_spec$n_latent

for (j in top_species[1]) {
  png(here("output", "plot", "lambda_jSDM.png"))
  par(mfrow = c(2, 2))
  for (l in 1:n_latent) {
    coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
    coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]), 
                   main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])
                                [np + l],", species : ", names(top_species[top_species == j])))
   
  }
  dev.off()
}

## Latent variables W_i for the first two sites

for (l in 1:n_latent) {
  png(here("output", "plot", paste0("W_", l, "_jSDM.png")))
  par(mfrow = c(2, 2))
  for (i in 1:2) {
    coda::traceplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                    main = paste0("Latent variable W_", l, ", site ", i))
    coda::densplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                   main = paste0("Latent variable W_", l, ", site ", i))
  }
  dev.off()
}

## alpha_i of the first two sites
png(here("output", "plot", "alpha_jSDM.png"))
plot(coda::as.mcmc(jSDM_binom_pro$mcmc.alpha[, 1:2]))
dev.off()

## V_alpha
png(here("output", "plot", "V_alpha_Deviance_jSDM.png"))
par(mfrow = c(2, 2))
coda::traceplot(jSDM_binom_pro$mcmc.V_alpha)
coda::densplot(jSDM_binom_pro$mcmc.V_alpha)
## Deviance
coda::traceplot(jSDM_binom_pro$mcmc.Deviance)
coda::densplot(jSDM_binom_pro$mcmc.Deviance)
dev.off()

## probit_theta
png(here("output", "plot", "theta_jSDM.png"))
par (mfrow = c(2, 1))
hist(jSDM_binom_pro$probit_theta_latent,
     main = "Predicted probit theta", xlab = "predicted probit theta")
hist(jSDM_binom_pro$theta_latent,
     main = "Predicted theta", xlab = "predicted theta")
dev.off()

##================
##
## Plotting occurences for each species
##
##================

species_to_plot <- "Calophyllum.caledonicum"
NC_PIPPN <- read.csv2(here("data_raw", "NCpippn", "NC_PIPPN_2022.csv"), sep = ",")
coord <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
coord$longitude <- as.numeric(coord$longitude)
coord$latitude <- as.numeric(coord$latitude)
coord$X <- NULL
coord[, species_to_plot] <- PA[, species_to_plot]

NC <- ne_countries(scale = 10, returnclass = "sf")
pres <- data.frame(coord[coord[, species_to_plot] == 1, 2:3])
abs <- data.frame(coord[coord[, species_to_plot] == 0, 2:3])
save(coord, species_to_plot, file = here("output", "RData", "obs_occ.RData"))
load(here("output", "RData", "obs_occ.RData"))

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = coord[,c("longitude", "latitude", species_to_plot)], 
             aes(x = longitude, y = latitude, shape = as.logical(coord[,species_to_plot]),
                 color = as.logical(coord[,species_to_plot])), size = 1) +
  ggtitle(paste0('Observed occurences of \n', species_to_plot)) +
  scale_shape_manual(labels = c("Absence", "Presence"), values = c(1, 3), name = "Species state") +
  scale_color_manual(labels = c("Absence", "Presence"), values = c("red", "black"), name = "Species state") +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "Observed_occurences.png"))

data_theta <- data.frame(coord[,c("longitude", "latitude")])
data_theta$theta <- jSDM_binom_pro$theta_latent[, species_to_plot]
save(data_theta, species_to_plot, file = here("output", "RData", "est_prob.RData"))
load(here("output", "RData", "est_prob.RData"))

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = data_theta, 
             aes(x = longitude, y = latitude, colour = theta), size = 1) +
  ggtitle(paste0("Estimated probabilities of presence for \n", species_to_plot)) +
  scale_colour_gradientn(colours = rev(rocket(10)), name = "Prediction", limits = c(0, 1)) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "estimated_probabilities.png"))

##================
##
##Plotting Species Richness
##
##================

# sites richness 
species_richness <- data.frame(richness = rowSums(PA))
max_richness_reel <- max(species_richness)
species_richness$latitude <- coord$latitude
species_richness$longitude <- coord$longitude
species_richness$estimate <- rowSums(jSDM_binom_pro$theta_latent)

save(species_richness, file = here("output", "RData", "species_richness.RData"))
load(here("output", "RData", "species_richness.RData"))

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = species_richness, 
             aes(x = longitude, y = latitude, colour = richness), size = 1) +
  ggtitle("Observed species richness") +
  scale_colour_gradientn(colours = rev(rocket(10)), name = "Number of species", 
                         limits = c(min(species_richness$richness), max(species_richness$richness))) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "species_richness_observed.png"))

# Estimated richness

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = species_richness, 
             aes(x = longitude, y = latitude, colour = estimate), size = 1) +
  ggtitle("Estimated species richness") +
  scale_colour_gradientn(colours = rev(rocket(10)), name = "Number of species", 
                        limits = c(min(species_richness$estimate), max(species_richness$estimate))) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "species_richness_estimated.png"))

png(here("output", "plot", "species_richness_observed_vs_estimated.png"))
par(mfrow = c(1,1))
plot(rowSums(PA), rowSums(jSDM_binom_pro$theta_latent), main = "Species richness for area \n with less than 50 differents species",
     xlab = "observed", ylab = "estimated", 
     cex.main = 1.4, cex.lab = 1.3, xlim = c(0, 50), ylim = c(0, 50))
abline(a = 0, b = 1, col = 'red')
dev.off()

##================
##
## Spatial Interpolation for alpha and latents variables
##
##================

## Initialize GRASS
setwd(here("output"))
Sys.setenv(LD_LIBRARY_PATH = paste("/usr/lib/grass80/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
# use a georeferenced raster
system(glue('grass -c {here("output", "jSDM_data_final.tif")} grassdata/plot'))
# connect to grass database
initGRASS(gisBase = "/usr/lib/grass80", 
          gisDbase = "grassdata", home = tempdir(), 
          location = "plot", mapset = "PERMANENT",
          override = TRUE)


# Import sites parameters
params_sites <- read.csv(here("data_raw", "NCpippn", "var_site.csv"), sep = ",")
coord <- read.csv(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
params_sites$latitude <- coord$latitude
params_sites$longitude <- coord$longitude
longlat <- SpatialPoints(params_sites[, c("longitude", "latitude")])
proj4string(longlat) <- CRS("+proj=longlat +ellps=GRS80 +units=m")
# lat-long to UTM58S projection 
xy = spTransform(longlat, CRS("EPSG:3163"))

#====
# alpha 
#====

alpha_sp <- terra::vect(x = data.frame(alpha = colMeans(jSDM_binom_pro$mcmc.alpha),
                                       x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                        crs = xy@proj4string@projargs)
rgrass7::write_VECT(alpha_sp, "alpha")

# Re-sample with RST
# for punctual data use function v.surf.rst
system("v.surf.rst --overwrite --verbose -t tension=3 input=alpha zcolumn=alpha \\
       smooth=0.0 elevation=alpha_rst")

# Export
system(glue('r.out.gdal --overwrite input=alpha_rst \\
             output={here("output", "alphas.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))

# Representation
forest <- read_sf(here("output", "moist_forest.shp"))
alpha_in <- st_crop(read_stars(here("output", "alphas.tif")), forest)
# mean to 0
alpha_in[[1]] <- alpha_in[[1]] - mean(alpha_in[[1]], na.rm = TRUE)

write_stars(alpha_in, dsn = here("output", "alphas.tif"), options = c("COMPRESS=LZW", "PREDICTOR=2"))
plot(alpha_in, main = "Site effect alpha interpolated by RST")
alpha_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  alpha_xy[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "alphas.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
save(jSDM_binom_pro, alpha_xy, file = here("output", "RData", "alpha.RData"))
load(here("output", "RData", "alpha.RData"))

png(here("output", "plot", "alpha_estimated_vs_interpolated_by_RST"))
plot(alpha_xy, colMeans(jSDM_binom_pro$mcmc.alpha), xlab = "alpha interpolated by RST", ylab = "alpha estimated by JSDM",
     main = "Random site effect")
abline(a = 0, b = 1, col = 'red')
dev.off()

#====
# W1 
#====

W1_sp <- terra::vect(x = data.frame(W1 = colMeans(jSDM_binom_pro$mcmc.latent$lv_1),
                                    x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                     crs = xy@proj4string@projargs)
rgrass7::write_VECT(W1_sp, "W1")

# Re-sample with RST
system("v.surf.rst --overwrite --verbose -t tension=3 input=W1 zcolumn=W1 \\
       smooth=0 elevation=W1_rst ")
# Export
system(glue('r.out.gdal --overwrite input=W1_rst \\
             output={here("output", "lv_W1.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
# Representation
W1_in <- st_crop(read_stars(here("output", "lv_W1.tif")), forest)
# scale to (0,1)
W1_in[[1]] <- (W1_in[[1]] - mean(W1_in[[1]], na.rm = TRUE)) / sd(W1_in[[1]], na.rm = TRUE)

write_stars(W1_in, dsn = here("output", "lv_W1.tif"), options = c("COMPRESS=LZW", "PREDICTOR=2"))
plot(W1_in, main = "Latent variable W1 interpolated by RST")
W1_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  W1_xy[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W1.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
save(jSDM_binom_pro, W1_xy, file = here("output", "RData", "W1.RData"))
load(here("output", "RData", "W1.RData"))

png(here("output", "plot", "W1_estimated_vs_interpolated_by_RST"))
plot(W1_xy, colMeans(jSDM_binom_pro$mcmc.latent$lv_1), xlab = "W1 interpolated by RST", ylab = "W1 estimated by JSDM", 
     main = "First latent axe")
abline(a = 0, b = 1, col = 'red')
dev.off()

#====
# W2 
#====

W2_sp <- terra::vect(x = data.frame(W2 = colMeans(jSDM_binom_pro$mcmc.latent$lv_2),
                                    x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                     crs = xy@proj4string@projargs)
rgrass7::write_VECT(W2_sp, "W2")

# Re-sample with RST
system("v.surf.rst --overwrite --verbose -t tension=3 input=W2 zcolumn=W2 \\
       smooth=0.0 elevation=W2_rst ")
# Export
system(glue('r.out.gdal --overwrite input=W2_rst \\
             output={here("output", "lv_W2.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
# Representation
W2_in <- st_crop(read_stars(here("output", "lv_W2.tif")), forest)
# scale to (0,1)
W2_in[[1]] <- (W2_in[[1]] - mean(W2_in[[1]], na.rm = TRUE)) / sd(W2_in[[1]], na.rm = TRUE)

write_stars(W2_in, dsn = here("output", "lv_W2.tif"), options = c("COMPRESS=LZW", "PREDICTOR=2"))
plot(W2_in, main = "Latent variable W2 interpolated by RST")
W2_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  W2_xy[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W2.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
save(jSDM_binom_pro, W2_xy, file = here("output", "RData", "W2.RData"))
load(here("output", "RData", "W2.RData"))

png(here("output", "plot", "W2_estimated_vs_interpolated_by_RST"))
plot(W2_xy, colMeans(jSDM_binom_pro$mcmc.latent$lv_2), xlab = "W2 interpolated by RST", ylab = "W2 estimated by JSDM", 
     main = "Second latent axe")
abline(a = 0, b = 1, col = 'red')
dev.off()

##================
##
## Compute probabilities for each species
##
##================

rst_alpha <- read_stars(here("output", "alphas.tif"))
rst_W1 <- read_stars(here("output", "lv_W1.tif"))
rst_W2 <- read_stars(here("output", "lv_W2.tif"))
knn_alpha <- read_stars(here("output", "alpha_knn.tif"))
knn_W1 <- read_stars(here("output", "W1_knn.tif"))
knn_W2 <- read_stars(here("output", "W2_knn.tif"))

load(here("output", "RData", "jSDM_binom_pro.RData"))
data_site <- read.csv2(here("output", "data_site.csv"), sep = ",")[, -1]

# Species parameters 
### fixed species effect lambdas and betas 

# Current climatic variables 
X <- data.frame(intercep = rep(1, nrow = data_site), data_site)
scaled_clim_var <- split(read_stars(here("output", "jSDM_data_final.tif")))
scaled_clim_var[[1]][is.na(scaled_clim_var[[1]])] <- 0
scaled_clim_var <- st_normalize(st_crop(scaled_clim_var, forest))
names(scaled_clim_var) <- c("ultramafic", "bio1", "bio4", "bio15", "cwd", "prec", "bio1^2", "bio4^2", "bio15^2", "cwd^2", "prec^2")

np <- st_dimensions(merge(scaled_clim_var))$attributes$to
nsp <- length(colnames(jSDM_binom_pro$theta_latent))
npart <- 30
species <- colnames(jSDM_binom_pro$model_spec$presence_data)
n_latent <- length(jSDM_binom_pro$mcmc.latent)
lambdas <- matrix(0, nsp, n_latent)
betas <- matrix(0, nsp, np + 2)
for (j in 1:nsp){
  for (l in 1:n_latent){
    lambdas[j, l] <- mean(jSDM_binom_pro$mcmc.sp[[j]][, np + 2 + l])
  }
  for (p in 1:(np + 2)){
    betas[j, p] <- mean(jSDM_binom_pro$mcmc.sp[[j]][, p])
  }
}
colnames(betas) <- colnames(jSDM_binom_pro$mcmc.sp[[1]])[1:(np + 2)]
params_species <- data.frame(species = species,
                             betas, lambda_1 = lambdas[, 1], lambda_2 = lambdas[, 2])
write.csv(params_species, file = here("output", "params_species.csv"), row.names = F)
params_species <- read.csv2(here("output", "params_species.csv"), sep = ",", dec = ".")

# Function to compute probit_theta at New Caledonia scale 
predfun <- function(scaled_clim_var, params_species, knn_alpha, knn_W1, knn_W2, species.range, Xbeta_aire){
  # Get lambda values for each species
  lambda_1 <- as.matrix(params_species[, "lambda_1"])
  lambda_2 <- as.matrix(params_species[, "lambda_2"])
  beta <- as.matrix(params_species[,2:13])
  
  ## Xbeta_1
  np <- st_dimensions(merge(scaled_clim_var))$attributes$to
  Xbeta_1 <- st_as_stars(matrix(beta[1,1][[1]], ncol = ncol(knn_alpha), nrow = nrow(knn_alpha)))

  st_crs(Xbeta_1) <- st_crs(knn_alpha)
  st_dimensions(Xbeta_1)$X1$delta <- 1000
  st_dimensions(Xbeta_1)$X2$delta <- -1000
  st_set_bbox(Xbeta_1, st_bbox(knn_alpha))
  st_dimensions(Xbeta_1)[["X1"]]$offset <- st_dimensions(knn_alpha)[['x']]$offset
  st_dimensions(Xbeta_1)[["X2"]]$offset <- st_dimensions(knn_alpha)[['y']]$offset
  for (p in 1:np) {
    Xbeta_1 <- Xbeta_1 + scaled_clim_var[[p]] * beta[1, p + 1] 
  }
  
  ## Wlambda_1
  Wlambda_1 <- knn_W1*lambda_1[1] + knn_W2*lambda_2[1]
  ## probit_theta_1
  probit_theta_1 <- Xbeta_1 + Wlambda_1 + knn_alpha + Xbeta_aire
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
    Wlambda_j <- knn_W1 * lambda_1[j] + knn_W2 * lambda_2[j] 
    
    ## probit_theta_j
    probit_theta_j <- Xbeta_j + Wlambda_j + knn_alpha
    probit_theta <- c(probit_theta, probit_theta_j)
    remove(list=c("probit_theta_j", "Xbeta_j", "Wlambda_j"))
  }
  names(probit_theta) <- params_species$species[species.range[1]:species.range[2]]
  return(probit_theta)
}

# Compute theta in n parts because it's too large

dir.create(here("output", "theta"))
first.species <- seq(1, nsp, by = floor(nsp / npart) + 1)
load(here("output", "RData", "mu_sd_aire.RData"))
Xbeta_aire <- (log(1000) - mu_aire) / sd_aire
for (n in 1:npart){
  # change RST by knn
  probit_theta <- predfun(scaled_clim_var, params_species, knn_alpha, knn_W1, knn_W2,
                          species.range = c(first.species[n], min(nsp, first.species[n] + floor(nsp / npart))),
                          Xbeta_aire = Xbeta_aire)
  write_stars(st_crop(merge(probit_theta), forest), options = c("COMPRESS=LZW", "PREDICTOR=2"),
              dsn = here("output", "theta", paste0("KNN_probit_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  
  # Compute probabilities of presence theta
  theta <- merge(probit_theta)
  remove(probit_theta)
  for (j in 1:st_dimensions(theta)$attributes$to) {
    theta[[1]][,,j] <- pnorm(theta[[1]][,,j])
  }
  write_stars(theta, options = c("COMPRESS=LZW", "PREDICTOR=2"), 
              dsn = here("output", "theta", paste0("KNN_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  write_stars(merge(st_crop(split(theta), forest)), options = c("COMPRESS=LZW", "PREDICTOR=2"), 
              dsn = here("output", "theta", paste0("KNN_theta_forest_", str_pad(n, width = 2, pad = "0"), ".tif")))
  remove(theta)
}

##================
##
## Predictive maps of presence probabilities
##
##================

species_to_plot <- "Acropogon.macrocarpus"
theta <- read_stars(here("output", "theta", "KNN_theta_forest_01.tif"))
latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
NC_PIPPN <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep = ",")


# Observed presence absence
id_pres <- which(PA[, species_to_plot] == 1)
obs_pres <- data.frame(lat = coord[id_pres, 2], long = coord[id_pres, 1])
obs_pres <- st_as_sf(obs_pres, coords = c("long", "lat"), crs = 3163)

# Plotting for one specie
save(theta, species_to_plot, file = here("output", "RData", "interpolated_current_specie.RData"))
load(here("output", "RData", "interpolated_current_specie.RData"))
GT <- readOGR(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))

ggplot() + 
  geom_polygon(data = GT, aes(x = long, y = lat), colour = "black", fill = "grey") +
  geom_stars(data = theta) +
  ggtitle(paste('Interpolated current probabilities of presence \n for', species_to_plot)) +
  # geom_sf(data = obs_pres, shape = 3, color = "white", size = 2) 
  scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
  coord_fixed() +
  theme_bw() +
  labs(fill = "Number of species") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "interpolated_presence_for_acropogon_NC.png"))

##================
##
## Predictive maps of species richness
##
##================

list_theta <- list.files(here("output", "theta"), pattern = "KNN_theta_forest", full.names = TRUE)
theta_sum <- sum(rast(list_theta[1]))
npart <- length(list_theta)
for (i in list_theta[2:npart])
{
 theta_sum <- theta_sum + sum(terra::rast(i))
}
terra::writeRaster(theta_sum, here("output", "theta_forest_sum.tif"), overwrite = TRUE)
theta_sum <- read_stars(here("output", "theta_forest_sum.tif"))
GT <- readOGR(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))

ggplot() + 
  geom_polygon(data = GT, aes(x = long, y = lat), colour = "black", fill = "grey") +
  geom_stars(data = theta_sum) +
  ggtitle("Estimated current species richness") +
  scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
  coord_fixed() +
  theme_bw() +
  labs(fill = "Number of species") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "estimated_species_richness_forest.png"))  

##================
##
## Init PCA and tSNE
## 
##================

# get all cells for PCA
theta_stars <- read_stars(here("output", "theta", "KNN_theta_forest_01.tif"))[[1]]
ultramafic <- st_crop(read_stars(here("output", "environ_allNC.tif"))[,,,15], forest)[[1]]
theta <- terra::rast(here("output", "theta", "KNN_theta_forest_01.tif"))
nb_cell <- dim(values(theta, na.rm = TRUE))
theta_matrix <- matrix(0, nrow = nb_cell[1], ncol = nb_cell[2])
cell <- matrix(0, ncol = 2, nrow = dim(theta_matrix)[1])
UM <- rep(-1, dim(theta_matrix)[1])
increment <- 1
for (x in 1:dim(theta_stars)[1])
{
  for (y in 1:dim(theta_stars)[2]) 
  {
    if ( sum(is.na(theta_stars[x,y,])) == 0)
    {
      cell[increment,] <- c(x, y)
      theta_matrix[increment,] <- theta_stars[x,y,]
      UM[increment] <- ultramafic[x, y, 1]
      increment <- increment + 1
    }
  }
}
save(cell, file = here("output","RData", "cell.RData"))

# get position of each cell
npart <- 30
for (file in 2:npart) {
  increment <- 1
  theta_stars <- read_stars(here("output", "theta", paste0("KNN_theta_forest_",str_pad(file, 2, pad = "0") , ".tif")))[[1]]
  theta_matrix_k <- matrix(0, nrow = nb_cell[1], ncol = dim(theta_stars)[3])
  for (x in 1:dim(theta_stars)[1])
  {
    for (y in 1:dim(theta_stars)[2]) 
    {
      if ( sum(is.na(theta_stars[x,y,])) == 0)
      {
        theta_matrix_k[increment,] <- theta_stars[x,y,]
        increment <- increment + 1
      }
    }
  }
  theta_matrix <- cbind(theta_matrix, theta_matrix_k)
}
theta_df <- data.frame(theta_matrix)
rm(theta_matrix_k)
save(theta_df, file = here("output", "RData", "theta_df.RData"))

# get environ variables for each cell 
forest <- read_sf(here("output", "moist_forest.shp"))
env <- st_crop(read_stars(here("output", "jSDM_data_final.tif")), forest)[[1]]
env_cell <- matrix(0, ncol = dim(env)[3], nb_cell[1])
increment <- 1
for (x in 1:dim(env)[1])
{
  for (y in 1:dim(env)[2]) 
  {
    if ( sum(is.na(env[x,y,])) == 0)
    {
      env_cell[increment,] <- env[x,y,]
      increment <- increment + 1
    }
  }
}
colnames(env_cell) <- c("ultramafic", "T_mean", "T_seas", "prec", "P_seas", "cwd", "T_mean^2", "T_seas^2", "prec^2", "P_seas^2", "cwd^2")
env_cell <- data.frame(env_cell)

##========
## t-SNE  
##========

# t SNE 
load(here("output", "RData", "theta_df.RData"))
tSNE_fit <- Rtsne(scale(theta_df), dims = 3, max_iter = 5000, num_threads = 0,
                  verbose = TRUE)
tSNE_df <- as.data.frame(tSNE_fit$Y)
ggplot(data = tSNE_df, aes(x = V1, y = V2)) +
  geom_point() +
  theme(legend.position = "bottom")

# CAH with tSNE coord
dist_tSNE <- dist(tSNE_df)
n = nrow(theta_df) # add line where cut dendo
nb_class <- 5
# nb_class_2 <- 9
MidPoint = (hclust(dist_tSNE)$height[n - nb_class] + hclust(dist_tSNE)$height[n - nb_class + 1]) / 2
# MidPoint2 = (hclust(dist_tSNE)$height[n - nb_class_2] + hclust(dist_tSNE)$height[n - nb_class_2 + 1]) / 2

png(here("output", "plot", "dendo_tSNE.png"))
plot(hclust(dist_tSNE), labels = FALSE, main = "", xlab = "", ylab = "", sub = "", ylim = "none") 
abline(h = MidPoint, lty = 2, col = "red")
abline(h = MidPoint2, lty = 2, col = "red")
dev.off()

##========
## PCA
##========

# PCA with explanatory variable in sup
pca_with_env <- PCA(cbind(theta_df, env_cell[, 1:6]), quanti.sup = dim(theta_df)[2] + 1:6, graph = FALSE)
save(pca_with_env, file = here("output", "RData", "pca_with_env.RData"))
load(here("output", "RData", "pca_with_env.RData"))

# keep axis while explained variance is more than 1% (7 in this case)
png(here("output", "plot", "PCA_explained_variances.png"))
fviz_eig(pca_with_env)
dev.off()
print(paste0("Explained variances by first seven axis : ", round(sum(fviz_eig(pca_with_env)$data[1:7,2])), "%"))

# CAH with PCA coord
dist_site <- dist(get_pca_ind(pca_with_env)$coord, method = "euclidian")
png(here("output", "plot", "dendo_ACP.png"))
plot(hclust(dist_site), labels = FALSE, main = "PCA dendrogram") 
dev.off()

##===========
##
## K-means with t-SNE coords
## plot KM groups on t-SNE
##===========

KM_class <- kmeans(tSNE_df, centers = nb_class, iter.max = 20, nstart = 1000)
tSNE_site_group <- KM_class$cluster
table(tSNE_site_group)

save(tSNE_df, tSNE_site_group, file = here("output", "RData", "tSNE.RData"))
load(here("output", "RData", "tSNE.RData"))

ggplot(data = cbind(tSNE_df, tSNE_site_group),
       aes(x = V1, y = V2, color = as.factor(tSNE_site_group))) +
  scale_color_manual(values = viridis(nb_class)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(here("output", "plot", "tSNE_axis_1_2_group.png"))

ggplot(data = cbind(tSNE_df, tSNE_site_group),
       aes(x = V2, y = V3, color = as.factor(tSNE_site_group))) +
  scale_color_manual(values = viridis(nb_class)) +
  geom_point() + 
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) 
ggsave(here("output", "plot", "tSNE_axis_2_3_group.png"))

##=====
##
## Plot PCA with KM groups from tSNE
##
##=====

save(pca_with_env, tSNE_site_group, file = here("output", "RData", "pca_tSNE.RData"))
load(here("output", "RData", "pca_tSNE.RData"))
# Axis 1 & 2 no var
fviz_pca_ind(pca_with_env, label = "none", fill = "orange", pointshape = 21, col.ind = "NA") + 
  theme(legend.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))
ggsave(here("output", "plot", "PCA_simple_axis_1_2.png"))

# Axis 1 and 2
fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), repel = TRUE, col.quanti.sup = "black",
                invisible = "var", fill.ind = as.factor(tSNE_site_group), pointshape = 21, labelsize = 7, axes = c(1, 2),
                col.ind = "NA", alpha = 0.8) +
  ggpubr::fill_palette(viridis(nb_class)) +
  ggtitle("PCA on species with environ variable in sup") +
  guides(fill = guide_legend(title = "Groups by t-SNE algo")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
ggsave(here("output", "plot", "PCA_with_environ_sup_axis_1_2_tSNE.png"))

# Axis 2 and 3
fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), repel = TRUE, col.quanti.sup = "black",
                invisible = "var", fill.ind = as.factor(tSNE_site_group), pointshape = 21, labelsize = 5, axes = c(2, 3),
                col.ind = "NA", alpha = 0.8) +
  ggpubr::fill_palette(viridis(nb_class)) +
  ggtitle("PCA on species with environ variable in sup") +
  guides(fill = guide_legend(title = "Groups by tSNE algo")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "PCA_with_environ_sup_axis_2_3_tSNE.png"))

##==============
##
## Plot PCA 2 & 3 with ultramafics soil
## Plot t-SNE with ultramafics soil
##===============

col <- UM
col[col == 1] <- "Ultramafic"
col[col == 0] <- "No Ultramafic"
save(col, file = here("output", "RData", "UM_col.RData"))
load(here("output", "RData", "UM_col.RData"))

# PCA & t-SNE Ultramafic distinction
load(here("output", "RData", "pca_with_env.RData"))
fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), repel = TRUE, col.quanti.sup = "black",
                invisible = "var", fill.ind = col, pointshape = 21, labelsize = 7, axes = c(2, 3), 
                 col.ind = "NA", alpha = 0.6, 
                legend.title = "Type of soil", title = "PCA : soil area") +
  ggpubr::fill_palette(viridis(2)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "PCA_Axis_2-3_UM_NUM.png"))

col3d_UM <- UM
col3d_UM[col3d_UM == 1] <- viridis(2)[1]
col3d_UM[col3d_UM == 0] <- viridis(2)[2]
par3d(windowRect = c(20, 30, 800, 800))
plot3d(tSNE_df$V1, tSNE_df$V2, tSNE_df$V3, col = col3d_UM, xlab = "x", ylab = "y", zlab = "altitude")
legend3d("bottomright", legend = c("Ultramafic", "No Ultramafic"), pch = 16, col = viridis(2),
         cex = 2, inset = c(0.02))

save(col3d_UM, tSNE_df, file = here("output", "RData", "plot3d_UM.RData"))
##=====
## Change the coordinate scale for [10.255].
## for plot color for each pixel
##=====
save(pca_with_env, theta_df, tSNE_df, file = here("output", "RData", "RGB.RData"))
load(here("output", "RData", "RGB.RData"))
load(here("output", "RData", "cell.RData"))

# init stars object with CRS, extent, NA,...
rastRGB <- read_stars(here("output", "theta", "KNN_theta_forest_01.tif"))[,,,1:3]
# PCA_site_coord <- get_pca_ind(pca_with_env)$coord[,1:3]
tSNE_site_coord <- theta_df
# RGB_value_PCA <- PCA_site_coord
RGB_value_tSNE <- tSNE_df

for (l in 1:3) { # set to [0,255]
  # RGB_value_PCA[, l] <- PCA_site_coord[, l] - min(PCA_site_coord[, l])
  # RGB_value_PCA[, l] <- RGB_value_PCA[, l] / max(RGB_value_PCA[, l]) * 255
  RGB_value_tSNE[, l] <- tSNE_site_coord[, l] - min(tSNE_site_coord[, l])
  RGB_value_tSNE[, l] <- RGB_value_tSNE[, l] / max(RGB_value_tSNE[, l]) * 255
}

# put RGB color in raster and one col by group 
# rastRGB_PCA <- rastRGB
rastRGB_tSNE <- rastRGB
for (i in 1:dim(get_pca_ind(pca_with_env)$coord[,1:3])[1]){
  # rastRGB_PCA[[1]][cell[i, 1], cell[i, 2],] <- unlist(RGB_value_PCA[i,])
  rastRGB_tSNE[[1]][cell[i, 1], cell[i, 2],] <- unlist(RGB_value_tSNE[i,])
}

# Coloration RGB
# write_stars(rastRGB_PCA, dsn = here("output", "RGB_forest_PCA.tif"), 
#             options = c("COMPRESS=LZW", "PREDICTOR=2"))
write_stars(rastRGB_tSNE, dsn = here("output", "RGB_forest_tSNE.tif"), 
            options = c("COMPRESS=LZW", "PREDICTOR=2"))

##================
##
## Plot R_G_B and RGB with PCA coord & t-SNE coord
## Plot KM groups with t-SNE coord
##================

# prediction on forest area

# RGB_forest_PCA <- terra::rast(here("output", "RGB_forest_PCA.tif"))
RGB_forest_tSNE <- terra::rast(here("output", "RGB_forest_tSNE.tif"))
colorize(RGB_forest_tSNE, value = 1:3, to = "col", stretch = "hist", filename = here("output", "RGB_forest_tSNE_one.tif"), overwrite = TRUE)
# RGB_forest_PCA <- read_stars(here("output", "RGB_forest_PCA.tif"))
RGB_forest_tSNE <- read_stars(here("output", "RGB_forest_tSNE.tif"))
R_tSNE <- split(RGB_forest_tSNE)[1,,]
names(R_tSNE) <- "Axis 1 tSNE"
G_tSNE <- split(RGB_forest_tSNE)[2,,]
names(G_tSNE) <- "Axis 2 tSNE"
B_tSNE <- split(RGB_forest_tSNE)[3,,]
names(B_tSNE) <- "Axis 3 tSNE"
GT <- readOGR(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))

R_tSNE_plot <- ggplot() + 
  geom_polygon(data = GT, aes(x = long, y = lat), colour = "black", fill = "grey") +
  geom_stars(data = R_tSNE) +
  scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"), na.value = "transparent") +
  coord_fixed() +
  theme_bw()
G_tSNE_plot <- ggplot() + 
  geom_polygon(data = GT, aes(x = long, y = lat), colour = "black", fill = "grey") +
  geom_stars(data = G_tSNE) +
  scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Greens"), na.value = "transparent") +
  coord_fixed() +
  theme_bw()
B_tSNE_plot <- ggplot() + 
  geom_polygon(data = GT, aes(x = long, y = lat), colour = "black", fill = "grey") +
  geom_stars(data = B_tSNE) +
  scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Blues"), na.value = "transparent") +
  coord_fixed() +
  theme_bw()
R_G_B_plot <- plot_grid(R_tSNE_plot, G_tSNE_plot, B_tSNE_plot, ncol = 1, nrow = 3)
R_G_B_plot
save_plot(filename = here("output", "plot", "R_G_B_forest_tSNE.png"), plot = R_G_B_plot, nrow = 3)

# get mean color of pixels for each group
load(here("output","RData", "cell.RData"))
load(here("output", "RData", "pca_tSNE.RData"))
RGB_forest_tSNE <- read_stars(here("output", "RGB_forest_tSNE_one.tif"))
col_palette <- attr(RGB_forest_tSNE$RGB_forest_tSNE_one.tif, "color")
pixel_col_hist <- matrix(0, ncol = 4, nrow = dim(cell)[1])
for (n in 1:dim(cell)[1]) {
    pixel_col_hist[n, 1:3] <- t(col2rgb(col_palette[RGB_forest_tSNE[[1]][cell[n, 1], cell[n, 2]]]))
}
pixel_col_hist[,4] <- tSNE_site_group
col_mean_group <- matrix(0, ncol = 3, nrow = nb_class) 
for (group in 1:nb_class) {
  col_mean_group[group, ] <- colMeans(pixel_col_hist[pixel_col_hist[, 4] == group, 1:3]) # mean col for each group
}
RGB_forest_group <- read_stars(here("output", "RGB_forest_tSNE.tif"))
for (i in 1:dim(cell)[1]) {
  RGB_forest_group[[1]][cell[i,1], cell[i,2],] <- col_mean_group[pixel_col_hist[i, 4], ]
}
write_stars(RGB_forest_group, here("output", "RGB_forest_group.tif"))
# plotRGB(terra::rast(here("output", "RGB_forest_group.tif")))

RGB_forest_tSNE <- terra::rast(here("output", "RGB_forest_tSNE.tif"))
png(here("output", "plot", "RGB_forest_tSNE.png"))
par(mfrow = c(1, 1))
plotRGB(RGB_forest_tSNE, stretch = "hist")
title(main = "Predictives species community \n on forest area with tSNE",
      cex.main = 1.5, line = -3)
dev.off()

RGB_forest_group <- terra::rast(here("output", "RGB_forest_group.tif"))
color_group_hex <- rep(0, nb_class)
for (col  in 1:nb_class) {
  color_group_hex[col] <- rgb2hex(r = round(col_mean_group[col, 1]), g = round(col_mean_group[col, 2]), b = round(col_mean_group[col, 3]))
} 
ggplot() +
  geom_polygon(data = GT, aes(x = long, y = lat), colour = "black", fill = "white") +
  geom_tile() +
  geom_spatial_rgb(data = RGB_forest_group, 
                   aes(x = x, y = y, r = red, g = green, b = blue, asp = 1)) +
  coord_fixed() +
  geom_label(aes(x = 280000, y = 200000, label = "Group 1", size = 5), fill = color_group_hex[1]) + # 
  geom_label(aes(x = 280000, y = 215000, label = "Group 2", size = 5), fill = color_group_hex[2]) +
  geom_label(aes(x = 280000, y = 230000, label = "Group 3", size = 5), fill = color_group_hex[3]) +
  geom_label(aes(x = 280000, y = 245000, label = "Group 4", size = 5), fill = color_group_hex[4]) +
  geom_label(aes(x = 280000, y = 260000, label = "Group 5", size = 5), fill = color_group_hex[5]) +
  # geom_label(aes(x = 280000, y = 275000, label = "Group 6", size = 5), fill = color_group_hex[6]) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "longitude", y = "latitude")
ggsave(here("output", "plot", "RGB_forest_group_tSNE.png"))

# 20 species more likely on site of each group
load(here("output", "RData", "RGB.RData"))
max_species_group <- matrix(0, ncol = nb_class * 3, nrow = 20) 
min_species_group <- matrix(0, ncol = nb_class * 3, nrow = 20) 
data_clear <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep =",")
params_species <- read.csv2(here("output", "params_species.csv"), sep = ",", dec = ".")
names(theta_df) <- params_species$species
colmax <- rep(0, 20)
colmin <- rep(0, 20)
for (h in 1:nb_class)
{
  max_species_group[,3 * h - 2] <- names(sort(colMeans(theta_df[tSNE_site_group == h,]), decreasing = TRUE)[1:20])
  min_species_group[,3 * h - 2] <- names(sort(colMeans(theta_df[tSNE_site_group == h,]), decreasing = FALSE)[1:20])
  for (hh in 1:dim(max_species_group)[1]) 
  {
    colmax[hh] <- which(names(theta_df) == max_species_group[hh, 3 * h - 2])
    colmin[hh] <- which(names(theta_df) == min_species_group[hh, 3 * h - 2])
    max_species_group[hh, 3 * h - 2] <- gsub("\ \ ", "\\.\ ", gsub("\\.", "\ ", max_species_group[hh, 3 * h - 2]))
    min_species_group[hh, 3 * h - 2] <- gsub("\ \ ", "\\.\ ", gsub("\\.", "\ ", min_species_group[hh, 3 * h - 2]))
  }
  max_species_group[, 3 * h - 1] <- colSds(data.matrix(theta_df)[tSNE_site_group == h,], cols = colmax)
  min_species_group[, 3 * h - 1] <- colSds(data.matrix(theta_df)[tSNE_site_group == h,], cols = colmin)
  max_species_group[, 3 * h] <- colMeans(data.matrix(theta_df)[tSNE_site_group == h,])[colmax]
  min_species_group[, 3 * h] <- colMeans(data.matrix(theta_df)[tSNE_site_group == h,])[colmin]
  writeLines(max_species_group[, 3 * h - 2], here("output", "classification", paste0("Group_", h, "_species_max.txt")))
  writeLines(min_species_group[, 3 * h - 2], here("output", "classification", paste0("Group_", h, "_species_min.txt")))
}
colnames(max_species_group) <- c(matrix(c(paste("Group", 1:nb_class), paste("Sd", 1:nb_class),
                                          paste("Mean", 1:nb_class)), ncol = nb_class, byrow = TRUE))
colnames(min_species_group) <- c(matrix(c(paste("Group", 1:nb_class), paste("Sd", 1:nb_class), 
                                          paste("Mean", 1:nb_class)), ncol = nb_class, byrow = TRUE))
write_csv2(data.frame(max_species_group), file = here("output", "classification", "max_species_group.csv"))
write_csv2(data.frame(min_species_group), file = here("output", "classification", "min_species_group.csv"))

##================
##
## Plot random effect alpha with
## size of the place
##================

alpha <- read_stars(here("output", "alphas.tif"))
latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
alpha_site <- rep(0, nrow = dim(coord)[1])
for (i in 1:dim(coord)[1]){
  alpha_site[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "alphas.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
data <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep = ",")
aire_site <- rep(0, dim(coord)[1])
k <- 1
for (i in unique(data$plot_name)){
  tempo <- data$AREA[data$plot_name == i][1]
  if (length(tempo) > 1){print(i)}
  aire_site[k] <- tempo[1]
  k <- k +1
}
save(aire_site, alpha_site, file = here("output", "RData", "alpha_vs_area.RData"))
load(here("output", "RData", "alpha_vs_area.RData"))

png(here("output", "plot", "alpha_vs_area.png"))
par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
plot(aire_site, alpha_site, main = "Site random effect \n as a function of site's area")
dev.off()

##================
##
## Specificity ~0.989
## Sensibility ~0.512
##
##================

npart <- 30
PA <-  read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
load(here("output", "RData", "jSDM_binom_pro.RData"))
theta <- jSDM_binom_pro$theta_latent

# Accuracy 
# Compute the proportion of species predicted as present 
# among the species  observed at each site for the different data-sets.
Sensitivity <- function(PA, theta){
  n_sites <- nrow(PA)
  score <- rep(0, n_sites)
  for(i in 1:n_sites){
    # Sensitivity 
    # True positive
    obs_sp <- which(PA[i,] > 0)
    nobs_sp <- length(obs_sp)
    pred_sp <- which(theta[i,] >= sort(theta[i,], decreasing = TRUE)[nobs_sp])
    score[i] <- sum(pred_sp %in% obs_sp) / ifelse(nobs_sp != 0, nobs_sp, 1)
  }
  return(score)
}
Specificity  <- function(PA, theta){
  n_sites <- nrow(PA)
  score <- rep(0, n_sites)
  for(i in 1:n_sites){
    # Specificity
    # True negative
    abs_sp <- which(PA[i,] == 0)
    nabs_sp <- length(abs_sp)
    pred_abs_sp <- which(theta[i,] <= sort(theta[i,])[nabs_sp])
    score[i] <- sum(pred_abs_sp %in% abs_sp) / ifelse(nabs_sp !=0, nabs_sp, 1)
  }
  return(score)
}

sens <- Sensitivity(PA, theta)
speci <- Specificity(PA, theta)
coord <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
df_sens <- data.frame(cbind(sens = sens, latitude = as.numeric(coord[,3]), longitude = as.numeric(coord[,4])))
df_speci <- data.frame(cbind(speci = speci, latitude = as.numeric(coord[,3]), longitude = as.numeric(coord[,4])))
df_TSS <- data.frame(cbind(TSS = speci + sens - 1, latitude = as.numeric(coord[,3]), longitude = as.numeric(coord[,4])))
save(df_sens, df_speci, file = here("output", "RData", "sens_speci.RData"))

load(here("output", "RData", "sens_speci.RData"))
NC <- ne_countries(scale = 10, returnclass = "sf")

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = df_sens, aes(latitude, longitude, colour = sens)) +
  scale_colour_gradientn(colours = rev(rocket(10)), name = "sensitivity", limits = c(0, 1)) +
  ggtitle("Sensitivity for each site") +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "sensitivity_by_site.png"))

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = df_speci, aes(latitude, longitude, colour = speci)) +
  scale_colour_gradientn(colours = rev(rocket(10)), name = "specificity", limits = c(0.9, 1)) +
  ggtitle("Specificity for each site") +
  theme(legend.position = c(0.2, 0.3)) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = here("output", "plot", "specificity_by_site.png"))

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = df_TSS, aes(latitude, longitude, colour = TSS)) +
  scale_colour_gradientn(colours = rev(rocket(10)), name = "TSS") +
  # ggtitle("TSS for each site") +
  theme(legend.position = c(0.2, 0.3)) +
  coord_sf(xlim = c(163.5, 168.5 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = here("output", "plot", "TSS_by_site.png"))

##================
##
## Check if variables are usefull
##
##================
load(here("output", "RData", "jSDM_binom_pro.RData"))
beta_check <- matrix(0, ncol = dim(jSDM_binom_pro$model_spec$site_data)[2] + 1 + jSDM_binom_pro$model_spec$n_latent, nrow = dim(PA)[2])
names(beta_check) <- c("beta_(Intercept)", "beta_ultramafic", "beta_bio1", "beta_bio4", "beta_bio12", "beta_bio15", "beta_cwd", "beta_bio1^2",
                       "beta_bio4^2", "beta_bio12^2", "beta_bio15^2", "beta_cwd^2", "beta_log(aire)", "lambda_1", "lambda_2")
for (i in 1:dim(PA)[2]){
  beta_check[i, ] <- colMeans(jSDM_binom_pro$mcmc.sp[[i]])
}
# variable's coef are differents for each species ie all variables have effect on prediction
for (j in 1:ncol(beta_check)) {
  print(names(beta_check)[j])
  print(summary(beta_check[,j]))
}

##================
##
## Plot group number on inventory site
##
##================

coord = read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")[,3:4]
group_site = matrix(0, ncol = 3, nrow = dim(coord)[1])
for (i in 1:dim(coord)[1]) # some Nan because not in "forest"
{
  group_site[i,] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "RGB_forest_group.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
coord$group <- rowSums(round(group_site))
coord$group[is.nan(coord$group)] <- 0
coord$longitude <- as.numeric(coord$longitude)
coord$latitude <- as.numeric(coord$latitude)
for (i  in 1:6) {
  coord$group[coord$group == sum(col2rgb(color_group_hex[i]))] <- i
}
coord$group <- as.factor(coord$group)
# none site are on group 6
NC <- ne_countries(scale = 10, returnclass = "sf")
nb_class <- 5

ggplot() +
  geom_sf(data = NC) +
  geom_point( aes(x = coord$longitude, y = coord$latitude, color = coord$group, shape = coord$group)) +
  scale_color_manual(values = c("#000000", color_group_hex[1:5]), 
                     label = c("0" = "Hors forêt", "1" = "Groupe 1", "2" = "Groupe 2","3" = "Groupe 3","4" = "Groupe 4", "5" = "Groupe 5")) +
  scale_shape_manual(values = c(3, rep(16, 5)),
                     label = c("0" = "Hors forêt", "1" = "Groupe 1", "2" = "Groupe 2","3" = "Groupe 3","4" = "Groupe 4", "5" = "Groupe 5")) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme_bw() +
  labs(y = "latitude", x = "longitude", color = "Groupe", shape = "Groupe") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
ggsave(here("output", "plot", "Observed_occurences.png"))
Sys.time()
