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
library(tidyterra)

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

reproject_sf_stars <- function(EPSG = 3163, sf_stars_object){
  #' Reproject shapefile with right EPSG
  #' 
  #' @param EPSG target EPSG, default 3163
  #' @param sf_stars_object sf or stars object to reproject
  #' @return same object then  `sf_stars_object` with `EPSG` projection
  
  sf_stars_object <- st_transform(sf_stars_object, EPSG)
 
  return(sf_stars_object)
}
 
create_sf_forest <- function(EPSG = 3163, altitude_min = 0, percentage_min = 50, elevation, forest){
  #' Create sf object of forest with minimal altitude and percentage of forest
  #'
  #' @param EPSG int. target EPSG, default 3163
  #' @param altitude_min float. minimal altitude to consider forest, default 0
  #' @param percentage_min int. percentage minimal of forest to consider pixels as forest, default 50
  #' @param elevation stars object. output of "read_stars" of "stars" library with elevation in meters
  #' @param forest stars object. output of "read_stars" of "stars" library with percentage of forest
  #' @return sf object. with only forest higher than `altitude_min` & `percentage_min`
  
  elevation <- reproject_sf_stars(EPSG = EPSG, sf_object = elevation)
  elevation[[1]][elevation[[1]] < 10] <- NA
  forest <- reproject_sf_stars(EPSG = EPSG, sf_object = forest)
  forest[forest < 50] <- NA
  forest <- st_crop(forest, elevation)
  forest <- st_as_sf(forest)
}

create_sf_moist_forest <- function(EPSG = 3163, forest, monthly_precipitation){
  #' Create sf object of moist forest according to FAO's criteria. 
  #'  MAP >=1500 mm
  #' <5 months where precipitations < 100 mm
  #'
  #' @param EPSG int. target EPSG, default 3163
  #' @param forest sf object. output of "read_sf" of "sf" library with boolean values of forest.
  #' @param monthly_precipitation multilayers stars object. output of "read_stars" of "stars" library with montly mean of precipitation in mm.
  #' @return sf object. with only moist forest area.
  
  forest <- reproject_sf_stars(EPSG = EPSG, sf_object = forest)
  monthly_precipitation <- reproject_sf_stars(EPSG = EPSG, sf_object = monthly_precipitation)
  if (length(names(st_dimensions(monthly_precipitation))) != 2){
    monthly_precipitation <- split(monthly_precipitation)
  }
  annual_precipitation <- monthly_precipitation[1,,]
  for (i in 2:12) { 
    annual_precipitation <- annual_precipitation + monthly_precipitation[i,,]
  }
  MAP <- monthly_precipitation[1,,]
  MAP[[1]] <- 0
  MAP[[1]][annual_precipitation[[1]] < 1500] <- 1 # dry
  moist_month <- monthly_precipitation[1,,]
  moist_month[[1]] <- 0
  for (i in 1:12) {
    moist_month[[1]] <- moist_month[[1]] + monthly_precipitation[[1]][,,i] < 100
  }
  moist_forest <- moist_month * MAP
  moist_forest[[1]][moist_forest[[1]] > 5] <- NA # remove no moist forest
  moist_forest <- st_as_sf(moist_forest)
  return(moist_forest)
}

# Merge climate and environ in one tif file
write_stars(c(st_normalize(st_crop(read_stars(here("output", "environ_allNC.tif")), forest)), read_stars(here("output", "current_chelsaNC.tif")), 
              along = "band"), dsn = here("output", "dataNC.tif"), 
            options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)

extract_var_JSDM <- function(stars_object, variables_names, power_variable = rep(1, length(variables_names)), 
                             scale = rep(TRUE, length(variables_names))){
  #' Extract layers and create power wanted of each variable 
  #' 
  #' @param stars_object stars object. output of "read_stars" of "stars" library with layer's names.
  #' @param variables_names a character vector. of layer's name needed
  #' @param power_variable int vector same length then `variables_names` . maximal power needed for each variable 
  #' @param scale boolean vector same length then `variables_names` . scale variables and its power
  #' @return stars object with only variables and power asked.
  
  if( length(variables_names) != length(power_variable)){
    print("variables_names length is different then power_variable length")
    break
  }
  if (length(names(st_dimensions(stars_object))) != 2){
    stars_object <- split(stars_object)
  }
  stars_output <- stars_file[variable_names,,]
  for (i in 1:length(variables_names)) {
    for (power in 1:power_variable[i]) {
      if ( power != 1){
        var_power <- stars_output[variables_names[i],,]^power
        if (scale[i]){
          var_scale <- scale(c(var_power[[1]]), na.rm = TRUE)
          var_power[[1]] <- (var_power[[1]] - attr(var_scale, "scaled:center")) / attr(var_scale, "scaled:scale")
        }
        stars_output <- c(stars_output,var_power)
        names(stars_output) <- names(stars_output, paste0(variables_names[i], "^", power))
      }
    }
    if (scale[i]){
      var_scale <- scale(c(stars_output[variables_names[i],,][[1]]), na.rm = TRUE)
      stars_output[variables_names[i],,][[1]] <- (stars_output[variables_names[i],,][[1]] - attr(var_scale, "scaled:center")) / attr(var_scale, "scaled:scale")
    }
  }
  return(stars_output)
}

data_JSDM <- function(EPSG = 3163, latlon_site, area_site, path_tiff, log_area_site = TRUE){
  #' Create dataframe with values of each variable on each site
  #'
  #' @param EPSG int. target EPSG, default 3163
  #' @param latlon_site float matrix with 2 columns. Columns are latitude and longitude in this order.
  #' @param area_site float vector same dimension then number of row `latlon`. area of each site in m²
  #' @param path_tiff character. path to Tiff file
  #' @param log_area_site boolean. Replace area of each site by log area of each site. Anyway, this variable is scale.
  #' @return matrix with latitude, longitude and all variables. If `latlon_site` has row names, output has same row name. 
  
  if (dim(latlon_site)[1] != length(area_site)){
    print("number of rows in latlon is different than lengh of area_site")
    break
  }
  stars_object <- read_stars(path_tiff)
  if (length(names(st_dimensions(stars_object))) != 2){
    stars_object <- split(stars_object)
  }
  if (log_area_site){
    area_site <- scale(log(area_site))
    name_area <- "log_area_site"
  }else{
    area_site <- scale(area_site)
    name_area <- "area_site"
  }
  data_site <- matrix(0, ncol = length(names(stars_object)), nrow = length(area_site))
  row.names(data_site) <- row.names(latlon_site)
  proj.t <- paste0("EPSG:", EPSG) 
  for (i in 1:dim(latlon_site)[1]) {
    data_site[i, 1:length(names(stars_object))] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {path_tiff} \\
                                                                       -wgs84 {latlon_site[i,1]} {latlon_site[i,2]}'), intern = TRUE))
  }
  data_site <- cbind(data_site, area_site)
  names(data_site) <- c(names(stars_object), name_area)
  return(data_site)
}

##================
##
## jSDM Binomial Probit
## about 1h to run
##================

JSDM_bino_pro <- function(presence_data, site_data, n_latent = 0, V_beta, mu_beta, nb_species_plot = 2,
                          display_plot = TRUE, save_plot = FALSE){
  #' Run jSDM_binomial_probit from jSDM package and plot results for some species. Plots can be display and/or save.
  #'
  #' @param presence_data int matrix. Presence Absence for each species (col) on each site (row). Avoid empty row or column
  #' @param site_data float matrix. Explanatories variables for each site without missing values. Same number of row than number of row in `presence_data`
  #' @param n_latent int. number of latent variables to use in the model, default is 0.
  #' @param V_beta float vector. Variances of Normal priors for the beta parameters.
  #' @param mu_beta float vector. Means of Normal priors for the beta parameters.
  #' @param nb_species_plot number of species whose results are display.
  #' @param display_plot boolean. If TRUE, display all plot, default is TRUE.
  #' @param save_plot boolean. Save in local in .png all plot in folder plot, default is FALSE.
  
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
  
  top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[nb_species_plot])
  np <- nrow(jSDM_binom_pro$model_spec$beta_start)
  
  
  for (j in top_species[1]) {
    for (p in 1:np) {
      par(mfrow = c(1, 2))
      coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
      coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
      mtext(paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p], ", species : ", names(top_species[top_species == j])),
            side = 3, line = - 2, outer = TRUE)
      z = recordPlot()
      dev.off((1 - display_plot) * 2)
      if (save_plot){
        dir.create(here("plot"), showWarnings = FALSE)
        png(here("plot", paste0("beta_jSDM_", p, ".png")))
        replayPlot(z)
        dev.off()
      }
    }
  }
  
  
  ## lambda_j of the top five species
  n_latent <- jSDM_binom_pro$model_spec$n_latent
  
  for (j in top_species[1]) {
    
    par(mfrow = c(1, 2))
    for (l in 1:n_latent) {
      coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
      coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
      mtext(paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[np + l],", species : ", names(top_species[top_species == j])),
            side = 3, line = - 2, outer = TRUE)
      z = recordPlot()
      dev.off((1 - display_plot) * 2)
      if (save_plot){
        png(here("plot", "lambda_jSDM.png"))
        replayPlot(z)
        dev.off()
      }
    }
  }
  
  ## Latent variables W_i for the first two sites
  
  for (l in 1:n_latent) {
    par(mfrow = c(2, 2))
    for (i in 1:2) {
      coda::traceplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                      main = paste0("Latent variable W_", l, ", site ", i))
      coda::densplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                     main = paste0("Latent variable W_", l, ", site ", i))
    }
    z = recordPlot()
    dev.off((1 - display_plot) * 2)
    if (save_plot){
      png(here("plot", paste0("W_", l, "_jSDM.png")))
      replayPlot(z)
      dev.off()
    }
  }
  
  ## alpha_i of the first two sites
  plot(coda::as.mcmc(jSDM_binom_pro$mcmc.alpha[, 1:2]))
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "alpha_jSDM.png"))
    replayPlot(z)
    dev.off()
  }
  
  ## V_alpha
  par(mfrow = c(2, 2))
  coda::traceplot(jSDM_binom_pro$mcmc.V_alpha)
  coda::densplot(jSDM_binom_pro$mcmc.V_alpha)
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "V_alpha_Deviance_jSDM.png"))
    replayPlot(z)
    dev.off()
  }
  
  ## Deviance
  par(mfrow = c(1, 2))
  coda::traceplot(jSDM_binom_pro$mcmc.Deviance)
  coda::densplot(jSDM_binom_pro$mcmc.Deviance)
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "Deviance_jSDM.png"))
    replayPlot(z)
    dev.off()
  }
  
  ## probit_theta
  par (mfrow = c(2, 1))
  hist(jSDM_binom_pro$probit_theta_latent,
       main = "Predicted probit theta", xlab = "predicted probit theta")
  hist(jSDM_binom_pro$theta_latent,
       main = "Predicted theta", xlab = "predicted theta")
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "theta_jSDM.png"))
    replayPlot(z)
    dev.off()
  }
}

##================
##
## Plotting occurences for each species
##
##================

plot_pres_est_one_species <- function(jSDM_binom_pro, species_to_plot = colnames(jSDM_binom_pro$model_spec$presence_data)[1],
                                      coord_site, country_name = NULL, latlon_output = NULL, display_plot = TRUE, 
                                      save_plot = FALSE){
  #' Plot presence absence  and estimated probabilities of presence for one species. Possibility to save plots
  #'
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @param species_to_plot character. species to plot, default is the first species of the list.
  #' @param coord_site dataframe. columns name must contain latitude and longitude as names.
  #' @param country_name character. English name of the country where are inventory sites, default is NULL.
  #' @param latlon_output float vector. inconsistent with `country_name` coord of output box is this format {lon_min, lat_min, lon_max, lat_max}, default is NULL.
  #' @param display_plot boolean. show plot, default is TRUE.
  #' @param save_plot boolean. Write plot in plot folder as .png files, default is FALSE.
 
  if (is.null(country_name) & is.null(latlon_output)){
    print("Need latlon coordinates")
    break
  }
  if (!is.null(country_name)){
    latlon_ouput <- st_bbox(ne_countries(scale = 10, country = country_name, returnclass = "sf"))
  }
  map <- ne_countries(scale = 10, returnclass = "sf")
  df_PA <- jSDM_binom_pro$model_spec$presence_data
  coord_site[, species_to_plot] <- df_PA[, species_to_plot]
  pres <- data.frame(coord_site[coord_site[, species_to_plot] == 1, c("longitude", "latitude")])
  abs <- data.frame(coord_site[coord_site[, species_to_plot] == 0, c("longitude", "latitude")])
  ggplot(data = map) +
    geom_sf() +
    geom_point(data = coord_site[,c("longitude", "latitude", species_to_plot)], 
               aes(x = as.numeric(longitude), y = as.numeric(latitude), shape = as.logical(coord_site[,species_to_plot]),
                   color = as.logical(coord_site[,species_to_plot])), size = 1) +
    ggtitle(paste0('Observed occurences of \n', species_to_plot)) +
    scale_shape_manual(labels = c("Absence", "Presence"), values = c(1, 3), name = "Species state") +
    scale_color_manual(labels = c("Absence", "Presence"), values = c("red", "black"), name = "Species state") +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    dir.create(here("plot"))
    ggsave(here("plot", paste0("Observed_occurences_", species_to_plot, ".png")))
  }
  
  coord_site$theta <- jSDM_binom_pro$theta_latent[, species_to_plot]
  ggplot(data = map) +
    geom_sf() +
    geom_point(data = data_theta, 
               aes(x = as.numeric(longitude), y = as.numeric(latitude), colour = theta), size = 1) +
    ggtitle(paste0("Estimated probabilities of presence for \n", species_to_plot)) +
    scale_colour_gradientn(colours = rev(rocket(10)), name = "Prediction", limits = c(0, 1)) +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    ggsave(here("plot", paste0("estimated_probabilities_", species_to_plot, ".png")))
  }
}

##================
##
##Plotting Species Richness
##
##================

plot_species_richness <- function(jSDM_binom_pro, coord_site, country_name = NULL, 
                                  latlon_output = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Plot species richnessobserved and  estimated probabilities of presence for all sites. Possibility to save plots
  #'
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @param coord_site dataframe. columns name must contain latitude and longitude as names.
  #' @param country_name character. English name of the country where are inventory sites, default is NULL.
  #' @param latlon_output float vector. inconsistent with `country_name` coord of output box is this format {lon_min, lat_min, lon_max, lat_max}, default is NULL.
  #' @param display_plot boolean. show plot, default is TRUE.
  #' @param save_plot boolean. Write plot in plot folder as .png files, default is FALSE.
  
  if (is.null(country_name) & is.null(latlon_output)){
    print("Need latlon coordinates")
    break
  }
  if (!is.null(country_name)){
    latlon_ouput <- st_bbox(ne_countries(scale = 10, country = country_name, returnclass = "sf"))
  }
  df_PA <- jSDM_binom_pro$model_spec$presence_data
  species_richness <- data.frame(richness = rowSums(df_PA))
  species_richness$latitude <- coord_site$latitude
  species_richness$longitude <- coord_site$longitude
  species_richness$estimate <- rowSums(jSDM_binom_pro$theta_latent)
  map <- ne_countries(scale = 10, returnclass = "sf")
  ggplot(data = map) +
    geom_sf() +
    geom_point(data = species_richness, 
               aes(x = as.numeric(longitude), y = as.numeric(latitude), colour = richness), size = 1) +
    ggtitle("Observed species richness") +
    scale_colour_gradientn(colours = rev(rocket(10)), name = "Number of species", 
                           limits = c(floor(min(species_richness[, c(1,4)])), ceiling(max(species_richness[, c(1,4)])))) +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    dir.create(here("plot"))
    ggsave(here("plot", "species_richness_observed.png"))
  }
  
  # Estimated richness
  
  ggplot(data = map) +
    geom_sf() +
    geom_point(data = species_richness, 
               aes(x = as.numeric(longitude), y = as.numeric(latitude), colour = estimate), size = 1) +
    ggtitle("Estimated species richness") +
    scale_colour_gradientn(colours = rev(rocket(10)), name = "Number of species", 
                           limits = c(floor(min(species_richness[, c(1,4)])), ceiling(max(species_richness[, c(1,4)])))) +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    ggsave(here("plot", "species_richness_estimated.png"))
  }
  
  par(mfrow = c(1,1))
  plot(rowSums(PA), rowSums(jSDM_binom_pro$theta_latent), 
       main = "Species richness for area \n with less than 50 differents species",
       xlab = "observed", ylab = "estimated", 
       cex.main = 1.4, cex.lab = 1.3, xlim = c(0, 50), ylim = c(0, 50))
  abline(a = 0, b = 1, col = 'red')
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("output", "plot", "species_richness_observed_vs_estimated.png"), width = 480, height = 480, units = "px")
    replayPlot(z)
    dev.off()
  }
}

##================
##
## Spatial Interpolation for alpha and latents variables with Knn
##
##================

knn_interpolation_jSDM <- function(jSDM_binom_pro, k = 5, coord_site, sf_forest){
  #' Interpole alpha and all latent variables with knn algorithm. Then save interpolation on disk in output folder.
  #'
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @param k int. number of neighbours to consider in knn algorithm, default is 5.
  #' @param coord_site dataframe. columns name must contain latitude and longitude as names.
  #' @param sf_forest sf object. area where interpolation will be done.
  
  dir.create(here("output"))
  coord_pixel_latlon <- coordinates(spTransform(as_Spatial(forest), CRS("+proj=longlat +datum=WGS84")))
  # create stars object with forest's arguments
  stars_object <- matrix(0, nrow = (st_bbox(forest)[3] - st_bbox(forest)[1]) / 1000, 
                         ncol = (st_bbox(forest)[4] - st_bbox(forest)[2]) / 1000)
  stars_object <- st_as_stars(stars_object)
  stars_object <- st_set_crs(stars_object, st_crs(forest))
  stars_object <- st_set_bbox(stars_object, st_bbox(forest))
  st_dimensions(stars_object)$X1$delta <- 1000
  st_dimensions(stars_object)$X2$delta <- 1000
  stars_object <- st_crop(stars_object, forest)
  distance <- stars_object
  coord_pixel_df <- which(!is.na(stars_object[[1]]), arr.ind = TRUE)
  for (site in 1:dim(coord_site)[1]) # distance between pixels & each site ~16min
  {
    for (pixel in 1:dim(coord_pixel_latlon)[1]) 
    {
      distance[[1]][coord_pixel_df[pixel, 1], 
                    coord_pixel_df[pixel, 2]] <- sqrt((coord_pixel_latlon[pixel, 1] - as.numeric(coord_site[site, 1]))^2 + 
                                                        (coord_pixel_latlon[pixel, 2] - as.numeric(coord_site[site, 2]))^2)
    }
    if (site == 1)
    {
      all_distance <- distance
    }else{
      all_distance <- c(all_distance, distance)
    }
  }
  
  for (name_inter in c("alpha", paste0("lv_", 1:jSDM_binom_pro$model_spec$n_latent))) { # alpha + all latent variables
    interpolated_stars <- stars_object
    names(interpolated_stars) <- name_inter
    if ( name_inter == "alpha"){
      coord_site[, ncol(coord_site) + 1] <- colMeans(jSDM_binom_pro$mcmc.alpha)
      
    }else{
      coord_site[, ncol(coord_site) + 1] <- colMeans(jSDM_binom_pro$mcmc.latent[[name_inter]])
    }
    names(coord_site) <- c(names(coord_site), name_inter)
    
    for (site in 1:dim(coord_pixel_latlon)[1]) {
      dis_pixel <- all_distance[[1]][coord_pixel_stars[site, 1], coord_pixel_stars[site, 2],]
      knn_site_nb <- which(dis_pixel >= sort(dis_pixel, decreasing = TRUE)[knn]) # site number
      knn_dis <- dis_pixel[knn_site_nb] # distance to each site
      inter_value <- sum(coord_site[knn_site_nb, ncol(coord_site)] * (1 - knn_dis / sum(knn_dis)) / (length(knn_site_nb) - 1))
      interpolated_stars[[1]][coord_pixel_stars[site, 1], coord_pixel_stars[site, 2]] <- inter_value
    }
    if (name_inter == "alpha"){
      interpolated_stars[[1]] <- interpolated_stars[[1]] - mean(interpolated_stars[[1]], na.rm = TRUE)
    }else{
      interpolated_stars[[1]] <- scale(interpolated_stars[[1]])
    }
    write_stars(interpolated_stars, here("output", paste0(name_inter, "_knn.tif")), options = c("COMPRESS=LZW","PREDICTOR=2"))
  }
}

##================
##
## Compute probabilities for each species
##
##================

prob_est_species_forest <- function(alpha_stars, latent_var_stars, jSDM_binom_pro, data_stars, area_stars = NULL){
  #' Create tif files with probabilities of presences for each species on a given area. Results are save in output/theta
  #'
  #' @param alpha_stars stars object. with centered values in 0.
  #' @param latent_var_stars multilayer stars object. with as layer as latent variables. Make sure each latent variable is scale.
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library.
  #' @param data_stars multilayer stars object. with same explanatories variables whom in `jSDM_binom_pro`. Make sure your explanatories variable are scale.
  #' @param area_stars stars object. with values of pixel size with same scale parameters used in `jSDM_binom_pro`.
  
  dir.create(here("output"))
  dir.create(here("output", "theta"))
  if ( length(data_stars) == 1){
    data_stars <- split(data_stars)
  }
  if ( length(latent_var_stars) == 1){
    latent_var_stars <- split(latent_var_stars)
  }
  data_stars <- c(data_stars, area_stars)
  n_var <- length(data_stars)
  n_species <- length(colnames(jSDM_binom_pro$theta_latent))
  npart <- ceiling(n_species / 30) # number of species in each file set to 30
  
  n_latent <- length(jSDM_binom_pro$mcmc.latent)
  lambdas <- matrix(0, n_species, n_latent)
  betas <- matrix(0, n_species, n_var + 1 )
  for (j in 1:n_species){
    lambdas[j, ] <- colMeans(jSDM_binom_pro$mcmc.sp[[j]][, n_var + 2 + 1:n_latent ])
    betas[j, ] <- colMeans(jSDM_binom_pro$mcmc.sp[[j]][, 1:(n_var + 1)]) # n_var + intercept 
  }
  colnames(betas) <- names(data_stars)
  params_species <- data.frame(betas, lambda = lambdas)
  rownames(params_species) <- colnames(jSDM_binom_pro$model_spec$presence_data)
  predfun <- function(data_stars, params_species, alpha_stars, latent_var_stars, species.range){
    #' give prediction for each species
    #' in two step 
    #' First init for first species
    #' Second same process but stack on first species output (stars object)
    
    # Xbeta_1
    n_var <- length(data_stars)
    Xbeta1 <- alpha_stars
    Xbeta1[[1]] <- params_species[1, "beta_(Intercept)"]
    for (p in 1:n_var) {
      Xbeta_1 <- Xbeta_1 + data_stars[[p]] * params_species[1, p + 1] 
    }
    
    # Wlambda_1
    Wlambda_1 <- alpha_stars
    Wlambda_1[[1]] <- 0
    for (layer in 1: length(latent_var_stars)) {
      Wlambda_1 <- Wlambda_1 + latent_var_stars[layer,,] * params_species[1, paste0("lambda.", layer)] 
    }
    
    # probit_theta_1
    probit_theta_1 <- Xbeta_1 + Wlambda_1 + alpha_stars + area_stars
    probit_theta <- probit_theta_1
    
    ## Other species
    for (j in (species.range[1] + 1):species.range[2]) {
      
      ## Xbeta_j
      Xbeta_j <- Xbeta_1
      Xbeta_j[[1]] <- params_species[j, "beta_(Intercept)"]
      for (p in 1:(n_var - 1)) {
        Xbeta_j <- Xbeta_j + data_stars[[p]] * params_species[j, p + 1] 
      }
      
      ## Wlambda_j
      Wlambda_j <- alpha_stars
      Wlambda_j[[1]] <- 0
      for (layer in 1: length(latent_var_stars)) {
        Wlambda_j <- Wlambda_j + latent_var_stars[layer,,] * params_species[j, paste0("lambda.", layer)] 
      }
      
      ## probit_theta_j
      probit_theta_j <- Xbeta_j + Wlambda_j + alpha_stars + area_stars
      probit_theta <- c(probit_theta, probit_theta_j)
    }
    names(probit_theta) <- rownames(params_species)[species.range[1]:species.range[2]]
    return(probit_theta)
  }
  
  first.species <- seq(1, n_species, by = floor(n_species / npart) + 1)
  for (n in 1:npart){
    probit_theta <- predfun(data_stars, params_species, alpha_stars, latent_var_stars,
                            species.range = c(first.species[n], min(n_species, first.species[n] + floor(n_species / npart))))
    theta <- probit_theta
    for (j in 1:length(theta)) {
      theta[[1]][j,,] <- pnorm(theta[[1]][j,,])
    }
    dir.create(here("output"), showWarnings = FALSE)
    dir.create(here("output", "theta"), showWarnings = FALSE)
    write_stars(theta, options = c("COMPRESS=LZW", "PREDICTOR=2"), 
                dsn = here("output", "theta", paste0("KNN_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  }
}

##================
##
## Predictive maps of presence probabilities
##
##================

plot_prob_pres_interp <- function(theta_stars, species_to_plot = names(theta_stars)[1], country_name = NULL, country_sf = NULL,
                                  display_plot = TRUE, save_plot = FALSE){
  #' Create plot with probabilities of presence interpolated for a choosen species. Plot can be hide and/or save in plot folder
  #' 
  #' @param theta_stars multilayer stars_object. from files save in prob_est_species_forest function.
  #' @param species_to_plot character. name of one layer of `theta_stars`, default is first one.
  #' @param country_name character. optional, for display border of country on map, default is NULL.
  #' @param country_sf sf object. optional, inconsistent with `country_name`. Display borders on map, default is NULL.
  #' @param display_plot boolean. Display plot, default is TRUE.
  #' @param save_plot boolean. Save plot in .png format in folder plot, default is FALSE.
  
  if (!is.null(country_name) & !is.null(country_sf)){
    print("Only country_name or country_sf")
    break
  }
  gplot <- ggplot()  
  if(!is.null(country_name)){
    gplot <- gplot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                             colour = "black", fill = "grey") 
  }
  if(!is.null(country_sf)){
    gplot <- gplot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
  }
  gplot <- gplot + geom_stars(data = theta) +
    ggtitle(paste('Interpolated current probabilities of presence \n for', species_to_plot)) +
    scale_fill_gradientn(colours = rev(rocket(9)), na.value = "transparent") +
    coord_fixed() +
    theme_bw() +
    labs(fill = "Number of species") +
    theme(plot.title = element_text(hjust = 0.5))
  gplot
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(plot = gplot, here("plot", paste0("interpolated_presence_for", species_to_plot, ".png")))
  }
}

##================
##
## Predictive maps of species richness
##
##================

plot_species_richness_interpolated <- function(list_theta_path, country_name = NULL, 
                                               country_sf = NULL, display_plot = TRUE, save_plot = FALSE, save_tif = FALSE){
  #' Create plot with species richness with optional sf border.
  #'
  #' @param list_theta_path character vector. with full path of Tiff files who contain probabilities of presence of each species to consider.
  #' @param country_name character. optional, for display border of country on map, default is NULL.
  #' @param country_sf sf object. optional, inconsistent with `country_name`. Display borders on map, default is NULL.
  #' @param display_plot boolean. Display plot, default is TRUE.
  #' @param save_plot boolean. Save plot in .png format in folder plot, default is FALSE.
  #' @param save_tif boolean. Save species richness in .tif file in folder output, default is FALSE.
  #' @return terra object. species richness in designated area.
  
  theta_sum <- sum(terra::rast(list_theta_path[1]))
  n_files <- length(list_theta_path)
  if (n_files > 1){
    for (i in list_theta_path[2:n_files])
    {
      theta_sum <- theta_sum + sum(terra::rast(i))
    }
  }
  if (save_tif){
    dir.create(here("output"), showWarnings = FALSE)
    terra::writeRaster(theta_sum, here("output", "species_richness.tif"), overwrite = TRUE)
  }
  gplot <- ggplot()
  if(!is.null(country_name)){
    gplot <- gplot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                             colour = "black", fill = "grey") 
  }
  if(!is.null(country_sf)){
    gplot <- gplot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
  }
  gplot <- gplot + geom_stars(data = theta_sum) +
    ggtitle("Estimated current species richness") +
    scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
    coord_fixed() +
    theme_bw() +
    labs(fill = "Number of species") + 
    theme(plot.title = element_text(hjust = 0.5))
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(here("plot", "estimated_species_richness_forest.png"))
  }
  return(theta_sum)
}

##================
##
## Init PCA and tSNE
## 
##================
PCA_or_tSNE_on_pro_pres_est <- function(tSNE = TRUE, PCA = FALSE, list_theta_path, var_exp_stars = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Run PCA or tSNE on pixels with as coordinates, probabilities of presence of species. Plot Axis 1 vs axis 2, axis 2 vs axis 3 and dendrogram. Possibility to hide plot and save it in .png.
  #'
  #' @param tSNE boolean. compute tSNE reduction of dimensions algorithm. If TRUE, `PCA` must be FALSE, default is TRUE.
  #' @param PCA boolean. compute PCA reduction of dimensions algorithm. If TRUE, `tSNE` must be FALSE, default is FALSE.
  #' @param list_theta_path character vector.  with full path of Tiff files who contain probabilities of presence of each species to consider.
  #' @param var_exp_stars multilayer stars object. optional, only consider if `PCA` = TRUE. Plot as supplementary variables, default is NULL.
  #' @param display_plot boolean. Show plot, default is TRUE.
  #' @param save_plot boolean. Save plots in plot folder with png format.
  #' @return float matrix. with as much row as number of pixels and 3 columns for tSNE algorithm and as columns as number of axis with more than 1% of explained variance for PCA.
  
  if (tSNE * PCA){
    print("Chose only one among tSNE and PCA")
    break
  }
  cell <- which(!is.na(split(read_stars(list_theta_path[1]))[[1]]), arr.ind = TRUE)[, 1:2] # keep only x & y axis
  matrix_values <- NULL
  for (i in list_theta_path) { 
    # create matrix with probabilities for all species on each pixels
    matrix_values <- cbind(matrix_values, values(rast(i)))
  }
  if (tSNE){
    tSNE_fit <- Rtsne(scale(matrix_values), dims = 3, max_iter = 5000, num_threads = 0,
                      verbose = TRUE)
    pixel_df <- as.data.frame(tSNE_fit$Y)
    dist_pixel <- dist(pixel_df)
    axis1_2 <- ggplot(data = tSNE_df, aes(x = V1, y = V2)) +
      geom_point() +
      theme(legend.position = "bottom")
    axis2_3 <- ggplot(data = tSNE_df, aes(x = V2, y = V3)) +
      geom_point() +
      theme(legend.position = "bottom")
    if (display_plot){
      axis1_2
      axis2_3
    }
    if (save_plot){
      dir.create(here("plot"), showWarnings = FALSE)
      ggsave(here("plot", "tSNE_axis1_axis2.png"), plot = axis1_2)
      ggsave(here("plot", "tSNE_axis2_axis3.png"), plot = axis2_3)
    }
  }
  
  if (PCA){
    if (!is.null(var_exp_stars)){
      if (length(var_exp_stars) != 1){
        var_exp_stars <- merge(var_exp_stars)
      }
      matrix_var_exp <- matrix(var_exp_stars[[1]][!is.na(var_exp_stars[[1]])], ncol = length(split(var_exp_stars)))
    }
    pca_with_env <- PCA(cbind(matrix_values, matrix_var_exp), graph = FALSE, ncp = 10,
                        quanti.sup = dim(matrix_values)[2] + 1:length(split(var_exp_stars)))
    dist_pixel <- dist(get_pca_ind(pca_with_env)$coord, method = "euclidian")
    pixel_df <- get_pca_ind(pca_with_env)$coord[, 1: which.min(fviz_eig(pca_with_env)$data[,2][fviz_eig(pca_with_env)$data[,2] > 1])]
    exp_variances <- fviz_eig(pca_with_env)
    PCA_1_2 <- fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), fill.ind = "orange", repel = TRUE, 
                               col.quanti.sup = "black", invisible = "var", pointshape = 21, labelsize = 7, axes = c(1, 2),
                               col.ind = "NA", alpha = 0.8) +
        ggtitle("PCA on pixels with explanatories variables in sup") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    PCA_2_3 <- fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), fill.ind = "orange", repel = TRUE, 
                               col.quanti.sup = "black", invisible = "var", pointshape = 21, labelsize = 7, axes = c(2, 3),
                               col.ind = "NA", alpha = 0.8) +
      ggtitle("PCA on pixels with explanatories variables in sup") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    
    if (display_plot){
    exp_variances 
      print(paste0("Explained variances by first ", which.min(fviz_eig(pca_with_env)$data[,2][fviz_eig(pca_with_env)$data[,2] > 1]),
                   " axis : ", 
                   round(sum(fviz_eig(pca_with_env)$data[,2][fviz_eig(pca_with_env)$data[,2] > 1])), 
                   "%"))
      PCA_1_2
      PCA_2_3
    }
    if (save_plot){
      png(here("plot", "PCA_explained_variances.png"))
      exp_variances
      dev.off()
      ggsave(here("plot", "PCA_with_environ_sup_axis_1_2.png"), plot = PCA_1_2)
      ggsave(here("plot", "PCA_with_environ_sup_axis_1_2.png"), plot = PCA_2_3)
    }
  }
  
  # Plot CAH for choosen method
  plot(hclust(dist_pixel), labels = FALSE, main = "", xlab = "", ylab = "", sub = "", ylim = "none") 
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    method = c(PCA, tSNE)
    names(method) <- c("PCA", "tSNE")
    png(here("plot", paste0("dendro_", names(method)[max(method)], ".png")))
    plot(hclust(dist_tSNE), labels = FALSE, main = "", xlab = "", ylab = "", sub = "", ylim = "none") 
    dev.off()
  }
  return(pixel_df)
}

##===========
##
## K-means with t-SNE coords
## plot KM groups on t-SNE
##===========

plot_HCA_EM_Kmeans <- function(method = "KM", pixel_df, nb_group, plot_3d = FALSE, display_plot_2d = TRUE, save_plot_2d = FALSE){
  #' Display plots in 2 & 3D with number of groups set and a choosen method. Possibility to hide 2D plots and save it.
  #'
  #' @param method character. clustering method to chose among c("HCA", EM, "KM"). See `details` for more explications. Default is KM.
  #' @param pixel_df dataframe. output of "dist" function in "stats" library or "PCA_or_tSNE_on_pro_pres_est" function.
  #' @param nb_group int. number of groups for clustering. Please set nb_group to high than 1.
  #' @param plot_3d boolean. Display new window with 3d plot, default is FALSE.
  #' @param display_plot_2d boolean. Display 2d plots, default is TRUE.
  #' @param save_plot_2d boolean. Save in plot folder as .png file all 2d plots, default is FALSE.
  #' @return int vector. Group number for each pixel.  
  #' @details `method` values are HCA : Hierarchical Cluster Analysis; EM : Expectation Maximisation; KM : K-means algorithm
  
  if (nb_group == 1){
    print("Please select more than 1 group")
    break
  }
  if (method == "KM"){
    KM_class <- kmeans(pixel_df, centers = nb_group, iter.max = 20, nstart = 1000)
    pixel_group <- KM_class$cluster
  }
  if (method == "EM"){
    init_EM <- rand.EM(pixel_df, nclass = nb_group, min.n = 10)
    EM_class <- assign.class(pixel_df, emcluster(pixel_df, init_EM))
    pixel_group <- EM_class$class
  }
  if (method == "CAH"){
    dist_pixel <- dist(pixel_df)
    pixel_group <- cutree(hclust(dist_pixel), k = nb_group)
  }
  colnames(pixel_df) <- paste0("V", 1:dim(pixel_df)[2])
  plot1_2 <- ggplot(data = data.frame(pixel_df, pixel_group),
                    aes(x = V1, y = V2, color = as.factor(pixel_group))) +
    scale_color_manual(values = viridis(nb_group)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
  plot2_3 <- ggplot(data = data.frame(pixel_df, pixel_group),
                    aes(x = V2, y = V3, color = as.factor(pixel_group))) +
    scale_color_manual(values = viridis(nb_group)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
  if (display_plot_2d){
    plot1_2
    plot2_3
  }
  if (save_plot_2d){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(here("plot", paste0("axis_1_2_", method, "_groups.png")), plot = plot1_2)
    ggsave(here("plot", paste0("axis_2_3_", method, "_groups.png")), plot = plot2_3)
  }
  if (plot_3d){
    col_3d <- viridis(nb_group)[pixel_group]
    par3d(windowRect = c(20, 30, 800, 800))
    plot3d(pixel_df$V1, pixel_df$V2, pixel_df$V3, col = col_3d, xlab = "x", ylab = "y", zlab = "z", size = 5, pch = 16)
    legend3d("bottomright", legend = paste0("Group_", 1:length(unique(pixel_group))), pch = 16, col = viridis(nb_group),
             cex = 2, inset = c(0.02))
  }
  return(pixel_group)
}

##=====
## Change the coordinate scale for [10.255].
## for plot color for each pixel
##=====
plot_RGB_group_by_color <- function(stars_pixels_group, coord_pixel, country_name = NULL, 
                                    country_sf = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Create plots with coordinates from dimensional reduction algorithm and plot pixel group with multiple colors. Plot can be hide and/or save.
  #'
  #' @param stars_pixel_group stars object. with group number in each pixel usefull, NA in others pixels.
  #' @param coord_pixel flot matrix. coordinates of each pixel, obtained with `PCA_or_tSNE_on_pro_pres_est` function.
  #' @param country_name character. optional, add country's border to plots , default is NULL.
  #' @param country_sf sf object. optional, add border to plots, default is NULL.
  #' @param display_plot boolean. Show plots, default is TRUE.
  #' @param save_plot boolean. Save plots in plot folder as .png file, default is FALSE.
  
  nb_group <- length(unique(c(stars_pixels_group[[1]])))
  coord_pixel <- coord_pixel[, 1:3]
  cells_pixel <- which( !is.na(stars_pixels_group), arr.ind = TRUE)[, 1:2]
  R_stars <- G_stars <- B_stars <- stars_pixels_group
  for (l in 1:3) { 
    # set to [0,255]
    coord_pixel[, l] <- coord_pixel[, l] - min(coord_pixel[, l])
    coord_pixel[, l] <- coord_pixel[, l] / max(coord_pixel[, l]) * 255
  }
  for (cell in 1:dim(cells_pixel)[1]) { 
    # affect values in [0, 255] to each pixel in stars object 
    R_stars[[1]][cells_pixel[cell, 1], cells_pixel[cell, 2]] <- coord_pixel[cell, 1]
    G_stars[[1]][cells_pixel[cell, 1], cells_pixel[cell, 2]] <- coord_pixel[cell, 2]
    B_stars[[1]][cells_pixel[cell, 1], cells_pixel[cell, 2]] <- coord_pixel[cell, 3]
  }
  # RGB plot 
  RGB_terra <- terra::rast(c(R_stars, G_stars, B_stars))
  RGB_terra <- colorize(RGB_terra, value = 1:3, to = "col", stretch = "hist")
  col_palette <- coltab(RGB_terra)[[1]][, 1:4]
  col_hex <- rep(0, dim(col_palette)[1])
  for (i in 1:dim(col_palette)[1]) {
    col_hex[i] <- rgb2hex(r = col_palette[i, 2], g = col_palette[i, 3], b = col_palette[i, 4])
  }
  RGB_stars <- st_as_stars(RGB_terra)
  # plot groups with mean color of each group
  col_mean_group <- rep(0, nb_group)
  for (group in 1:nb_group){
    col_mean_group[group] <- rbg2hex(r = mean(R_stars[[1]][which(stars_pixels_group[[1]] == group, arr.ind = TRUE)]),
                                     g = mean(G_stars[[1]][which(stars_pixels_group[[1]] == group, arr.ind = TRUE)]),
                                     b = mean(B_stars[[1]][which(stars_pixels_group[[1]] == group, arr.ind = TRUE)]))
  }
  names(col_mean_group) <- 1:length(col_mean_group) 
  # init ggplot
  R_plot <- ggplot()
  G_plot <- ggplot()
  B_plot <- ggplot()
  RGB_plot <- ggplot()
  group_plot <- ggplot()
  if(!is.null(country_name)){
    R_plot <- R_plot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                               colour = "black", fill = "grey") 
    G_plot <- G_plot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                               colour = "black", fill = "grey") 
    B_plot <- B_plot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                               colour = "black", fill = "grey") 
    RGB_plot <- RGB_plot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                                   colour = "black", fill = "grey") 
    group_plot <- group_plot + geom_sf(data =  ne_countries(scale = 10, returnclass = "sf", country = country_name), 
                                       colour = "black", fill = "grey") 
  }
  if(!is.null(country_sf)){
    R_plot <- R_plot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
    G_plot <- G_plot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
    B_plot <- B_plot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
    RGB_plot <- RGB_plot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
    group_plot <- group_plot + geom_sf(data = country_sf, colour = "black", fill = "grey") 
  }
  R_plot <- R_plot + geom_stars(data = R_stars) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"), na.value = "transparent") +
    coord_fixed() +
    theme_bw()
  G_plot <- G_plot +geom_stars(data = G_stars) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Greens"), na.value = "transparent") +
    coord_fixed() +
    theme_bw()
  B_plot <- B_plot +geom_stars(data = B_stars) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Blues"), na.value = "transparent") +
    coord_fixed() +
    theme_bw()
  title <- ggdraw() + 
    draw_label("Pixel coordinates of the first 3 axis \n with color variations")
  R_G_B_plot <- plot_grid(title, R_plot, G_plot, B_plot, ncol = 1, nrow = 3)
  RGB_plot <- RGB_plot + geom_stars(data = RGB_stars) +
    scale_fill_gradientn(colors = col_hex, na.value = "transparent") +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text("Pixel coordinates in the first 3 axis \n with RGB color variations"))
  group_plot <- group_plot + geom_stars(data = stars_pixels_group) + 
    scale_fill_gradientn(colors = col_mean_group, na.value = "transparent") +
    coord_fixed() + 
    theme_bw() +
    theme(title = element_text("Group of pixels \n filled by mean color of each group"))
  if (display_plot){
    R_G_B_plot
    RGB_plot
    group_plot
  }
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    save_plot(filename = here("plot", "R_G_B.png"), plot = R_G_B_plot, nrow = 3)
    ggsave(filename = here("plot", "RGB.png"), plot = RGB_plot)
    ggsave(here("plot", "color_group.png"), plot = group_plot)
  }
}

max_min_species_group <- function(list_theta_path, stars_pixels_group, nb_top_species = 20, save_csv = TRUE){
  #' Get species with estimated probabilities of presence highest for each group. Same thing with lowest probabilities.
  #' 
  #' @param list_theta_path character vector. full path to Tiff files containing species probabilities of presence to considers.
  #' @param stars_pixel_group stars object. each pixel contains values of its group number.
  #' @param nb_top_species int. Number of species to display for each group, for lists of lowest and highest probabilities, default is 20.
  #' @param save_csv boolean. Allows to save return in output folder in .csv file, default is TRUE.
  #' @return dataframe. For each group, species name with highest probabilities of presence, mean probabilities, sd probabilities. Same for species with lowest probabilities of presence.
  
  stars_pixels_group <- rast(stars_pixels_group)
  matrix_values <- values(stars_pixels_group)
  names(matrix_values) <- "group"
  for (i in list_theta_path) { 
    # create matrix with probabilities for all species on each pixels
    theta_terra <- rast(i)
    matrix_values <- cbind(matrix_values, values(theta_terra))
    names(matrix_values) <- c(names(matrix_values), names(theta_terra))
  }
  nb_group <- length(unique(matrix_values[, 1]))
  max_species_group <- matrix(0, ncol = nb_group * 3, nrow = nb_top_species) 
  min_species_group <- matrix(0, ncol = nb_group * 3, nrow = nb_top_species) 
  name_min_max <- NULL
  for (group in 1:nb_group) {
    max_species_group[, 3 * group - 2] <- names(sort(colMeans(matrix_values[matrix_values[, 1] == group, ]), decreasing = TRUE))[1:nb_top_species]
    min_species_group[, 3 * group - 2] <- names(sort(colMeans(matrix_values[matrix_values[, 1] == group, ]), decreasing = FALSE))[1:nb_top_species]
    max_species_group[, 3 * group - 1] <- colMeans(matrix_values[matrix_values[, 1] == group, max_species_group[, 3 * group - 2]])
    min_species_group[, 3 * group - 1] <- colMeans(matrix_values[matrix_values[, 1] == group, min_species_group[, 3 * group - 2]])
    max_species_group[, 3 * group] <- colSds(matrix_values[matrix_values[, 1] == group, max_species_group[, 3 * group - 2]])
    min_species_group[, 3 * group] <- colSds(matrix_values[matrix_values[, 1] == group, min_species_group[, 3 * group - 2]])
    name_min_max <- c(name_min_max, paste0(c("Group_", "Mean_", "Sd_"), group))
  }
  colnames(max_species_group) <- paste("Max_", name_min_max)
  colnames(min_species_group) <- paste("Min_", name_min_max)
  min_max_species_group <- data.frame(max_species_group, min_species_group)
  if (save_csv){
    dir.create(here("output"), showWarnings = FALSE)
    write_csv(min_max_species_group, col_names = TRUE, file = here("output", "min_max_species_group.csv"))
  }
  return(min_max_species_group)
}
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
  #' Process Sensitivity (species rightly predicted as present) on inventory sites.
  #'
  #' @param PA dataframe. containing only 0 and 1 with colnames with species names.
  #' @param jSDM_bino_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @return float vector. percentage of species rightly predicted as present on each inventory site.
  
  theta <- JSDM_bino_pro$theta_latent
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
  #' Process Specificity (species rightly predicted as absent) on inventory sites.
  #'
  #' @param PA dataframe. containing only 0 and 1 with colnames with species names.
  #' @param jSDM_bino_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @return float vector. percentage of species rightly predicted as absent on each inventory site.
  
  theta <- JSDM_bino_pro$theta_latent 
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
