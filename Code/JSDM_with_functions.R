library(ade4)
library(cowplot) # multiple ggplot ie par(mfrow)
library(EMCluster)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(ggtern)
library(glue)
library(grid)
library(here)
library(jSDM)
library(matrixStats)
library(RColorBrewer)
library(readr)
library(rgdal)
library(rgl) # plot 3d
library(rgrass7)
library(rnaturalearth) # plotting maps
library(rnaturalearthdata)
library(rnaturalearthhires)
library(Rtsne) # alternative to PCA
library(sp)
library(stars)
library(stringr)
library(terra)
library(terrainr)
library(tidyterra)
library(viridis)

source(here("Code", "JSDM_NC_function.R"))

##=================
##
## Init
##
##=================

set.seed(1234)
EPSG <- 3163

GT <- read_sf(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))

GT <- reproject_sf_stars(EPSG = EPSG, sf_stars_object = GT)
write_sf(GT, here("data_raw", "Grande_Terre", "Grande_Terre.shp"))

write_stars(c(st_normalize(st_crop(read_stars(here("output", "environ_allNC.tif")), GT)),
              st_normalize(st_crop(read_stars(here("output", "current_chelsaNC.tif")), GT)),
              along = "band"),
            here("output", "dataNC.tif"),
            options = c("COMPRESS=LZW","PREDICTOR=2"))

elevation <- split(read_stars(here("output", "dataNC.tif")))["elevation"]
forest_per <- split(read_stars(here("output", "dataNC.tif")))["forest"]

##=================
##
## Create moist forest mask
##
##=================

forest <- create_sf_forest(EPSG = EPSG, altitude_min = 10, percentage_min = 50, elevation = elevation, forest = forest_per)
monthly_precipitation <- split(read_stars(here("output", "current_chelsaNC.tif")))[paste0("pr", 1:12)]

moist_forest_sf <- create_sf_moist_forest(EPSG = EPSG, forest = forest, monthly_precipitation = monthly_precipitation)
write_sf(moist_forest_sf, here("output", "moist_forest.shp"))

##=================
##
## Create explanorities variables in stars files
##
##=================

data_NC <- read_stars(here("output", "dataNC.tif"))
var_jSDM <- extract_var_JSDM(stars_object = data_NC, variables_names = c("peridotites", "bio1", "bio4", "bio12", "bio15", "cwd"),
                 power_variable = c(1, rep(2,5)), scale = c(FALSE, rep(TRUE, 5)))
write_stars(merge(var_jSDM), here("output", "var_jSDM.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))

coord_site <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
coord_site[, "AREA"] <- 0
NC_PIPPN <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep = ",")

for (site in 1:dim(coord_site)[1]) {
  coord_site$AREA[site] <- as.numeric(NC_PIPPN$AREA[coord_site$plot_name[site] == NC_PIPPN$plot_name][1])
}
var_site <- data_JSDM(EPSG = EPSG, latlon_site = coord_site[, 2:3], area_site = coord_site$AREA, 
                      path_tiff = here("output", "var_jSDM.tif"), log_area_site = TRUE)


##=================
##
## jSDM binomial probit
##
##=================

PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
jSDM_bin_pro <- JSDM_bino_pro(presence_data = PA, site_data = var_site, n_latent = 2, V_beta = c(rep(1, 12), 0.001),
              mu_beta = c(rep(0, 12), 0.25), nb_species_plot = 2, display = TRUE, save_plot = TRUE)

save(jSDM_bin_pro, here("jSDM_bin_pro.RData"))

plot_pres_est_one_species(jSDM_binom_pro = jSDM_bin_pro, coord_site = coord_site, country_name = "New Caledonia",
                          display_plot = TRUE, save_plot = TRUE)

plot_species_richness(jSDM_binom_pro = jSDM_bin_pro, coord_site = coord_site, country_name = "New Caledonia",
                      display_plot = TRUE, save_plot = TRUE)
##=================
##
## Knn interpolation
##
##=================

knn_interpolation_jSDM(jSDM_binom_pro = jSDM_bin_pro, k = 5, coord_site = coord_site, sf_forest = moist_forest_sf)
# 25min

##=================
##
## Probabilities on forest area and plots
##
##=================

alpha <- read_stars(here("output", "alpha_knn.tif"))
W1 <- read_stars(here("output", "lv_1_knn.tif"))
W2 <- read_stars(here("output", "lv_2_knn.tif"))
W <- c(W1, W2, along = "band")
sc <- scale(log(coord_site$AREA))
log_area_scale <- (log(1000 * 1000) - attr(sc, "scaled:center")) / attr(sc, "scaled:scale")
area_stars <- alpha
area_stars[[1]][which(!is.na(area_stars[[1]]), arr.ind = TRUE)] <- log_area_scale / 2 ###############
names(area_stars) <- "log_area"

prob_est_species_forest(alpha_stars = alpha, latent_var_stars = W, jSDM_binom_pro = jSDM_bin_pro,
                        data_stars = var_jSDM, area_stars = area_stars)

plot_prob_pres_interp(theta_stars = read_stars(here("output", "theta", "KNN_theta_01.tif")), country_sf = GT, 
                      display_plot = TRUE, save_plot = TRUE)

species_richness_terra <- plot_species_richness_interpolated(list_theta_path = list.files(here("output", "theta"), full.names = TRUE),
                                                             country_sf = GT, display_plot = TRUE, save_plot = TRUE, save_tif = TRUE)

##=================
##
## Compute PCA or tSNE and clustering with KM/EM/CAH
##
##=================
Sys.time()
coord_PCA_or_tSNE <- PCA_or_tSNE_on_pro_pres_est(tSNE = TRUE, PCA = FALSE, list_theta_path = list.files(here("output", "theta"), full.names = TRUE),
                                          display_plot = TRUE, save_plot = TRUE)
Sys.time()

KM_groups <- plot_HCA_EM_Kmeans(method = "KM", pixel_df = coord_PCA_or_tSNE, nb_group = 5, display_plot_3d = TRUE,
                   display_plot_2d = TRUE, save_plot_2d = TRUE)

##=================
##
## Plot with colors and groups
##
##=================

pixels_group <- split(read_stars(here("output", "theta", "KNN_theta_01.tif")))[1,,]
pixels_group[[1]][!is.na(pixels_group[[1]])] <- KM_groups
names(pixels_group) <- "groups"

plot_RGB_group_by_color(stars_pixels_group = pixels_group, coord_pixel = coord_PCA_or_tSNE, country_sf = GT,
                        display_plot = TRUE, save_plot = TRUE)

max_min_df <- max_min_species_group(list_theta_path = list.files(here("output", "theta"), full.names = TRUE), 
                                    stars_pixels_group = pixels_group, nb_top_species = 20, save_csv = TRUE)

##=================
##
## Sensibility & Specificity
##
##=================

Sensitivity(PA = PA, jSDM_binom_pro = jSDM_bin_pro)
Specificity(PA = PA, jSDM_binom_pro = jSDM_bin_pro)
