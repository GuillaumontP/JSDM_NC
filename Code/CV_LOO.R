library(here)
library(stars)
library(sf)
library(sp)
library(ggplot2)
library(glue)
library(terra)
library("rnaturalearth") # plotting maps
library("rnaturalearthdata")
library("rnaturalearthhires")
library(rgrass7)
library(jSDM)
library(RColorBrewer)

dir.create(here("output", "CV"))
data <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
data$X <- NULL
border <- read_sf(here("output", "moist_forest.shp"))

# Create groups for Cross Validation Leave One (Group) Out
dist_site <- dist(data[,2:3], method = "canberra")
plot(hclust(dist_site), labels = FALSE)
nb_class <- 20
PCA_site_group <- cutree(hclust(dist_site), k = nb_class)
sort(table(PCA_site_group))
data$group <- PCA_site_group

# create color for each group  
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:nb_class]

# plot groups 
NC <- ne_countries(scale = 10, returnclass = "sf")
ggplot(data = NC) +
  geom_sf() +
  geom_point(data = data, aes(x = as.numeric(longitude), y = as.numeric(latitude),color = factor(group)), size = 1) +
  scale_color_manual(values = col_vector) +
  ggtitle("Groups of sites for CV") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
ggsave(here("output", "plot", "group_for_CV_LOO.png"))
# import files for jSDM and CV

set.seed(1234)
EPSG <- 3163
nodat <- -9999
ISO_country_code <- "NCL"
proj.t <- paste0("EPSG:", EPSG) 
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
data_site <- read.csv2(here("output", "data_site.csv"), sep =",")
data_site$X <- NULL
data_site <- cbind(data_site, data)

sens <- rep(0, nb_class)
spe <- sens
TSS <- sens
# CV with prediction only on half class
# for reduce spatial autocorrelation but keeping good extrapolation of alpha and W_i
for (i in 10:nb_class) {
  var_jSDM <- data.matrix(data_site[data$group != i, 1:12]) # n - 1 class
  nb_group <- dim(data.matrix(data_site[data$group == i, 1:12]))[1] 
  half_group <- sample(1:nb_group, nb_group / 2)
  # add half class 
  var_jSDM <- rbind(var_jSDM, data.matrix(data_site[data$group == i, 1:12])[half_group, ])
  var_jSDM[,1] <- as.numeric(var_jSDM[,1] == 1)
  for (j in 2:12) {
    var_jSDM[, j] <- scale(var_jSDM[, j])
  }
  PA_group <- PA[data$group != i, ]
  PA_group <- rbind(PA_group, PA[data$group == i, ][half_group,])
  
  jSDM_binom_pro <- jSDM_binomial_probit(
    burnin = 5000,
    mcmc = 10000,
    thin = 10,
    presence_data = data.matrix(PA_group),
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
  save(jSDM_binom_pro, file = here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  load(here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  
  ##====================
  ##
  ## Spatial Interpolation for alpha and latents variables
  ##
  ##====================
  
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
  params_sites <- data.frame(var_jSDM)
  coord <- rbind(data_site[data_site$group != i, 13:15], data_site[data_site$group == i, 13:15][half_group, ])
  params_sites$latitude <- as.numeric(coord$latitude)
  params_sites$longitude <- as.numeric(coord$longitude)
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
             output={here("output", "CV", paste0("alphas_", i, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
  
  # Representation
  alpha_in <- read_stars(here("output", "CV", paste0("alphas_", i, ".tif")))
  plot(st_crop(alpha_in, border), main = "Site effect alpha interpolated by RST")
  alpha_xy <- rep(0, dim(coord)[1])
  for (j in 1:dim(coord)[1]) 
  {
    alpha_xy[j] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("alphas_", i, ".tif"))} \\
                                             -wgs84 {coord[j,"longitude"]} {coord[j,"latitude"]}'), intern = TRUE))
  }
  plot(alpha_xy, colMeans(jSDM_binom_pro$mcmc.alpha), xlab = "alpha interpolated by RST", ylab = "alpha estimated by JSDM",
       main = "Random site effect")
  abline(a = 0, b = 1, col = 'red')
  
  #====
  # W1 
  #====
  
  W1_sp <- terra::vect(x = data.frame(W1 = colMeans(jSDM_binom_pro$mcmc.latent$lv_1),
                                      x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                       crs = xy@proj4string@projargs)
  rgrass7::write_VECT(W1_sp, "W1")
  
  # Re-sample with RST
  system("v.surf.rst --overwrite --verbose -t tension=3 input=W1 zcolumn=W1 \\
       smooth=0.0 elevation=W1_rst ")
  # Export
  system(glue('r.out.gdal --overwrite input=W1_rst \\
             output={here("output", "CV", paste0("lv_W1_", i, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
  # Representation
  W1_in <- read_stars(here("output", "CV", paste0("lv_W1_", i, ".tif")))
  plot(st_crop(W1_in, border), main = "Latent variable W1 interpolated by RST")
  W1_xy <- rep(0, dim(coord)[1])
  for (j in 1:dim(coord)[1]) 
  {
    W1_xy[j] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("lv_W1_", i, ".tif"))} \\
                                             -wgs84 {coord[j,"longitude"]} {coord[j,"latitude"]}'), intern = TRUE))
  }
  plot(W1_xy, colMeans(jSDM_binom_pro$mcmc.latent$lv_1), xlab = "W1 interpolated by RST", ylab = "W1 estimated by JSDM", 
       main = "First latent axe")
  abline(a = 0, b = 1, col = 'red')
  
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
             output={here("output", "CV", paste0("lv_W2_", i, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
  # Representation
  W2_in <- read_stars(here("output", "CV", paste0("lv_W2_", i, ".tif")))
  plot(st_crop(W2_in, border), main = "Latent variable W2 interpolated by RST")
  W2_xy <- rep(0, dim(coord)[1])
  for (j in 1:dim(coord)[1]) 
  {
    W2_xy[j] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("lv_W2_", i, ".tif"))} \\
                                             -wgs84 {coord[j,"longitude"]} {coord[j,"latitude"]}'), intern = TRUE))
  }
  plot(W2_xy, colMeans(jSDM_binom_pro$mcmc.latent$lv_2), xlab = "W2 interpolated by RST", ylab = "W2 estimated by JSDM", 
       main = "Second latent axe")
  abline(a = 0, b = 1, col = 'red')
  
  # Lambda & Beta
  n_latent <- 2
  np <- 11
  nsp <- length(colnames(jSDM_binom_pro$theta_latent))
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
  
  var_pred <- data_site[data_site$group == i, ]
  result <- matrix(0, ncol = nsp, nrow = dim(var_pred)[1])
  for (k in 1:dim(var_pred)[1]) {
    result[k,] <- pnorm(betas[, 2:13] %*% as.numeric(var_pred[k, 1:12]) + alpha_xy[k] + betas[k, 1] + W1_xy[k] * lambdas[,1] + W2_xy[k] * lambdas[,2])
  }
  
  Sensitivity_CV <- function(PA, theta){
    n_sites <- nrow(theta)
    score <- rep(0, n_sites)
    for(i in 1:n_sites){
      # Sensitivity 
      obs_sp <- which(PA[i,] > 0)
      nobs_sp <- length(obs_sp)
      pred_sp <- which(theta[i,] >= sort(theta[i,], decreasing = TRUE)[nobs_sp])
      score[i] <- sum(pred_sp %in% obs_sp) / ifelse(nobs_sp != 0, nobs_sp, 1)
    }
    return(score)
  }
  Specificity_CV <- function(PA, theta){
    n_sites <- nrow(theta)
    score <- rep(0, n_sites)
    for(i in 1:n_sites){
      # Specificity 
      abs_sp <- which(PA[i,] == 0)
      nabs_sp <- length(abs_sp)
      pred_abs_sp <- which(theta[i,] <= sort(theta[i,])[nabs_sp])
      score[i] <- sum(pred_abs_sp %in% abs_sp) / ifelse(nabs_sp !=0, nabs_sp, 1)
    }
    return(score)
  }
  sens[i] = sum(Sensitivity_CV(PA, result)) / nrow(result)
  spe[i] = sum(Specificity_CV(PA, result)) / nrow(result)
  TSS[i] = sens[i] + spe[i] - 1
  Sys.time()
}
save(sens, spe, TSS, file = here("output", "RData", "TSS_CV.RData"))
Sys.time()

##===================
##
## Analisys of each jSDM model 
##
##===================

for (i in 1:nb_class){
  # readline()
  load(here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  
  top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[5])
  np <- nrow(jSDM_binom_pro$model_spec$beta_start)
  
  # ## beta_j of the top five species
  # par(mfrow = c(3, 2))
  # for (j in top_species) {
  #   for (p in 1:np) {
  #     coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
  #     coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]), 
  #                    main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p],
  #                                 ", species : ", names(top_species[top_species == j])))
  #   }
  # }
  
  ## lambda_j of the top five species
  n_latent <- jSDM_binom_pro$model_spec$n_latent
  par(mfrow = c(2, 2))
  for (j in top_species) {
    for (l in 1:n_latent) {
      coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
      coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]), 
                     main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])
                                  [np + l],", species : ", names(top_species[top_species == j])))
    }
  }
  
  ## Latent variables W_i for the first two sites
  par(mfrow = c(2, 2))
  for (l in 1:n_latent) {
    for (i in 1:2) {
      coda::traceplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                      main = paste0("Latent variable W_", l, ", site ", i))
      coda::densplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                     main = paste0("Latent variable W_", l, ", site ", i))
    }
  }
  
  ## alpha_i of the first two sites
  plot(coda::as.mcmc(jSDM_binom_pro$mcmc.alpha[, 1:2]))
  
  ## V_alpha
  par(mfrow = c(2, 2))
  coda::traceplot(jSDM_binom_pro$mcmc.V_alpha)
  coda::densplot(jSDM_binom_pro$mcmc.V_alpha)
  ## Deviance
  coda::traceplot(jSDM_binom_pro$mcmc.Deviance)
  coda::densplot(jSDM_binom_pro$mcmc.Deviance)
  
  ## probit_theta
  par (mfrow = c(2, 1))
  hist(jSDM_binom_pro$probit_theta_latent,
       main = "Predicted probit theta", xlab = "predicted probit theta")
  hist(jSDM_binom_pro$theta_latent,
       main = "Predicted theta", xlab = "predicted theta")
  
}


##===================
##
## difference between models
##
##===================

# alpha <- read_stars(here("output", "alphas.tif"))
# W1 <- read_stars(here("output", "lv_W1.tif"))
# W2 <- read_stars(here("output", "lv_W2.tif"))
coord <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")[,3:4]
alpha_coord <- rep(0, dim(coord)[1]) 
W1_coord <- alpha_coord
W2_coord <- alpha_coord
for (i in 1:length(coord[,1])) {
  alpha_coord[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "alphas.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
  W1_coord[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W1.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
  W2_coord[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W2.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
Sys.time()
alpha_mat <- matrix(0, ncol = dim(coord)[1], nrow = 20)
W1_mat <- alpha_mat
W2_mat <- alpha_mat
for (j in 1:20) {
  alpha_coord_cv <- rep(0, dim(coord)[1]) 
  W1_coord_cv <- alpha_coord_cv
  W2_coord_cv <- alpha_coord_cv
  for (i in 1:length(coord[,1])) {
    alpha_coord_cv[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("alphas_", j, ".tif"))} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
    W1_coord_cv[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("lv_W1_", j, ".tif"))} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
    W2_coord_cv[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("lv_W2_", j, ".tif"))} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
  }
  alpha_mat[j, ] <- alpha_coord - alpha_coord_cv
  W1_mat[j, ] <- W1_coord - W1_coord_cv
  W2_mat[j, ] <- W2_coord - W2_coord_cv
}
Sys.time()
##============
##
## plot diff btw all dataset & CV
##
##============

matplot(t(alpha_mat), type = "p", col = "black", pch = 3, main = "alpha ds Vs alpha CV")
abline(a = sd(alpha_coord), b = 0, col = "red")
abline(a = -sd(alpha_coord), b = 0, col = "red")
mean(abs(alpha_mat) > sd(alpha_coord))

matplot(t(W1_mat), type = "p", col = "black", pch = 3, main = "W1 ds Vs alpha CV")
abline(a = sd(W1_coord), b = 0, col = "red")
abline(a = -sd(W1_coord), b = 0, col = "red")
mean(abs(W1_mat) > sd(W1_coord))

matplot(t(W2_mat), type = "p", col = "black", pch = 3, main = "W2 ds Vs alpha CV")
abline(a = sd(W2_coord), b = 0, col = "red")
abline(a = -sd(W2_coord), b = 0, col = "red")
mean(abs(W2_mat) > sd(W2_coord))

## 
ggplot() +
  geom_point(aes(x = 1:length(colMeans(alpha_mat)), y = colMeans(alpha_mat))) +
  geom_line(aes(x = 1:length(colMeans(alpha_mat)), y = alpha_coord, col = "all database"), alpha = 0.8) + 
  labs(title = "Alpha")

ggplot() +
  geom_point(aes(x = 1:length(colMeans(W1_mat)), y = colMeans(W1_mat))) +
  geom_line(aes(x = 1:length(colMeans(W1_mat)), y = W1_coord, col = "all database"), alpha = 0.8) +
  labs(title = "W1")

ggplot() +
  geom_point(aes(x = 1:length(colMeans(W2_mat)), y = colMeans(W2_mat))) +
  geom_line(aes(x = 1:length(colMeans(W2_mat)), y = W2_coord, col = "all database"), alpha = 0.8) +
  labs(title = "W2")

beta_check <- matrix(0, ncol = dim(jSDM_binom_pro$model_spec$site_data)[2] + 1 +
                       jSDM_binom_pro$model_spec$n_latent, nrow = dim(PA)[2])
beta_reel <- beta_check

for (i in 1:nb_class) {
  load(here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  for (j in 1:dim(PA)[2]){
    beta_check[j, ] <- beta_check[j, ] + colMeans(jSDM_binom_pro$mcmc.sp[[j]])
  }
}

beta_check <- beta_check / nb_class
load(here("output", "RData", "jSDM_binom_pro.RData"))
for (j in 1:dim(PA)[2]){
  beta_reel[j, ] <- beta_reel[j, ] + colMeans(jSDM_binom_pro$mcmc.sp[[j]])
}
beta_check <- data.frame(beta_check)
beta_reel <- data.frame(beta_reel)
colnames(beta_check) <- c("beta_(Intercept)", "beta_ultramafic", "beta_bio1", "beta_bio4", "beta_bio12", "beta_bio15", "beta_cwd", "beta_bio1^2",
                       "beta_bio4^2", "beta_bio12^2", "beta_bio15^2", "beta_cwd^2", "beta_log(aire)", "lambda_1", "lambda_2")
colnames(beta_reel) <- colnames(beta_check)
for (k in 1:13) {
  beta_check_plot <- beta_check[[colnames(beta_check)[k]]]
  beta_reel_plot <- beta_reel[[colnames(beta_reel)[k]]]
  gg <- ggplot() +
    geom_point(aes(x = 1:dim(PA)[2], y = beta_check_plot[order(beta_reel_plot)], color = "modèle estimé en validation croisée"), alpha = 0.8) +
    geom_point(aes(x = 1:dim(PA)[2], y = sort(beta_reel_plot), color = "modèle avec tout les sites"), alpha = 0.8) +
    theme_bw() +
    ggtitle(paste0("Comparaison of ", colnames(beta_check)[k], "\n between reel model and CV")) +
    labs(y = "valeur de beta", x = "espèces", color = "") +
    scale_color_manual(values = c("red", "black")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  print(gg)
  ggsave(here("output", "plot", paste0("comparaison_CV_all_data_", colnames(beta_check)[k], ".png")))
}
