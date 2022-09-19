library(stats)
library(here)
library(stars)
library(terra)
library(stringr)
library(EMCluster)
library(factoextra)
library(terrainr)
library(grid)
library(combinat) # use for CV of methods
# permn(1:6)

set.seed(1234)
dir.create(here("output", "classification"))
##================
##
## Dendogram on tSNE output
## determine number of class 
## CAH
##================

load(here("output", "RData", "tSNE.RData"))
dist_col <- dist(tSNE_df)
plot(hclust(dist_col), labels = FALSE) # 4 or 6 groups
nb_class <- 6
tSNE_group <- cutree(hclust(dist_col), k = nb_class)
for (j in 1:nb_class) {
  group <- which(tSNE_group == j)
  write.csv2(sort(group), here("output", "classification", paste0("group_tSNE_CAH_", j, ".csv")))
}
print(table(tSNE_group))

##================
##
## Expectation Maximisation
##
##================

init_EM <- rand.EM(tSNE_df, nclass = nb_class, min.n = 10)
EM_class <- assign.class(tSNE_df, emcluster(tSNE_df, init_EM))
EM_group <- EM_class$class
for (j in 1:nb_class) {
  group <- which(EM_group == j)
  write.csv2(sort(group), here("output", "classification", paste0("group_tSNE_EM_", j, ".csv")))
}
print(table(EM_group))

##================
##
## K-means
##
##================

KM_class <- kmeans(tSNE_df, centers = nb_class, iter.max = 20, nstart = 1000)
KM_group <- KM_class$cluster
for (j in 1:nb_class) {
  group <- which(KM_group == j)
    write.csv2(sort(group), here("output", "classification", paste0("group_tSNE_KM_", j, ".csv")))
}
print(KM_class$size)

##================
##
## Comparaison between EM, KM and CAH
##
##================

## EM -> KM
similarity_EM_KM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_EM_KM) <- paste0("EM_", 1:nb_class)
rownames(similarity_EM_KM) <- c("number of pixels", "number of pixels in common", " percentage of EM group", "number of KM group")
for (j in 1:nb_class) {
  EM <- read.csv2(here("output", "classification", paste0("group_tSNE_EM_", j, ".csv")))$x
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    KM <- read.csv2(here("output", "classification", paste0("group_tSNE_KM_", k, ".csv")))$x
    nb_common[k] <- as.numeric(max(length(intersect(EM, KM)), 0))
  }
  similarity_EM_KM[1, j] <- length(EM)
  similarity_EM_KM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_EM_KM[3, j] <- as.integer(similarity_EM_KM[2, j] / length(EM) * 100) 
  similarity_EM_KM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_EM_KM)

## KM -> EM
similarity_KM_EM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_KM_EM) <- paste0("KM_", 1:nb_class)
rownames(similarity_KM_EM) <- c("number of pixels", "number of pixels in common", " percentage of KM group", "number of EM group")
for (j in 1:nb_class) {
  KM <- read.csv2(here("output", "classification", paste0("group_tSNE_KM_", j, ".csv")))$x
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    EM <- read.csv2(here("output", "classification", paste0("group_tSNE_EM_", k, ".csv")))$x
    nb_common[k] <- as.numeric(max(length(intersect(EM, KM)), 0))
  }
  similarity_KM_EM[1, j] <- length(KM)
  similarity_KM_EM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_KM_EM[3, j] <- as.integer(similarity_KM_EM[2, j] / length(KM) * 100) 
  similarity_KM_EM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_KM_EM)

## CAH -> EM
similarity_CAH_EM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_CAH_EM) <- paste0("CAH_", 1:nb_class)
rownames(similarity_CAH_EM) <- c("number of pixels", "number of pixels in common", " percentage of CAH group", "number of EM group")
for (j in 1:nb_class) {
  CAH <- read.csv2(here("output", "classification", paste0("group_tSNE_CAH_", j, ".csv")))$x
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    EM <- read.csv2(here("output", "classification", paste0("group_tSNE_EM_", k, ".csv")))$x
    nb_common[k] <- as.numeric(max(length(intersect(EM, CAH)), 0))
  }
  similarity_CAH_EM[1, j] <- length(CAH)
  similarity_CAH_EM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_CAH_EM[3, j] <- as.integer(similarity_CAH_EM[2, j] / length(CAH) * 100) 
  similarity_CAH_EM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_CAH_EM)

## EM -> CAH
similarity_EM_CAH <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_EM_CAH) <- paste0("EM_", 1:nb_class)
rownames(similarity_EM_CAH) <- c("number of pixels", "number of pixels in common", " percentage of EM group", "number of CAH group")
for (j in 1:nb_class) {
  EM <- read.csv2(here("output", "classification", paste0("group_tSNE_EM_", j, ".csv")))$x
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    CAH <- read.csv2(here("output", "classification", paste0("group_tSNE_CAH_", k, ".csv")))$x
    nb_common[k] <- as.numeric(max(length(intersect(EM, CAH)), 0))
  }
  similarity_EM_CAH[1, j] <- length(EM)
  similarity_EM_CAH[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_EM_CAH[3, j] <- as.integer(similarity_EM_CAH[2, j] / length(EM) * 100) 
  similarity_EM_CAH[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_EM_CAH)

## CAH -> KM
similarity_CAH_KM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_CAH_KM) <- paste0("CAH_", 1:nb_class)
rownames(similarity_CAH_KM) <- c("number of pixels", "number of pixels in common", " percentage of CAH group", "number of KM group")
for (j in 1:nb_class) {
  CAH <- read.csv2(here("output", "classification", paste0("group_tSNE_CAH_", j, ".csv")))$x
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    KM <- read.csv2(here("output", "classification", paste0("group_tSNE_KM_", k, ".csv")))$x
    nb_common[k] <- as.numeric(max(length(intersect(KM, CAH)), 0))
  }
  similarity_CAH_KM[1, j] <- length(CAH)
  similarity_CAH_KM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_CAH_KM[3, j] <- as.integer(similarity_CAH_KM[2, j] / length(CAH) * 100) 
  similarity_CAH_KM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_CAH_KM)

## KM -> CAH
similarity_KM_CAH <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_KM_CAH) <- paste0("KM_", 1:nb_class)
rownames(similarity_KM_CAH) <- c("number of pixels", "number of pixels in common", " percentage of KM group", "number of CAH group")
for (j in 1:nb_class) {
  KM <- read.csv2(here("output", "classification", paste0("group_tSNE_KM_", j, ".csv")))$x
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    CAH <- read.csv2(here("output", "classification", paste0("group_tSNE_CAH_", k, ".csv")))$x
    nb_common[k] <- as.numeric(max(length(intersect(KM, CAH)), 0))
  }
  similarity_KM_CAH[1, j] <- length(KM)
  similarity_KM_CAH[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_KM_CAH[3, j] <- as.integer(similarity_KM_CAH[2, j] / length(KM) * 100) 
  similarity_KM_CAH[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_KM_CAH)

##================
##
## get pixels who are always in same group
##
##================
print(similarity_EM_KM)
print(similarity_KM_EM)
print(similarity_CAH_KM)
print(similarity_KM_CAH)
print(similarity_CAH_EM)
print(similarity_EM_CAH)

nb_class <- 5
final_groups <- list(NULL)
nb_pixels_by_group <- rep(0, nb_class)
for (i in 1:5) {
  pca <- read.csv2(here("output", "classification", paste0("group_CAH_", i, ".csv")))
  em_group <- similarity_CAH_EM[4, i]
  km_group <- similarity_CAH_KM[4, i]
  em <- read.csv2(here("output", "classification", paste0("group_EM_", em_group, ".csv")))
  km <- read.csv2(here("output", "classification", paste0("group_KM_", km_group, ".csv")))
  common_pixels <- intersect(intersect(pca, em), km)
  nb_pixels_by_group[i] <- length(common_pixels)
  final_groups[[i]] <- common_pixels
  write.csv2(sort(final_groups[[i]]), here("output", "classification", paste0("final_group_", i, ".csv")))
}
nb_pixels_by_group
sum(nb_pixels_by_group)

##================
##
## Color of each group 
##
##================

npart <- 30
nb_pixel_color <- 600
col_group <- matrix(0, nrow = nb_class, ncol = 3)
hex_color <- rep(0, nb_class)
for (i in 1:nb_class) {
  group_i <- read.csv2(here("output", "classification", paste0("final_group_", i, ".csv")))
  for (j in 1:length(group_i)){
    group_i[j] <- paste0("X", data_clear$id_taxon_ref[data_clear$nom_taxon == group_i[j]][1])
  }
  theta_group <- terra::rast(here("output", "theta", "RST_theta_forest_01.tif"))
  sum_theta_group <- theta_group[["X1644"]] - theta_group[["X1644"]] # init rast empty with right crs, extent,...
  for (k in 1:npart) {
    theta_group <- terra::rast(here("output", "theta", paste0("RST_theta_forest_", str_pad(k, width = 2, pad = "0"), ".tif")))
    common <- intersect(names(theta_group), group_i)
    if (length(common) != 0)
    {
      # stack probabilities of each species of each group
      sum_theta_group <- sum_theta_group + sum(theta_group[[common]])
    }
  }
  terra::writeRaster(sum_theta_group, here("output", "classification", paste0("theta_group_", i, ".tif")), overwrite = TRUE)
  sum_theta_group <- st_as_stars(sum_theta_group)
  # value_min <- sort(sum_theta_group[[1]], decreasing = TRUE)[nb_pixel_color]
  nb_in_group <- length(group_i)
  sum_theta_group[[1]][sum_theta_group[[1]] < 0.6 * nb_in_group] <- NA
  sum_theta_group_sf <- st_as_sf(sum_theta_group)
  write_stars(st_crop(read_stars(here("output", "species_turnover.tif")), sum_theta_group_sf),
              dsn = here("output", "classification", "species_group.tif"))
  
  plotRGB(terra::rast(here("output", "classification", "species_group.tif")), main = paste0("Main area of species from group ", i))
  RGB_map <- read_stars(here("output", "species_turnover.tif"))
  R <- RGB_map[[1]][,,1]
  G <- RGB_map[[1]][,,2]
  B <- RGB_map[[1]][,,3]
  # get values for each band RGB
  R <- R[!is.na(sum_theta_group[[1]])]
  G <- G[!is.na(sum_theta_group[[1]])]
  B <- B[!is.na(sum_theta_group[[1]])]
  
  # Find nearest color of mean color
  mindist <- which.min(sum(colMeans(cbind(R, G, B))^2) - rowSums(cbind(R, G, B)^2))

  R_near <- R[mindist] 
  G_near <- G[mindist]
  B_near <- B[mindist] 
  col_group[i, ] <- c(R_near, G_near, B_near)
  hex_color[i] <- rgb(red = col_group[i,1], green = col_group[i,2], blue = col_group[i,3], maxColorValue = 255)
}

ggplot() +
  geom_spatial_rgb(data = terra::rast(here("output", "species_turnover.tif")), 
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  geom_label(aes(x = 200000, y = 260000, label = "Group 1", size = 5), fill = hex_color[1]) +
  geom_label(aes(x = 200000, y = 245000, label = "Group 2", size = 5), fill = hex_color[2]) +
  geom_label(aes(x = 200000, y = 230000, label = "Group 3", size = 5), fill = hex_color[3]) +
  geom_label(aes(x = 200000, y = 215000, label = "Group 4", size = 5), fill = hex_color[4]) +
  geom_label(aes(x = 200000, y = 200000, label = "Group 5", size = 5), fill = hex_color[5]) +
  theme_bw() +
  theme(legend.position = "none")

# sum(!is.na(read_stars(here("output", "classification", "theta_group_1.tif"))[[1]])) = 6548

##================
##
## CAH on probabilities for each species and each pixels
##
##================

list_theta <- list.files(here("output", "theta"), pattern = "RST_theta_forest", full.names = TRUE)
load(here("output", "cell.RData"))
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
prob_all <- matrix(0, ncol = dim(PA)[2], nrow = dim(cell)[1])
for (j in 1:length(list_theta)) {
  theta <- read_stars(list_theta[j])[[1]]
  for (i in 1:dim(cell)[1]) {
    prob_all[i, (1 + dim(theta)[3] * (j - 1)): min(dim(theta)[3] * j, dim(PA)[2])] <- theta[cell[j, 1], cell[j, 2], ]
  }
}
# CAH isn't working
nb_class <- 5
init_EM <- rand.EM(prob_all, nclass = nb_class, min.n = 10)
EM_class <- assign.class(prob_all, emcluster(prob_all, init_EM))
EM_group <- EM_class$class
names(EM_group) <- rownames(prob_all)
data_clear <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep =",")
for (j in 1:nb_class) {
  group <- names(EM_group)[EM_group == j]
}
print(table(EM_group))
Sys.time()
