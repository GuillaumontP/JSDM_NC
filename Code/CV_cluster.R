library(here)
library(EMCluster)
library(combinat)

##===== Initialization
set.seed(1234)
nb_class <- 6
load(here("output", "RData", "tSNE.RData"))
nb_cv <- 1000
nb_pixel <- dim(tSNE_df)[1]

##======
##
## CV for CAH, EM and KM methods
##
##======

# CAH
dist_cah <- dist(tSNE_df)
cah_group_all <- cutree(hclust(dist_cah), k = nb_class)
cah_group_mean <- matrix(0, nrow = nb_class, ncol = dim(tSNE_df)[2])
best_groups_cah_cv <- rep(0, nb_cv)

# EM
init_em <- rand.EM(tSNE_df, nclass = nb_class, min.n = 10)
em_class <- assign.class(tSNE_df, emcluster(tSNE_df, init_em))
em_group_all <- em_class$class
em_group_mean <- matrix(0, nrow = nb_class, ncol = dim(tSNE_df)[2])
best_groups_em_cv <- rep(0, nb_cv)

# KM 
km_class <- kmeans(tSNE_df, centers = nb_class, iter.max = 20, nstart = 1000)
km_group_all <- km_class$cluster
km_group_mean <- matrix(0, nrow = nb_class, ncol = dim(tSNE_df)[2])
best_groups_km_cv <- rep(0, nb_cv)
for (group in 1:nb_class) {
  # centroide of each group on all dataset
  cah_group_mean[group,] <- colMeans(tSNE_df[cah_group_all == group,])
  em_group_mean[group,] <- colMeans(tSNE_df[em_group_all == group,])
  km_group_mean[group,] <- colMeans(tSNE_df[km_group_all == group,])
}

cah_mean <- matrix(0, ncol = dim(tSNE_df)[2], nrow = nb_class)
em_mean <- matrix(0, ncol = dim(tSNE_df)[2], nrow = nb_class)
km_mean <- matrix(0, ncol = dim(tSNE_df)[2], nrow = nb_class)

Sys.time()
for (i in 1:nb_cv) {
  pixels <- sample(1:dim(tSNE_df)[1], size = nb_pixel, replace = TRUE)
  tSNE_cv <- tSNE_df[pixels, ]   
  # CAH
  
  cah_group_cv <- cutree(hclust(dist(tSNE_cv)), k = nb_class)
  # EM
  init_em <- rand.EM(tSNE_cv, nclass = nb_class, min.n = 10)
  em_class_cv <- assign.class(tSNE_cv, emcluster(tSNE_cv, init_em))
  em_group_cv <- em_class_cv$class
  # KM
  km_class <- kmeans(tSNE_cv, centers = nb_class, iter.max = 20, nstart = 1000)
  km_group_cv <- km_class$cluster
  for (group in 1:nb_class) {
    # CAH
    cah_mean[group, ] <- colMeans(tSNE_cv[cah_group_cv == group, ])
    # EM
    em_mean[group, ] <- colMeans(tSNE_cv[em_group_cv == group, ])
    # KM
    km_mean[group, ] <- colMeans(tSNE_cv[km_group_cv == group, ])
  }
  # get dist between group with all data and CV groups
  # CAH
  dist_cah_cv <- as.matrix(dist(rbind(cah_group_mean, cah_mean)))[1:nb_class, nb_class + 1:nb_class]
  sum_dist_cah_cv <- rep(0, length(permn(nb_class)))
  # EM
  dist_em_cv <- as.matrix(dist(rbind(em_group_mean, em_mean)))[1:nb_class, nb_class + 1:nb_class]
  sum_dist_em_cv <- rep(0, length(permn(nb_class)))
  # KM
  dist_km_cv <- as.matrix(dist(rbind(km_group_mean, km_mean)))[1:nb_class, nb_class + 1:nb_class]
  sum_dist_km_cv <- rep(0, length(permn(nb_class)))
  k = 1
  for (perm in permn(nb_class)) {
    sum_dist_cah_cv[k] <- sum(dist_cah_cv[cbind(1:nb_class, unlist(perm))]) # CAH
    sum_dist_em_cv[k] <- sum(dist_em_cv[cbind(1:nb_class, unlist(perm))]) # EM
    sum_dist_km_cv[k] <- sum(dist_km_cv[cbind(1:nb_class, unlist(perm))]) # KM
    k = k + 1   
  }
  best_groups_cah_cv[i] <- min(sum_dist_cah_cv) # CAH
  best_groups_em_cv[i] <- min(sum_dist_em_cv) # EM
  best_groups_km_cv[i] <- min(sum_dist_km_cv) # KM
}
mean(best_groups_cah_cv)
mean(best_groups_em_cv)
mean(best_groups_km_cv)
sd(best_groups_cah_cv)
sd(best_groups_em_cv)
sd(best_groups_km_cv)
Sys.time()
save(best_groups_cah_cv, best_groups_em_cv, best_groups_km_cv, file = here("output", "RData", "best_group_cah_em_km.RData"))

##====== 
##
## Expectation Maximisation 
##
##======

init_em <- rand.EM(tSNE_df, nclass = nb_class, min.n = 10)
em_class <- assign.class(tSNE_df, emcluster(tSNE_df, init_em))
em_group_all <- em_class$class
em_group_mean <- matrix(0, nrow = nb_class, ncol = dim(tSNE_df)[2])
best_groups_em_cv <- rep(0, nb_cv)
for (group in 1:nb_class) {
  # centroide of each group on all dataset
  em_group_mean[group,] <- colMeans(tSNE_df[em_group_all == group,])
}

em_mean <- matrix(0, ncol = dim(tSNE_df)[2], nrow = nb_class)
Sys.time()
for (i in 1:nb_cv) {
  # sample pixels 
  pixels <- sample(1:dim(tSNE_df)[1], size = nb_pixel, replace = TRUE)
  tSNE_em <- tSNE_df[pixels, ]  
  init_em <- rand.EM(tSNE_em, nclass = nb_class, min.n = 10)
  em_class_cv <- assign.class(tSNE_em, emcluster(tSNE_em, init_em))
  em_group_cv <- em_class_cv$class
  for (group in 1:nb_class) {
    em_mean[group, ] <- colMeans(tSNE_em[em_group_cv == group, ])
  }
  # get dist between group with all data and CV groups
  dist_em_cv <- as.matrix(dist(rbind(em_group_mean, em_mean)))[1:nb_class, nb_class + 1:nb_class]
  sum_dist_em_cv <- rep(0, length(permn(nb_class)))
  k = 1
  for (perm in permn(nb_class)) {
    sum_dist_em_cv[k] <- sum(dist_em_cv[cbind(1:nb_class, unlist(perm))])
    k = k + 1   
  }
  best_groups_em_cv[i] <- min(sum_dist_em_cv)
}
mean(best_groups_em_cv)
sd(best_groups_em_cv)
Sys.time()

##==========
##
## K-Means
##
##==========

km_class <- kmeans(tSNE_df, centers = nb_class, iter.max = 20, nstart = 1000)
km_group_all <- km_class$cluster
km_group_mean <- matrix(0, nrow = nb_class, ncol = dim(tSNE_df)[2])
best_groups_km_cv <- rep(0, nb_cv)
for (group in 1:nb_class) {
  # centroide of each group on all dataset
  km_group_mean[group,] <- colMeans(tSNE_df[km_group_all == group,])
}

km_mean <- matrix(0, ncol = dim(tSNE_df)[2], nrow = nb_class)
Sys.time()
for (i in 1:nb_cv) {
  # sample pixels 
  pixels <- sample(1:dim(tSNE_df)[1], size = nb_pixel, replace = TRUE)
  tSNE_km <- tSNE_df[pixels, ]  
  km_class <- kmeans(tSNE_km, centers = nb_class, iter.max = 20, nstart = 1000)
  km_group_cv <- km_class$cluster
  for (group in 1:nb_class) {
    km_mean[group, ] <- colMeans(tSNE_km[km_group_cv == group, ])
  }
  # get dist between group with all data and CV groups
  dist_km_cv <- as.matrix(dist(rbind(km_group_mean, km_mean)))[1:nb_class, nb_class + 1:nb_class]
  sum_dist_km_cv <- rep(0, length(permn(nb_class)))
  k = 1
  for (perm in permn(nb_class)) {
    sum_dist_km_cv[k] <- sum(dist_km_cv[cbind(1:nb_class, unlist(perm))])
    k = k + 1   
  }
  best_groups_km_cv[i] <- min(sum_dist_km_cv)
}
mean(best_groups_km_cv)
sd(best_groups_km_cv)
Sys.time()