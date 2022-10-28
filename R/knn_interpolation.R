knn_interpolation_jSDM <- function(jSDM_binom_pro, k = 5, coord_site, sf_forest){
  #' Interpole alpha and all latent variables with knn algorithm. Then save interpolation on disk in output folder.
  #'
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @param k int. number of neighbours to consider in knn algorithm, default is 5.
  #' @param coord_site dataframe. columns name must contain latitude and longitude as names.
  #' @param sf_forest sf object. area where interpolation will be done.
  #'
  #' @import terra
  #' @import sp
  #' @import stars
  #' @export

  dir.create(here("output"), showWarnings = FALSE)
  coord_pixel_latlon <- coordinates(spTransform(as_Spatial(sf_forest), CRS("+proj=longlat +datum=WGS84")))
  # create stars object with forest's arguments
  stars_object <- matrix(0, nrow = (st_bbox(sf_forest)[3] - st_bbox(sf_forest)[1]) / 1000,
                         ncol = (st_bbox(sf_forest)[4] - st_bbox(sf_forest)[2]) / 1000)
  stars_object <- st_as_stars(stars_object)
  stars_object <- st_set_crs(stars_object, st_crs(sf_forest))
  stars_object <- st_set_bbox(stars_object, st_bbox(sf_forest))
  st_dimensions(stars_object)$X1$delta <- 1000
  st_dimensions(stars_object)$X2$delta <- -1000
  stars_object <- st_set_bbox(stars_object, st_bbox(sf_forest))
  stars_object <- st_crop(stars_object, sf_forest)
  distance <- stars_object
  coord_pixel_df <- which(!is.na(stars_object[[1]]), arr.ind = TRUE)
  for (site in 1:dim(coord_site)[1]) # distance between pixels & each site ~16min
  {
    for (pixel in 1:dim(coord_pixel_latlon)[1])
    {
      distance[[1]][coord_pixel_df[pixel, 1],
                    coord_pixel_df[pixel, 2]] <- sqrt((coord_pixel_latlon[pixel, 1] - as.numeric(coord_site[site, "longitude"]))^2 +
                                                        (coord_pixel_latlon[pixel, 2] - as.numeric(coord_site[site, "latitude"]))^2)
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
    names_coord_site <- colnames(coord_site)
    if ( name_inter == "alpha"){
      coord_site[, ncol(coord_site) + 1] <- colMeans(jSDM_binom_pro$mcmc.alpha)
      all_distance <- merge(all_distance)
    }else{
      coord_site[, ncol(coord_site) + 1] <- colMeans(jSDM_binom_pro$mcmc.latent[[name_inter]])
    }
    colnames(coord_site) <- c(names_coord_site, name_inter)
    for (site in 1:dim(coord_pixel_latlon)[1]) {
      dis_pixel <- all_distance[[1]][coord_pixel_df[site, 1], coord_pixel_df[site, 2], ]
      knn_site_nb <- which(dis_pixel >= sort(dis_pixel, decreasing = FALSE)[k]) # site number
      knn_dis <- dis_pixel[knn_site_nb] # distance to each site
      inter_value <- sum(coord_site[knn_site_nb, ncol(coord_site)] * (1 - knn_dis / sum(knn_dis)) / (length(knn_site_nb) - 1))
      interpolated_stars[[1]][coord_pixel_df[site, 1], coord_pixel_df[site, 2]] <- inter_value
    }
    if (name_inter == "alpha"){
      interpolated_stars[[1]] <- interpolated_stars[[1]] - mean(interpolated_stars[[1]], na.rm = TRUE)
    }else{
      interpolated_stars[[1]] <-  (interpolated_stars[[1]] - mean(interpolated_stars[[1]], na.rm = TRUE)) / sd(interpolated_stars[[1]], na.rm = TRUE)
    }
    write_stars(interpolated_stars, here("output", paste0(name_inter, "_knn.tif")), options = c("COMPRESS=LZW","PREDICTOR=2"))
  }
}
