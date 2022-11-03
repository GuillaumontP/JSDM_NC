create_sf_moist_forest <- function(EPSG = 3163, forest, monthly_precipitation){
  #' Create sf object of moist forest according to FAO's criteria.
  #'  MAP >=1500 mm
  #' <5 months where precipitations < 100 mm
  #'
  #' @param EPSG int. target EPSG, default 3163
  #' @param forest sf object. output of "read_sf" of "sf" library with boolean values of forest.
  #' @param monthly_precipitation multilayers stars object. output of "read_stars" of "stars" library with montly mean of precipitation in mm.
  #' @return sf object. with only moist forest area.
  #'
  #' @import stars
  #' @import sf
  #' @import terra
  #' @export

  forest <- reproject_sf_stars(EPSG = EPSG, sf_stars_object = forest)
  monthly_precipitation <- reproject_sf_stars(EPSG = EPSG, terra_object = monthly_precipitation)
  annual_precipitation <- sum(monthly_precipitation)
  MAP <- annual_precipitation < 1500
  moist_month <- monthly_precipitation[[1]]
  values(moist_month) <- 0
  for (i in 1:12) {
    month_100 <- monthly_precipitation[[i]] < 100
    values(month_100) <- as.integer(values(month_100))
    moist_month <- moist_month + month_100
  }
  moist_forest <- moist_month * MAP
  values(moist_forest)[values(moist_forest) > 5] <- NA # remove no moist forest
  values(moist_forest)[values(moist_forest) > 5] <- NA
  moist_forest <- st_crop(st_as_stars(moist_forest), st_as_stars(forest))
  moist_forest <- st_as_sf(moist_forest)

  return(moist_forest)
}
