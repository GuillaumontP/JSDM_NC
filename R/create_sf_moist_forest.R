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
  #' @export

  forest <- reproject_sf_stars(EPSG = EPSG, sf_stars_object = forest)
  monthly_precipitation <- reproject_sf_stars(EPSG = EPSG, sf_stars_object = monthly_precipitation)
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
    moist_month[[1]] <- moist_month[[1]] + as.integer(merge(monthly_precipitation)[[1]][,,i] < 100)
  }
  moist_forest <- moist_month * MAP
  moist_forest[[1]][moist_forest[[1]] > 5] <- NA # remove no moist forest
  moist_forest <- st_crop(moist_forest, forest)
  moist_forest <- st_as_sf(moist_forest)

  return(moist_forest)
}
