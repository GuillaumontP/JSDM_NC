create_sf_forest <- function(EPSG = 3163, altitude_min = 0, percentage_min = 50, elevation, forest){
  #' Create sf object of forest with minimal altitude and percentage of forest
  #'
  #' @param EPSG int. target EPSG, default 3163
  #' @param altitude_min float. minimal altitude to consider forest, default 0
  #' @param percentage_min int. percentage minimal of forest to consider pixels as forest, default 50
  #' @param elevation terra object. output of "rast" of "terra" library with elevation in meters
  #' @param forest terra object. output of "rast" of "terra" library with percentage of forest
  #' @return sf object. with only forest higher than `altitude_min` & `percentage_min`
  #'
  #' @import stars
  #' @import sf
  #' @export

  elevation <- reproject_sf_stars(EPSG = EPSG, terra_object = elevation)
  values(elevation)[values(elevation) < 10] <- NA
  forest <- reproject_sf_stars(EPSG = EPSG, terra_object = forest)
  forest[forest < 50] <- NA
  forest <- st_crop(st_as_stars(forest), st_as_stars(elevation))
  forest <- st_as_sf(forest)

  return(forest)
}
