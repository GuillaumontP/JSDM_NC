reproject_sf_stars <- function(EPSG = 3163, sf_stars_object){
  #' Reproject shapefile with right EPSG
  #'
  #' @param EPSG target EPSG, default 3163
  #' @param sf_stars_object sf or stars object to reproject
  #' @return same object then  `sf_stars_object` with `EPSG` projection
  #' @import stars
  #' @import sf
  #' @export

  sf_stars_object <- st_transform(sf_stars_object, EPSG)

  return(sf_stars_object)
}
