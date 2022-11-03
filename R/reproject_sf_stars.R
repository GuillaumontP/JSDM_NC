reproject_sf_stars <- function(EPSG = 3163, sf_stars_object = NULL, terra_object = NULL){
  #' Reproject shapefile, stars or terra object with right EPSG
  #'
  #' @description
  #' Only one variable between, `sf_stars_object` and `terra_object`.
  #'
  #' @param EPSG target EPSG, default 3163
  #' @param sf_stars_object sf or stars object to reproject
  #' @param terra_object terra object to reproject
  #' @return return same type then input with `EPSG` projection
  #' @import stars
  #' @import sf
  #' @import terra
  #' @export

  if(is.null(sf_stars_object) == is.null(terra_object)){
    print("Remove or add one item bewteen sf_stars_object and terra_object")
    stop()
  }
  if(!is.null(sf_stars_object)){
    return_object <- st_transform(sf_stars_object, EPSG)
  }else{
    return_object <- project(terra_object, paste0("epsg:", EPSG))
  }

  return(return_object)
}
