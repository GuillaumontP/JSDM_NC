data_JSDM <- function(EPSG = 3163, latlon_site, area_site, path_tiff, log_area_site = TRUE){
  #' Create dataframe with values of each variable on each site
  #'
  #' @param EPSG int. target EPSG, default 3163
  #' @param latlon_site float matrix with 2 columns. Columns are longitude and latitude in this order.
  #' @param area_site float vector same dimension then number of row `latlon`. area of each site in m²
  #' @param path_tiff character. path to Tiff file with explanatories variables
  #' @param log_area_site boolean. Replace area of each site by log area of each site. Anyway, this variable is scale.
  #' @return matrix with latitude, longitude and all variables. If `latlon_site` has row names, output has same row name.
  #'
  #' @import stars
  #' @importFrom glue glue
  #' @export

  if (dim(latlon_site)[1] != length(area_site)){
    print("number of rows in latlon is different than lengh of area_site")
    stop()
  }
  stars_object <- split(read_stars(path_tiff))
  if (log_area_site){
    area_site <- scale(log(area_site))
    name_area <- "log_area_site"
  }else{
    area_site <- scale(area_site)
    name_area <- "area_site"
  }
  data_site <- matrix(0, ncol = length(names(stars_object)), nrow = length(area_site))
  row.names(data_site) <- row.names(latlon_site)
  proj.t <- paste0("EPSG:", EPSG)
  for (i in 1:dim(latlon_site)[1]) {
    data_site[i, ] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {path_tiff} \\
                                                                       -wgs84 {latlon_site[i,1]} {latlon_site[i,2]}'), intern = TRUE))
  }
  data_site <- cbind(data_site, area_site)
  colnames(data_site) <- c(names(stars_object), name_area)

  return(data_site)
}
