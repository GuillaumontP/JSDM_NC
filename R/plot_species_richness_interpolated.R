plot_species_richness_interpolated <- function(list_theta_path, country_name = NULL,
                                               country_sf = NULL, display_plot = TRUE, save_plot = FALSE, save_tif = FALSE){
  #' Create plot with species richness with optional sf border.
  #'
  #' @param list_theta_path character vector. with full path of Tiff files who contain probabilities of presence of each species to consider.
  #' @param country_name character. optional, for display border of country on map, default is NULL.
  #' @param country_sf sf object. optional, inconsistent with `country_name`. Display borders on map, default is NULL.
  #' @param display_plot boolean. Display plot, default is TRUE.
  #' @param save_plot boolean. Save plot in .png format in folder plot, default is FALSE.
  #' @param save_tif boolean. Save species richness in .tif file in folder output, default is FALSE.
  #' @return terra object. species richness in designated area.
  #'
  #' @import terra
  #' @import rnaturalearth
  #' @import stars
  #' @import sf
  #' @import ggplot2
  #' @import viridis
  #' @export

  theta_sum <- sum(terra::rast(list_theta_path[1]))
  n_files <- length(list_theta_path)
  if (n_files > 1){
    for (i in list_theta_path[2:n_files])
    {
      theta_sum <- theta_sum + sum(terra::rast(i))
    }
  }
  if (save_tif){
    dir.create(here("output"), showWarnings = FALSE)
    terra::writeRaster(theta_sum, here("output", "species_richness.tif"), overwrite = TRUE)
  }
  if(!is.null(country_name)){
    map <- ne_countries(scale = 10, returnclass = "sf", country = country_name)[1]
    map <- st_transform(map, st_crs(theta_sum))
    latlon_output <- st_bbox(map)
  }
  if(!is.null(country_sf)){
    latlon_output <- st_bbox(country_sf)
    map <- country_sf
  }
  theta_sum <- st_as_stars(theta_sum)
  theta_sum[[1]] <- pmin(theta_sum[[1]], 400) ########################################"
  gplot <- ggplot() +
    geom_sf(data = map, colour = "black", fill = "grey") +
    geom_stars(data = theta_sum) +
    ggtitle("Estimated current species richness") +
    scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE) +
    theme_bw() +
    labs(fill = "Number \n of species") +
    theme(plot.title = element_text(hjust = 0.5))
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(here("plot", "estimated_species_richness_forest.png"))
  }
  return(theta_sum)
}
