plot_prob_pres_interp <- function(theta_stars, species_to_plot = names(theta_stars)[1], country_name = NULL,
                                  country_sf = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Create plot with probabilities of presence interpolated for a choosen species. Plot can be hide and/or save in plot folder
  #'
  #' @param theta_stars multilayer stars_object. from files save in prob_est_species_forest function.
  #' @param species_to_plot character. name of one layer of `theta_stars`, default is first one.
  #' @param country_name character. optional, for display border of country on map, default is NULL.
  #' @param country_sf sf object. optional, inconsistent with `country_name`. Display borders on map, default is NULL.
  #' @param display_plot boolean. Display plot, default is TRUE.
  #' @param save_plot boolean. Save plot in .png format in folder plot, default is FALSE.
  #'
  #' @import rnaturalearth
  #' @import stars
  #' @import ggplot2
  #' @import viridis
  #' @export

  if (!is.null(country_name) & !is.null(country_sf)){
    print("Only country_name or country_sf")
    stop()
  }
  if(!is.null(country_name)){
    map <- ne_countries(scale = 10, returnclass = "sf", country = country_name)[1]
    map <- st_transform(map, st_crs(theta_stars))
    latlon_output <- st_bbox(map)
  }
  if(!is.null(country_sf)){
    latlon_output <- st_bbox(country_sf)
    map <- country_sf
  }
  if (length(theta_stars) == 1){
    if (species_to_plot == names(theta_stars)[1]){
      species_to_plot = names(split(theta_stars))[1]
    }
    theta_stars <- split(theta_stars)
  }
  gplot <- ggplot() +
    geom_sf(data = map, colour = "black", fill = "grey") +
    geom_stars(data = theta_stars[species_to_plot,,]) +
    ggtitle(paste('Interpolated current probabilities of presence \n for', species_to_plot)) +
    scale_fill_gradientn(colours = rev(rocket(9)), na.value = "transparent") +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE) +
    theme_bw() +
    labs(fill = "Number of species") +
    theme(plot.title = element_text(hjust = 0.5))
  gplot
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(plot = gplot, here("plot", paste0("interpolated_presence_for_", species_to_plot, ".png")))
  }
}
