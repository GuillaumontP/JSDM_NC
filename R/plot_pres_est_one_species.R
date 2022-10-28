plot_pres_est_one_species <- function(jSDM_binom_pro, species_to_plot = colnames(jSDM_binom_pro$model_spec$presence_data)[1],
                                      coord_site, country_name = NULL, latlon_output = NULL, display_plot = TRUE,
                                      save_plot = FALSE){
  #' Plot presence absence  and estimated probabilities of presence for one species. Possibility to save plots
  #'
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @param species_to_plot character. species to plot, default is the first species of the list.
  #' @param coord_site dataframe. columns name must contain latitude and longitude as names.
  #' @param country_name character. English name of the country where are inventory sites, default is NULL.
  #' @param latlon_output float vector. inconsistent with `country_name` coord of output box is this format {lon_min, lat_min, lon_max, lat_max}, default is NULL.
  #' @param display_plot boolean. show plot, default is TRUE.
  #' @param save_plot boolean. Write plot in plot folder as .png files, default is FALSE.
  #'
  #' @import stars
  #' @import rnaturalearth
  #' @import ggplot2
  #' @export

  if (is.null(country_name) & is.null(latlon_output)){
    print("Need latlon coordinates")
    stop()
  }
  if (!is.null(country_name)){
    latlon_output <- st_bbox(ne_countries(scale = 10, country = country_name, returnclass = "sf"))
  }
  map <- ne_countries(scale = 10, returnclass = "sf")
  df_PA <- jSDM_binom_pro$model_spec$presence_data
  names_coord_site <- colnames(coord_site)
  coord_site <- data.frame(cbind(coord_site, df_PA[, species_to_plot]))
  colnames(coord_site) <- c(names_coord_site, species_to_plot)
  pres <- data.frame(coord_site[coord_site[, species_to_plot] == 1, ][, c("longitude", "latitude")])
  abs <- data.frame(coord_site[coord_site[, species_to_plot] == 0, ][, c("longitude", "latitude")])
  obs_plot <- ggplot(data = map) +
    geom_sf() +
    geom_point(data = coord_site[, c("longitude", "latitude", species_to_plot)],
               aes(x = as.numeric(longitude), y = as.numeric(latitude), shape = as.factor(coord_site[, species_to_plot]),
                   color = as.factor(coord_site[, species_to_plot])),
               size = 1) +
    ggtitle(paste0('Observed occurences of \n', species_to_plot)) +
    scale_shape_manual(labels = c("Absence", "Presence"), values = c(1, 3), name = "Species state") +
    scale_color_manual(labels = c("Absence", "Presence"), values = c("red", "black"), name = "Species state") +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    ylab("Latitude") +
    xlab("Longitude") +
    theme(plot.title = element_text(hjust = 0.5))
  if (display_plot){
    obs_plot
  }
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(here("plot", paste0("Observed_occurences_", gsub("\\.", "_", species_to_plot), ".png")), plot = obs_plot)
  }
  coord_site$theta <- jSDM_binom_pro$theta_latent[, species_to_plot]
  est_plot <- ggplot() +
    geom_sf(data = map) +
    geom_point(data = coord_site[, c("longitude", "latitude", "theta")],
               aes(x = as.numeric(longitude), y = as.numeric(latitude), colour = theta), size = 1) +
    ggtitle(paste0("Estimated probabilities of presence for \n", species_to_plot)) +
    scale_colour_gradientn(colours = rev(rocket(10)), name = "Prediction", limits = c(0, 1)) +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    ylab("Latitude") +
    xlab("Longitude") +
    theme(plot.title = element_text(hjust = 0.5))

  if (display_plot){
    est_plot
  }
  if (save_plot){
    ggsave(here("plot", paste0("estimated_probabilities_", gsub("\\.", "_", species_to_plot), ".png")), plot = est_plot)
  }
}
