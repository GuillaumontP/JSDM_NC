plot_species_richness <- function(jSDM_binom_pro, coord_site, country_name = NULL,
                                  latlon_output = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Plot species richnessobserved and  estimated probabilities of presence for all sites. Possibility to save plots
  #'
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @param coord_site dataframe. columns name must contain latitude and longitude as names.
  #' @param country_name character. English name of the country where are inventory sites, default is NULL.
  #' @param latlon_output float vector. inconsistent with `country_name` coord of output box is this format {lon_min, lat_min, lon_max, lat_max}, default is NULL.
  #' @param display_plot boolean. show plot, default is TRUE.
  #' @param save_plot boolean. Write plot in plot folder as .png files, default is FALSE.
  #'
  #' @import stars
  #' @import rnaturalearth
  #' @import ggplot2
  #' @import viridis
  #' @export

  if (is.null(country_name) & is.null(latlon_output)){
    print("Need latlon coordinates")
    stop()
  }
  if (!is.null(country_name)){
    latlon_output <- st_bbox(ne_countries(scale = 10, country = country_name, returnclass = "sf"))
  }
  df_PA <- jSDM_binom_pro$model_spec$presence_data
  species_richness <- data.frame(richness = rowSums(df_PA))
  species_richness$latitude <- coord_site$latitude
  species_richness$longitude <- coord_site$longitude
  species_richness$estimate <- rowSums(jSDM_binom_pro$theta_latent)
  map <- ne_countries(scale = 10, returnclass = "sf")
  ggplot(data = map) +
    geom_sf() +
    geom_point(data = species_richness,
               aes(x = as.numeric(longitude), y = as.numeric(latitude), colour = richness), size = 1) +
    ggtitle("Observed species richness") +
    scale_colour_gradientn(colours = rev(rocket(10)), name = "Number \n of species",
                           limits = c(floor(min(species_richness[, c(1,4)])), ceiling(max(species_richness[, c(1,4)])))) +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    ylab("Latitude") +
    xlab("Longitude") +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(here("plot", "species_richness_observed.png"))
  }

  # Estimated richness

  ggplot(data = map) +
    geom_sf() +
    geom_point(data = species_richness,
               aes(x = as.numeric(longitude), y = as.numeric(latitude), colour = estimate), size = 1) +
    ggtitle("Estimated species richness") +
    scale_colour_gradientn(colours = rev(rocket(10)), name = "Number \n of species",
                           limits = c(floor(min(species_richness[, c(1,4)])), ceiling(max(species_richness[, c(1,4)])))) +
    coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = FALSE) +
    theme_bw() +
    ylab("Latitude") +
    xlab("Longitude") +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    ggsave(here("plot", "species_richness_estimated.png"))
  }

  par(mfrow = c(1,1))
  plot(rowSums(df_PA), rowSums(jSDM_binom_pro$theta_latent),
       main = "Species richness for area \n with less than 50 differents species",
       xlab = "observed", ylab = "estimated",
       cex.main = 1.4, cex.lab = 1.3, xlim = c(0, 50), ylim = c(0, 50))
  abline(a = 0, b = 1, col = 'red')
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "species_richness_observed_vs_estimated.png"), width = 480, height = 480, units = "px")
    replayPlot(z)
    dev.off()
  }
}
