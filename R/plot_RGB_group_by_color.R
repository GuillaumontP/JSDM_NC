plot_RGB_group_by_color <- function(stars_pixels_group, coord_pixel, country_name = NULL,
                                    country_sf = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Create plots with coordinates from dimensional reduction algorithm and plot pixel group with multiple colors. Plot can be hide and/or save.
  #'
  #' @param stars_pixels_group stars object. with group number in each pixel usefull, NA in others pixels.
  #' @param coord_pixel flot matrix. coordinates of each pixel, obtained with `PCA_or_tSNE_on_pro_pres_est` function.
  #' @param country_name character. optional, add country's border to plots , default is NULL.
  #' @param country_sf sf object. optional, add border to plots, default is NULL.
  #' @param display_plot boolean. Show plots, default is TRUE.
  #' @param save_plot boolean. Save plots in plot folder as .png file, default is FALSE.
  #'
  #' @import stars
  #' @import terra
  #' @import ggplot2
  #' @import rnaturalearth
  #' @import RColorBrewer
  #' @import gridExtra
  #' @importFrom cowplot save_plot
  #' @importFrom grid textGrob gpar
  #' @export

  nb_group <- length(unique(na.omit(c(stars_pixels_group[[1]]))))
  coord_pixel <- coord_pixel[, 1:3]
  cells_pixel <- which( !is.na(stars_pixels_group[[1]]), arr.ind = TRUE)[, 1:2]
  R_stars <- G_stars <- B_stars <- stars_pixels_group
  for (l in 1:3) {
    # set to [0,255]
    coord_pixel[, l] <- coord_pixel[, l] - min(coord_pixel[, l])
    coord_pixel[, l] <- coord_pixel[, l] / max(coord_pixel[, l]) * 255
  }
  for (cell in 1:dim(cells_pixel)[1]) {
    # affect values in [0, 255] to each pixel in stars object
    R_stars[[1]][cells_pixel[cell, 1], cells_pixel[cell, 2]] <- coord_pixel[cell, 1]
    G_stars[[1]][cells_pixel[cell, 1], cells_pixel[cell, 2]] <- coord_pixel[cell, 2]
    B_stars[[1]][cells_pixel[cell, 1], cells_pixel[cell, 2]] <- coord_pixel[cell, 3]
  }
  # RGB plot
  RGB_terra <- terra::rast(c(R_stars, G_stars, B_stars))
  RGB_terra <- colorize(RGB_terra, value = 1:3, to = "col", stretch = "hist")
  col_palette <- coltab(RGB_terra)[[1]][, 1:4]
  col_hex <- rep(0, dim(col_palette)[1])
  for (i in 1:dim(col_palette)[1]) {
    col_hex[i] <- ggtern::rgb2hex(r = col_palette[i, 2], g = col_palette[i, 3], b = col_palette[i, 4])
  }
  RGB_stars <- st_as_stars(RGB_terra)
  # plot groups with mean color of each group
  col_mean_group <- rep(0, nb_group)
  for (group in 1:nb_group){
    col_mean_group[group] <- ggtern::rgb2hex(r = round(mean(R_stars[[1]][which(stars_pixels_group[[1]] == group, arr.ind = TRUE)])),
                                     g = round(mean(G_stars[[1]][which(stars_pixels_group[[1]] == group, arr.ind = TRUE)])),
                                     b = round(mean(B_stars[[1]][which(stars_pixels_group[[1]] == group, arr.ind = TRUE)])))
  }
  names(col_mean_group) <- 1:length(col_mean_group)
  # init ggplot
  R_plot <- ggplot()
  G_plot <- ggplot()
  B_plot <- ggplot()
  RGB_plot <- ggplot()
  group_plot <- ggplot()
  if(!is.null(country_name)){
    map <- ne_countries(scale = 10, returnclass = "sf", country = country_name)[1]
    map <- st_transform(map, st_crs(stars_pixels_group))
    latlon_output <- st_bbox(map)
    R_plot <- R_plot +
      geom_sf(data =  map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    G_plot <- G_plot +
      geom_sf(data =  map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    B_plot <- B_plot +
      geom_sf(data =  map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    RGB_plot <- RGB_plot +
      geom_sf(data =  map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    group_plot <- group_plot +
      geom_sf(data =  map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
  }
  if(!is.null(country_sf)){
    latlon_output <- st_bbox(country_sf)
    map <- country_sf
    R_plot <- R_plot +
      geom_sf(data = map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    G_plot <- G_plot +
      geom_sf(data = map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    B_plot <- B_plot +
      geom_sf(data = map, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    RGB_plot <- RGB_plot +
      geom_sf(data = country_sf, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
    group_plot <- group_plot +
      geom_sf(data = country_sf, colour = "black", fill = "grey") +
      coord_sf(xlim = c(latlon_output[1], latlon_output[3]), ylim = c(latlon_output[2], latlon_output[4]), expand = TRUE)
  }
  R_plot <- R_plot + geom_stars(data = R_stars) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"), na.value = "transparent") +
    theme_bw()
  G_plot <- G_plot + geom_stars(data = G_stars) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Greens"), na.value = "transparent") +
    theme_bw()
  B_plot <- B_plot + geom_stars(data = B_stars) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Blues"), na.value = "transparent") +
    theme_bw()
  title <- textGrob("Pixel coordinates of the first 3 axis \n with color variations", gp = gpar(fontsize = 20, font = 3))
  R_G_B_plot <- arrangeGrob(R_plot, G_plot, B_plot, ncol = 2, nrow = 2, top = title)
  RGB_plot <- RGB_plot + geom_stars(data = RGB_stars) +
    scale_fill_gradientn(colors = col_hex, na.value = "transparent") +
    theme_bw() +
    labs(title = element_text("Pixel coordinates in the first 3 axis \n with RGB color variations")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  group_plot <-group_plot + geom_stars(data = stars_pixels_group, aes(x = x, y = y, fill = as.factor(stars_pixels_group[[1]]))) +
    scale_fill_manual(name = "Groups", values = col_mean_group, na.value = "transparent") +
    theme_bw() +
    labs(title = "Group of pixels \n filled by mean color of each group") +
    theme(plot.title = element_text(hjust = 0.5))
  if (display_plot){
    print(R_G_B_plot)
    print(RGB_plot)
    print(group_plot)
  }
  if (save_plot){
    dir.create(here("plot"), showWarnings = FALSE)
    save_plot(filename = here("plot", "R_G_B.png"), plot = R_G_B_plot, nrow = 3)
    ggsave(filename = here("plot", "RGB.png"), plot = RGB_plot)
    ggsave(here("plot", "color_group.png"), plot = group_plot)
  }
}
