max_min_species_group <- function(list_theta_path, stars_pixels_group, nb_top_species = 20, save_csv = TRUE){
  #' Get species with estimated probabilities of presence highest for each group. Same thing with lowest probabilities.
  #'
  #' @param list_theta_path character vector. full path to Tiff files containing species probabilities of presence to considers.
  #' @param stars_pixels_group stars object. each pixel contains values of its group number.
  #' @param nb_top_species int. Number of species to display for each group, for lists of lowest and highest probabilities, default is 20.
  #' @param save_csv boolean. Allows to save return in output folder in .csv file, default is TRUE.
  #' @return dataframe. For each group, species name with highest probabilities of presence, mean probabilities, sd probabilities. Same for species with lowest probabilities of presence.
  #'
  #' @import terra
  #' @import stars
  #' @import matrixStats
  #' @import readr
  #' @export

  stars_pixels_group <- rast(stars_pixels_group)
  matrix_values <- values(stars_pixels_group, na.rm = TRUE)
  names(matrix_values) <- "group"
  for (i in list_theta_path) {
    # create matrix with probabilities for all species on each pixels
    theta_terra <- rast(i)
    matrix_values <- cbind(matrix_values, values(theta_terra, na.rm = TRUE))
    names(matrix_values) <- c(names(matrix_values), names(theta_terra))
  }
  nb_group <- length(unique(na.omit(c(matrix_values[, 1]))))
  max_species_group <- matrix(0, ncol = nb_group * 3, nrow = nb_top_species)
  min_species_group <- matrix(0, ncol = nb_group * 3, nrow = nb_top_species)
  name_min_max <- NULL
  for (group in 1:nb_group) {
    max_species_group[, 3 * group - 2] <- names(sort(colMeans(matrix_values[matrix_values[, 1] == group, ]), decreasing = TRUE))[1:nb_top_species]
    min_species_group[, 3 * group - 2] <- names(sort(colMeans(matrix_values[matrix_values[, 1] == group, ]), decreasing = FALSE))[1:nb_top_species]
    max_species_group[, 3 * group - 1] <- colMeans(matrix_values[matrix_values[, 1] == group, max_species_group[, 3 * group - 2]])
    min_species_group[, 3 * group - 1] <- colMeans(matrix_values[matrix_values[, 1] == group, min_species_group[, 3 * group - 2]])
    max_species_group[, 3 * group] <- colSds(matrix_values[matrix_values[, 1] == group, max_species_group[, 3 * group - 2]])
    min_species_group[, 3 * group] <- colSds(matrix_values[matrix_values[, 1] == group, min_species_group[, 3 * group - 2]])
    name_min_max <- c(name_min_max, paste0(c("Group_", "Mean_", "Sd_"), group))
  }
  colnames(max_species_group) <- paste0("Max_", name_min_max)
  colnames(min_species_group) <- paste0("Min_", name_min_max)
  min_max_species_group <- data.frame(max_species_group, min_species_group)
  if (save_csv){
    dir.create(here("output"), showWarnings = FALSE)
    write_csv(min_max_species_group, col_names = TRUE, file = here("output", "min_max_species_group.csv"))
  }

  return(min_max_species_group)
}
