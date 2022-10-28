extract_var_JSDM <- function(stars_object, variables_names, power_variable = rep(1, length(variables_names)),
                             scale = rep(TRUE, length(variables_names))){
  #' Extract layers and create power wanted of each variable
  #'
  #' @param stars_object stars object. output of "read_stars" of "stars" library with layer's names.
  #' @param variables_names a character vector. of layer's name needed
  #' @param power_variable int vector same length then `variables_names` . maximal power needed for each variable
  #' @param scale boolean vector same length then `variables_names` . scale variables and its power
  #' @return stars object with only variables and power asked.
  #'
  #' @import stars
  #' @export

  if( length(variables_names) != length(power_variable)){
    print("variables_names length is different then power_variable length")
    stop()
  }
  if (length(names(st_dimensions(stars_object))) != 2){
    stars_object <- split(stars_object)
  }
  stars_output <- stars_object[variables_names,,]
  for (i in 1:length(variables_names)) {
    for (power in 1:power_variable[i]) {
      if ( power != 1){
        var_power <- stars_output[variables_names[i],,]^power
        if (scale[i]){
          var_scale <- scale(c(var_power[[1]]))
          var_power[[1]] <- (var_power[[1]] - attr(var_scale, "scaled:center")) / attr(var_scale, "scaled:scale")
        }
        names_stars <- names(stars_output)
        stars_output <- c(stars_output, var_power)
        names(stars_output) <- c(names_stars, paste0(variables_names[i], "^", power))
      }
    }
    if (scale[i]){
      var_scale <- scale(c(stars_output[variables_names[i],,][[1]]))
      variable_position <- which(names(stars_output) == variables_names[i])
      stars_output <- merge(stars_output)
      stars_output[[1]][,, variable_position] <- (stars_output[[1]][,, variable_position] - attr(var_scale, "scaled:center")) / attr(var_scale, "scaled:scale")
      stars_output <- split(stars_output)
    }
  }

  return(stars_output)
}
