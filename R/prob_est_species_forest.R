prob_est_species_forest <- function(alpha_stars, latent_var_stars, jSDM_binom_pro, data_stars){
  #' Create tif files with probabilities of presences for each species on a given area. Results are save in output/theta
  #'
  #' @param alpha_stars stars object. with centered values in 0.
  #' @param latent_var_stars multilayer stars object. with as layer as latent variables. Make sure each latent variable is scale.
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library.
  #' @param data_stars multilayer stars object. with same explanatories variables whom in `jSDM_binom_pro`. Make sure your explanatories variable are scale.
  #'
  #' @import stars
  #' @import stringr
  #' @export

  dir.create(here("output"), showWarnings = FALSE)
  dir.create(here("output", "theta"), showWarnings = FALSE)
  if ( length(data_stars) == 1){
    data_stars <- split(data_stars)
  }
  if ( length(latent_var_stars) == 1){
    latent_var_stars <- split(latent_var_stars)
  }
  n_var <- length(data_stars)
  n_species <- length(colnames(jSDM_binom_pro$theta_latent))
  npart <- ceiling(n_species / 30) # number of species in each file set to 30

  n_latent <- length(jSDM_binom_pro$mcmc.latent)
  lambdas <- matrix(0, n_species, n_latent)
  betas <- matrix(0, n_species, n_var + 1 )
  for (j in 1:n_species){
    lambdas[j, ] <- colMeans(jSDM_binom_pro$mcmc.sp[[j]])[n_var + 1 + 1:n_latent ]
    betas[j, ] <- colMeans(jSDM_binom_pro$mcmc.sp[[j]])[1:(n_var + 1)] # n_var + intercept
  }
  if (length(data_stars) == 1 ){
    data_stars <- split(data_stars)
  }
  colnames(betas) <- colnames(jSDM_binom_pro$mcmc.sp[[1]])[1:(n_var + 1)]
  params_species <- data.frame(betas, lambda = lambdas)
  rownames(params_species) <- colnames(jSDM_binom_pro$model_spec$presence_data)
  predfun <- function(data_stars, params_species, alpha_stars, latent_var_stars, species.range){
    # give prediction for each species
    # in two step
    # First init for first species
    # Second same process but stack on first species output (stars object)

    # Xbeta1
    n_var <- length(data_stars)
    Xbeta1 <- alpha_stars
    Xbeta1[[1]] <- betas[1, "beta_(Intercept)"]
    for (p in 1:n_var) {
      Xbeta1 <- Xbeta1 + data_stars[[p]] * params_species[1, p + 1]
    }

    # Wlambda_1
    Wlambda_1 <- alpha_stars
    Wlambda_1[[1]] <- 0
    for (layer in 1: length(latent_var_stars)) {
      Wlambda_1 <- Wlambda_1 + latent_var_stars[layer,,] * params_species[1, paste0("lambda.", layer)]
    }

    # probit_theta_1
    probit_theta_1 <- Xbeta1 + Wlambda_1 # + alpha_stars
    probit_theta <- probit_theta_1

    ## Other species
    for (j in (species.range[1] + 1):species.range[2]) {

      ## Xbeta_j
      Xbeta_j <- Xbeta1
      Xbeta_j[[1]] <- betas[j, "beta_(Intercept)"]
      for (p in 1:(n_var - 1)) {
        Xbeta_j <- Xbeta_j + data_stars[[p]] * params_species[j, p + 1]
      }

      ## Wlambda_j
      Wlambda_j <- alpha_stars
      Wlambda_j[[1]] <- 0
      for (layer in 1: length(latent_var_stars)) {
        Wlambda_j <- Wlambda_j + latent_var_stars[layer,,] * params_species[j, paste0("lambda.", layer)]
      }

      ## probit_theta_j
      probit_theta_j <- Xbeta_j + Wlambda_j # + alpha_stars
      probit_theta <- c(probit_theta, probit_theta_j)
    }
    names(probit_theta) <- rownames(params_species)[species.range[1]:species.range[2]]
    return(probit_theta)
  }

  first.species <- seq(1, n_species, by = floor(n_species / npart) + 1)
  for (n in 1:npart){
    probit_theta <- predfun(data_stars = data_stars,
                            params_species = params_species,
                            alpha_stars = alpha_stars,
                            latent_var_stars = latent_var_stars,
                            species.range = c(first.species[n], min(n_species, first.species[n] + floor(n_species / npart))))
    theta <- merge(probit_theta)
    for (j in 1:length(split(theta))) {
      theta[[1]][,,j] <- pnorm(theta[[1]][,,j])
    }
    dir.create(here("output"), showWarnings = FALSE)
    dir.create(here("output", "theta"), showWarnings = FALSE)
    write_stars(theta, options = c("COMPRESS=LZW", "PREDICTOR=2"),
                dsn = here("output", "theta", paste0("KNN_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  }
}
