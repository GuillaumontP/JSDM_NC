JSDM_bino_pro <- function(PA, site_data, n_latent = 0, V_beta, mu_beta, nb_species_plot = 2,
                          display_plot = TRUE, save_plot = FALSE){
  #' Run jSDM_binomial_probit from jSDM package and plot results for some species. Plots can be display and/or save.
  #'
  #' @param PA int matrix. Presence Absence for each species (col) on each site (row). Avoid empty row or column
  #' @param site_data float matrix. Explanatories variables for each site without missing values. Same number of row than number of row in `PA`
  #' @param n_latent int. number of latent variables to use in the model, default is 0.
  #' @param V_beta float vector. Variances of Normal priors for the beta parameters.
  #' @param mu_beta float vector. Means of Normal priors for the beta parameters.
  #' @param nb_species_plot number of species whose results are display.
  #' @param display_plot boolean. If TRUE, display all plot, default is TRUE.
  #' @param save_plot boolean. Save in local in .png all plot in folder plot, default is FALSE.
  #' @return jSDM object. output of 'jSDM_binomial_probit' from 'jSDM' library
  #'
  #' @import stars
  #' @import jSDM
  #' @importFrom coda traceplot densplot as.mcmc
  #' @import here
  #' @import graphics
  #' @import grDevices
  #' @export

  jSDM_binom_pro <- jSDM_binomial_probit(
    burnin = 5000,
    mcmc = 10000,
    thin = 10,
    presence_data = data.matrix(PA),
    site_formula = ~.,
    site_data = site_data,
    n_latent = 2,
    site_effect = "random",
    V_lambda = 1,
    V_beta = V_beta,
    mu_beta =mu_beta,
    shape = 0.1,
    rate = 0.1,
    seed = 1234,
    verbose = 1
  )

  top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[nb_species_plot])
  np <- nrow(jSDM_binom_pro$model_spec$beta_start)


  for (j in top_species[1]) {
    for (p in 1:np) {
      par(mfrow = c(1, 2))
      traceplot(as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
      densplot(as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
      mtext(paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p], ", species : ", names(top_species[top_species == j])),
            side = 3, line = - 2, outer = TRUE)
      z = recordPlot()
      dev.off((1 - display_plot) * 2)
      if (save_plot){
        dir.create(here("plot"), showWarnings = FALSE)
        png(here("plot", paste0("beta_jSDM_", p, ".png")))
        replayPlot(z)
        dev.off()
      }
    }
  }


  ## lambda_j of the top five species
  n_latent <- jSDM_binom_pro$model_spec$n_latent

  for (j in top_species[1]) {

    par(mfrow = c(1, 2))
    for (l in 1:n_latent) {
      traceplot(as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
      densplot(as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
      mtext(paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[np + l],", species : ", names(top_species[top_species == j])),
            side = 3, line = - 2, outer = TRUE)
      z = recordPlot()
      dev.off((1 - display_plot) * 2)
      if (save_plot){
        png(here("plot", "lambda_jSDM.png"))
        replayPlot(z)
        dev.off()
      }
    }
  }

  ## Latent variables W_i for the first two sites

  for (l in 1:n_latent) {
    par(mfrow = c(2, 2))
    for (i in 1:2) {
      traceplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                      main = paste0("Latent variable W_", l, ", site ", i))
      densplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                     main = paste0("Latent variable W_", l, ", site ", i))
    }
    z = recordPlot()
    dev.off((1 - display_plot) * 2)
    if (save_plot){
      png(here("plot", paste0("W_", l, "_jSDM.png")))
      replayPlot(z)
      dev.off()
    }
  }

  ## alpha_i of the first two sites
  plot(as.mcmc(jSDM_binom_pro$mcmc.alpha[, 1:2]))
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "alpha_jSDM.png"))
    replayPlot(z)
    dev.off()
  }

  ## V_alpha
  par(mfrow = c(2, 2))
  traceplot(jSDM_binom_pro$mcmc.V_alpha)
  densplot(jSDM_binom_pro$mcmc.V_alpha)
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "V_alpha_Deviance_jSDM.png"))
    replayPlot(z)
    dev.off()
  }

  ## Deviance
  par(mfrow = c(1, 2))
  traceplot(jSDM_binom_pro$mcmc.Deviance)
  densplot(jSDM_binom_pro$mcmc.Deviance)
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "Deviance_jSDM.png"))
    replayPlot(z)
    dev.off()
  }

  ## probit_theta
  par(mfrow = c(2, 1))
  hist(jSDM_binom_pro$probit_theta_latent,
       main = "Predicted probit theta", xlab = "predicted probit theta")
  hist(jSDM_binom_pro$theta_latent,
       main = "Predicted theta", xlab = "predicted theta")
  z = recordPlot()
  dev.off((1 - display_plot) * 2)
  if (save_plot){
    png(here("plot", "theta_jSDM.png"))
    replayPlot(z)
    dev.off()
  }
  return(jSDM_binom_pro)
}
