sensitivity <- function(PA, jSDM_binom_pro){
  #' Process Sensitivity (species rightly predicted as present) on inventory sites.
  #'
  #' @param PA dataframe. containing only 0 and 1 with colnames with species names.
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @return float vector. percentage of species rightly predicted as present on each inventory site.
  #' @export

  theta <- jSDM_binom_pro$theta_latent
  n_sites <- nrow(PA)
  score <- rep(0, n_sites)
  for(i in 1:n_sites){
    # Sensitivity
    # True positive
    obs_sp <- which(PA[i,] > 0)
    nobs_sp <- length(obs_sp)
    pred_sp <- which(theta[i,] >= sort(theta[i,], decreasing = TRUE)[nobs_sp])
    score[i] <- sum(pred_sp %in% obs_sp) / ifelse(nobs_sp != 0, nobs_sp, 1)
  }
  return(score)
}
