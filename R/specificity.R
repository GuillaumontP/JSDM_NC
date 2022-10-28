specificity  <- function(PA, jSDM_binom_pro){
  #' Process Specificity (species rightly predicted as absent) on inventory sites.
  #'
  #' @param PA dataframe. containing only 0 and 1 with colnames with species names.
  #' @param jSDM_binom_pro object of class jSDM. output of "jSDM_binomial_probit" of "jSDM" library
  #' @return float vector. percentage of species rightly predicted as absent on each inventory site.
  #' @export

  theta <- jSDM_binom_pro$theta_latent
  n_sites <- nrow(PA)
  score <- rep(0, n_sites)
  for(i in 1:n_sites){
    # Specificity
    # True negative
    abs_sp <- which(PA[i,] == 0)
    nabs_sp <- length(abs_sp)
    pred_abs_sp <- which(theta[i,] <= sort(theta[i,])[nabs_sp])
    score[i] <- sum(pred_abs_sp %in% abs_sp) / ifelse(nabs_sp !=0, nabs_sp, 1)
  }
  return(score)
}
