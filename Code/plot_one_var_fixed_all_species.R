load(here("output", "RData", "jSDM_binom_pro.RData"))
beta <- matrix(0, nrow = 12, ncol = dim(jSDM_binom_pro$model_spec$presence_data)[2])
colnames(beta) <- colnames(jSDM_binom_pro$model_spec$presence_data)
for (i in 1:dim(jSDM_binom_pro$model_spec$presence_data)[2]) {
  beta[, i] <- colMeans(jSDM_binom_pro$mcmc.sp[[i]])[1:12]
}
mins <- rowMins(beta) # min of each var
maxs <- rowMaxs(beta) # max of each var
names(mins) <- names(maxs) <- colnames(jSDM_binom_pro$mcmc.sp[[i]])[1:12]
rownames(beta) <- names(maxs)
for (species in 1:dim(beta)[2]) {
  for (var in 3:(dim(beta)[1] / 2 + 1)) { # don't use intercep & ultramafic
    values <- seq(mins[var], maxs[var], length = 100)
    plot(values, values * beta[var, species] + values^2 * beta[var + 5, species])
    ggplot() + 
      geom_point(aes(x = values, y = values * beta[var, species] + values^2 * beta[var + 5, species]), color = "salmon") +
      theme_bw() +
      labs(x = paste(colnames(jSDM_binom_pro$mcmc.sp[[i]])[var], "values"),
           y = "",
           title = paste0(colnames(beta)[species], " with only \n", colnames(jSDM_binom_pro$mcmc.sp[[i]])[var], " moving")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
}
