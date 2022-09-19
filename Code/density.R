library(here)
library(ggplot2)
alpha_forest <- read_stars(here("output", "alpha_forest.tif"))
alpha_forest_values <- alpha_forest[[1]][!is.na(alpha_forest[[1]])]
diversity_alpha <- read_stars(here("output", "theta_forest_sum.tif"))
diversity_alpha_value <- diversity_alpha[[1]][!is.na(diversity_alpha[[1]])]
alpha_div_alpha_300 <- alpha_forest_values[diversity_alpha_value >= 300]
alpha_div_alpha_600 <- alpha_forest_values[diversity_alpha_value >= 600]
load(here("output", "RData", "jSDM_binom_pro.RData"))

## Alpha sites & forest & div alpha > 300 Density
colors <- viridis(4)
names(colors) <- c( "alpha inventory site", "alpha NC forest", "alpha values div alpha > 300", "alpha values div alpha > 600")
ggplot() +
  geom_density(aes(x = colMeans(jSDM_binom_pro$mcmc.alpha), color = "alpha inventory site"), size = 1.2) +
  geom_vline(aes(xintercept = mean(colMeans(jSDM_binom_pro$mcmc.alpha), na.rm = T)),
             color = colors[1], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = alpha_forest_values, color = "alpha NC forest"), size = 1.2) +
  geom_vline(aes(xintercept = mean(alpha_forest_values, na.rm = T)),
             color = colors[2], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = alpha_div_alpha_300, color = "alpha values div alpha > 300"), size = 1.2) +
  geom_vline(aes(xintercept = mean(alpha_div_alpha_300, na.rm = T)),
             color = colors[3], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = alpha_div_alpha_600, color = "alpha values div alpha > 600"), size = 1.2) +
  geom_vline(aes(xintercept = mean(alpha_div_alpha_600, na.rm = T)),
             color = colors[4], size = 1, linetype = "dashed" ) +
  labs(x = "Alpha values", y = "Density", colors = "Legend", title = "Alpha values on inventory sites & NC forest & div alpha > 300") +
  scale_color_manual(values = colors) + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "density_alpha_site_vs_forest.png"))

## min diversity alpha 
min_div_alpha <- min(diversity_alpha_value)

## W_1 sites & forest & div alpha > 300 Density 
W1_forest <- read_stars(here("output", "lv_W1_forest.tif"))
W1_forest_values <- W1_forest[[1]][!is.na(W1_forest[[1]])]
diversity_alpha <- read_stars(here("output", "theta_forest_sum.tif"))
diversity_alpha_value <- diversity_alpha[[1]][!is.na(diversity_alpha[[1]])]
W1_div_alpha_300 <- W1_forest_values[diversity_alpha_value >= 300]
W1_div_alpha_600 <- W1_forest_values[diversity_alpha_value >= 600]

colors <- viridis(4)
names(colors) <- c( "W1 inventory site", "W1 NC forest", "W1 values div alpha > 300", "W1 values div alpha > 600")
ggplot() +
  geom_density(aes(x = colMeans(jSDM_binom_pro$mcmc.latent$lv_1), color = "W1 inventory site"), size = 1.2) +
  geom_vline(aes(xintercept = mean(colMeans(jSDM_binom_pro$mcmc.latent$lv_1), na.rm = T)),
             color = colors[1], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = W1_forest_values, color = "W1 NC forest"), size = 1.2) +
  geom_vline(aes(xintercept = mean(W1_forest_values, na.rm = T)),
             color = colors[2], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = W1_div_alpha_300, color = "W1 values div alpha > 300"), size = 1.2) +
  geom_vline(aes(xintercept = mean(W1_div_alpha_300, na.rm = T)),
             color = colors[3], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = W1_div_alpha_600, color = "W1 values div alpha > 600"), size = 1.2) +
  geom_vline(aes(xintercept = mean(W1_div_alpha_600, na.rm = T)),
             color = colors[4], size = 1, linetype = "dashed" ) +
  labs(x = "W1 values", y = "Density", colors = "Legend", title = "W1 values on inventory sites & NC forest & div alpha > 300") +
  scale_color_manual(values = colors) + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "density_W1_site_vs_forest.png"))

## W_2 sites & forest & div alpha > 300 Density 
W2_forest <- read_stars(here("output", "lv_W2_forest.tif"))
W2_forest_values <- W2_forest[[1]][!is.na(W2_forest[[1]])]
diversity_alpha <- read_stars(here("output", "theta_forest_sum.tif"))
diversity_alpha_value <- diversity_alpha[[1]][!is.na(diversity_alpha[[1]])]
W2_div_alpha_300 <- W2_forest_values[diversity_alpha_value >= 300]
W2_div_alpha_600 <- W2_forest_values[diversity_alpha_value >= 600]

colors <- viridis(4)
names(colors) <- c( "W2 inventory site", "W2 NC forest", "W2 values div alpha > 300", "W2 values div alpha > 600")
ggplot() +
  geom_density(aes(x = colMeans(jSDM_binom_pro$mcmc.latent$lv_1), color = "W2 inventory site"), size = 1.2) +
  geom_vline(aes(xintercept = mean(colMeans(jSDM_binom_pro$mcmc.latent$lv_1), na.rm = T)),
             color = colors[1], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = W2_forest_values, color = "W2 NC forest"), size = 1.2) +
  geom_vline(aes(xintercept = mean(W2_forest_values, na.rm = T)),
             color = colors[2], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = W2_div_alpha_300, color = "W2 values div alpha > 300"), size = 1.2) +
  geom_vline(aes(xintercept = mean(W2_div_alpha_300, na.rm = T)),
             color = colors[3], size = 1, linetype = "dashed" ) +
  geom_density(aes(x = W2_div_alpha_300, color = "W2 values div alpha > 600"), size = 1.2) +
  geom_vline(aes(xintercept = mean(W2_div_alpha_600, na.rm = T)),
             color = colors[4], size = 1, linetype = "dashed" ) +
  labs(x = "W2 values", y = "Density", colors = "Legend", title = "W2 values on inventory sites & NC forest & div alpha > 300") +
  scale_color_manual(values = colors) + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "density_W2_site_vs_forest.png"))

colors <- viridis(2)
names(colors) <- c("lambda 1", "lambda 2")
lambda_mean <- matrix(0, ncol = 2, nrow = length(jSDM_binom_pro$mcmc.sp))
for (i in 1:length(jSDM_binom_pro$mcmc.sp)) {
  lambda_mean[i,] <- colMeans(jSDM_binom_pro$mcmc.sp[[i]][,14:15])
}
ggplot() +
  geom_density(aes(x = lambda_mean[,1], color = "lambda 1"), size = 1.2) +
  geom_vline(aes(xintercept = mean(lambda_mean[,1])),
            color = colors[1], size = 1, linetype = "dashed") +
  geom_density(aes(x = lambda_mean[,2], color = "lambda 2"), size = 1.2) +
  geom_vline(aes(xintercept = mean(lambda_mean[,2])),
             color =  colors[2], size = 1, linetype = "dashed") +
  labs(x = "Lambda", y = "Density", title = "Density of lambda 1 & lambda 2", color = "Legend") +
  scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here("output", "plot", "density_lambda_1_&_2.png"))

var_site <- read.csv2(here("data_raw", "NCpippn", "var_site.csv"), sep = ",", dec = ".")[,2:12]
moist_forest <- read_sf(here("output", "moist_forest.shp"))
colors <- viridis(2)
names(colors) <- c("inventory site", "NC forest")
for (i in 1:dim(var_site)[2]) {
  var_forest <- st_crop(read_stars(here("output", "jSDM_data_final.tif")), moist_forest)
  var_forest <- var_forest[[1]][,, i][!is.na(var_forest[[1]][,,i])]
  print(ggplot() + 
    geom_density(aes(x = var_forest, color = "NC forest"), size = 1.2) +
    geom_vline(aes(xintercept = mean(var_forest)),
               color = colors[1], size = 1, linetype = "dashed") +
    geom_density(aes(x = var_site[, i], color = "inventory site"), size = 1.2) +
    geom_vline(aes(xintercept = mean(var_site[, i])),
               color =  colors[2], size = 1, linetype = "dashed") +
    labs(x = colnames(var_site)[i], y = "Density", color = "Legend", 
         title = paste0("Density of ", colnames(var_site)[i]," on inventory site & NC forest")) +
    scale_color_manual(values = colors) +
    theme(plot.title = element_text(hjust = 0.5)))
  ggsave(here("output", "plot", paste0("Density on ", colnames(var_site)[i], " inventory_site_&_NC_forest.png")))
 }

