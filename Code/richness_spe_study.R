library(stars)
library(ggplot2)
library(here)
library(glue)
library(viridisLite)
# Init
val_max <- 400
richness_spe <- read_stars(here("output", "species_richness.tif"))
coord_high <- which(richness_spe[[1]] >= val_max, arr.ind = TRUE)
dir.create(here("plot", paste0("study_spe_richness_", val_max)), showWarnings = FALSE)

# Alpha
alpha <- read_stars(here("output", "alpha_knn.tif"))
alpha_high <- alpha[[1]][coord_high]
alpha_all <- alpha[[1]][!is.na(alpha[[1]])]
merge_alpha <- data.frame(rbind(cbind(alpha_all, 1), cbind(alpha_high, 0)))

# latent variables 

files <- list.files(here("output"), pattern = "lv_", full.names = TRUE)
for (i in 1:length(files)) {
  lv <- read_stars(files[i])
  lv_high <- lv[[1]][coord_high]
  lv_all <- lv[[1]][!is.na(lv[[1]])]
  assign(paste0("merge_lv_", i), data.frame(rbind(cbind(lv_all, 1), cbind(lv_high, 0))))
}

# explanatories variables

var <- read_stars(here("output", "var_jSDM.tif"))
UM <- split(var[,,,1])
UM_high <- UM[[1]][coord_high]
UM_all <- UM[[1]][!is.na(UM[[1]])]
merge_UM <- data.frame(rbind(cbind(UM_all, "All_pixels"), cbind(UM_high, paste0("Pixels higher than ", val_max))))
merge_UM[,1][merge_UM[,1] == 0] <- "NUM"
merge_UM[,1][merge_UM[,1] == 1] <- "UM"
for (i in 2:((length(names(split(var))) - 1) / 2)) {
  var_high <- split(var[,,, i + 1])[[1]][coord_high]
  var_all <- split(var[,,, i + 1])[[1]][!is.na(split(var[,,, i + 1])[[1]])]
  var_high2 <- split(var[,,, i + 6])[[1]][coord_high]
  var_all2 <- split(var[,,, i + 6])[[1]][!is.na(split(var[,,, i + 6])[[1]])]
  assign(paste0("merge_", names(split(var))[i]), data.frame(rbind(cbind(var_all, 1), cbind(var_high, 0))))
  assign(paste0("merge_", names(split(var))[i], ".2"), data.frame(rbind(cbind(var_all2, 1), cbind(var_high2, 0))))
}

##================
##
## Plot density for multiple argument
##
##================

# alpha

ggplot() + 
  geom_density(data = merge_alpha, aes(x = alpha_all, group = as.factor(V2), color = as.factor(V2)), lwd = 1,) + 
  scale_color_manual(labels = c(glue('Pixels with \n richness >= {val_max}'), 'All pixels'), values = viridis(2)) + 
  theme_bw() +
  labs(x = "Alpha values", color = "Type of pixel", shape = "Type of pixel") +
  ggtitle(glue("Distribution of alpha values for pixels \n with more than {val_max} species predicted \n and all values in New Caledonia")) + 
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5))
ggsave(here("plot", paste0("study_spe_richness_", val_max), "alpha.png"))
# latent variables

for (i in 1:length(files)) {
  p <- ggplot() + 
    geom_density(data = eval(parse(text = paste0("merge_lv_", i))), aes(x = lv_all, group = as.factor(V2), color = as.factor(V2)), lwd = 1,) + 
    scale_color_manual(labels = c(glue('Pixels with \n richness >= {val_max}'), 'All pixels'), values = viridis(2)) + 
    theme_bw() +
    labs(x = paste0("lv ", i, " values"), color = "Type of pixel", shape = "Type of pixel") +
    ggtitle(glue("Distribution of {paste0('latent values ', i, ' values')} for pixels \n with more than {val_max} species predicted \n and all values in New Caledonia")) + 
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5))
  print(p)
  ggsave(filename = here("plot", paste0("study_spe_richness_", val_max), paste0("lv_", i, ".png")))
}

# explanatories variables 

ggplot() +
  geom_bar(position="fill", data = merge_UM, aes(x = as.factor(V2), fill = as.factor(UM_all)), lwd = 1) +
  scale_fill_manual(labels = c("No ultramafics soil", 'Ultramafics soil'), values = viridis(2), name = "Type of soil") +
  theme_bw() +
  labs(x = "ultramafics status") +
  ggtitle(glue("Distribution of ultramafics soil for pixels \n with more than {val_max} species predicted \n and all values in New Caledonia")) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5))
ggsave(filename = here("plot", paste0("study_spe_richness_", val_max), "UM.png"))

for (i in 2:((length(names(split(var))) - 1) / 2)) {
  p <- ggplot() + 
    geom_density(data = eval(parse(text = paste0("merge_", names(split(var))[i]))), aes(x = var_all, group = as.factor(V2), color = as.factor(V2)), lwd = 1,) + 
    scale_color_manual(labels = c(glue('Pixels with \n richness >= {val_max}'), 'All pixels'), values = viridis(2)) + 
    theme_bw() +
    labs(x = paste0(names(split(var))[i], " values"), color = "Type of pixel", shape = "Type of pixel") +
    ggtitle(glue("Distribution of {names(split(var))[i]} values for pixels \n with more than {val_max} species predicted \n and all values in New Caledonia")) + 
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5))
  print(p)
  ggsave(filename = here("plot", paste0("study_spe_richness_", val_max), paste0(names(split(var))[i], ".png")))  
  
  g <- ggplot() + 
    geom_density(data = eval(parse(text = paste0("merge_", names(split(var))[i], ".2"))), aes(x = var_all2, group = as.factor(V2), color = as.factor(V2)), lwd = 1,) + 
    scale_color_manual(labels = c(glue('Pixels with \n richness >= {val_max}'), 'All pixels'), values = viridis(2)) + 
    theme_bw() +
    labs(x = paste0(names(split(var))[i], ".2", " values"), color = "Type of pixel", shape = "Type of pixel") +
    ggtitle(glue("Distribution of {names(split(var))[i]}.2 values for pixels \n with more than {val_max} species predicted \n and all values in New Caledonia")) + 
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5))
  print(g)
  ggsave(filename = here("plot", paste0("study_spe_richness_", val_max), paste0(names(split(var))[i], "_2.png")))  
}





