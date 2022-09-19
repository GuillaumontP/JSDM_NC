library(here)
library(sp)
library(jSDM)
library(stars)
library(glue)
library(rgdal)
library(readr)
library("rnaturalearth") # plotting maps
library("rnaturalearthdata")
library("rnaturalearthhires")
library(ggplot2)
library(viridis)
library(rgrass7)
library(terra)
library(stringr)
library(ade4)
library(parallel)
library(doParallel)
library(factoextra)
library(RColorBrewer)
set.seed(1234)
EPSG <- 3163
nodat <- -9999
ISO_country_code <- "NCL"
proj.t <- paste0("EPSG:", EPSG) 
# a = rep(0,554)
# data_clear <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep = ",")
# k = 1
# for ( i in unique(data_clear$id_locality)){
#   a[k] = data_clear$plot[data_clear$id_locality == i][1]
#   k = k + 1
# }


coord = read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")[,3:4]
data_site = matrix(0, ncol = 19, nrow = 555)
for (i in 1:dim(coord)[1]) 
{
  data_site[i, 1:19] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("data_raw", "chelsa_v2_1", "bio_1km.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
bio1_12 = data.frame(data_site[, c(1, 12)])
names(bio1_12) <- c("bio1", "bio12")
bio1_12$bio1 = bio1_12$bio1 / 10

write_stars(read_stars(here("data_raw", "chelsa_v2_1", "bio_1km.tif"))[,800:1300,200:570,], dsn = here("output", "plot", "pr_VS_temp.tif"))
datajsdm = terra::rast(here("output", "plot", "pr_VS_temp.tif"))
samp_theta <- spatSample(datajsdm, size = 10000, cells = TRUE, na.rm = TRUE, method = "random")
vs = samp_theta[,c(2,13)] 
vs[,1] <- vs[,1] / 10

# bio 1 vs bio 12
ggplot() +
  geom_density_2d(data = vs, aes(x = bio1, y = bio12, col = " Nouvelle-Calédonie")) +
  geom_density2d(data = bio1_12, aes(x = bio1, y = bio12, col = "Sites d'inventaire")) +
  theme_bw() +
  labs(color = "") +
  xlab("Température (en °C)") +
  ylab("Précipitation (en kg/m²/an)") +
  xlim(c(18, 25)) +
  ylim(c(830, 3400)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20))
ggsave(here("output", "plot", "Temp_Pr_NC_site.png"))

## plot bio4 saisonnalité des températures site Vs NC

bio4 = data.frame(samp_theta[,5])
names(bio4) = "bio_4"
bio4_site = data.frame(data_site[,4])
names(bio4_site) = "bio_site_4"

# ggplot() +
#   geom_histogram(data = bio4_site, aes(x = bio_site_4),  binwidth = ) +
#   ggtitle("") +
#   theme_bw() +
#   scale_y_continuous(labels = scales::percent) +
#   theme(plot.title = element_text(size = 15))

ggplot() + 
  geom_density(data = bio4, aes(x = bio_4, color = "density", fill = "density")) +
  geom_histogram(data = bio4_site, aes(x = bio_site_4, y = ..density..,
                 color = "hist", fill = "hist"), binwidth = 7, alpha = 0.8, stat = "bin") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("Variance des précipitation (kg/m²/an)") +
  xlim(c(1500,2700)) +
  ylab("") + 
  scale_colour_manual(name = "", values = c("hist" = "#E83F3FFF", "density" = "#961C5BFF"),
                      labels = c("hist" = "Site d'inventaire", "density" = "Nouvelle-Calédonie")) +
  scale_fill_manual(name = "", values = c("hist" = "#E83F3FFF", "density" = "#961C5BFF"),
                    labels = c("hist" = "Site d'inventaire", "density" = "Nouvelle-Calédonie")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20))
ggsave(here("output", "plot", "Variance_prec.png"))


##========= plot pb CV method CAH EM KM
point <- data.frame(matrix(c(1, 2, 3, 2, 2.25, 2, 4, 2), ncol = 2, byrow = TRUE)) 
point$paired <- c(1, 1, 1, 2)
point$group <- c("Groupe de référence", "Groupe de référence", "Groupe estimé", "Groupe estimé")
ggplot() +
  geom_line(data = point, aes(x = X1, y = X2, group = paired), linetype = "dashed", color = "black", size = 2) +
  geom_point(data = point, aes(x = X1, y = X2, color = as.factor(group)), size = 7) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(6, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.text = element_text(size = 20))
ggsave(here("output", "plot", "min_simple_cluster.png"))

point$paired <- c(1, 2, 1, 2)

ggplot() +
  geom_line(data = point, aes(x = X1, y = X2, group = paired), linetype = "dashed", color = "black", size = 2) +
  geom_point(data = point, aes(x = X1, y = X2, color = as.factor(group)), size = 7) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(6, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.text = element_text(size = 20))
ggsave(here("output", "plot", "min_all_cluster.png"))

##============ plot variable clim & env 
# ndm, bio1
# peridotites, distSea
border <- st_read(here("data_raw", "Grande_Terre", "Grande_Terre.shp"))
cur <- read_stars(here("output", "current_chelsaNC.tif"))
env <- st_crop(read_stars(here("output", "environ_allNC.tif")), border)

png(here("output", "plot", "stars_bio1.png"))
cur_terra73 <- as(cur[,,,73], "Raster")
plot(rast(cur_terra73), axes = FALSE)
dev.off()

cur_terra93 <- as(cur[,,,93], "Raster")
png(here("output", "plot", "stars_ndm.png"))
plot(rast(cur_terra93), axes = FALSE)
dev.off()

env_terra15 <- as(env[,,,15], "Raster")
png(here("output", "plot", "stars_peridotites.png"))
plot(rast(env_terra15), axes = FALSE)
dev.off()

env_terra9 <- as(env[,,,9], "Raster")
png(here("output", "plot", "stars_distSea.png"))
plot(rast(env_terra9), axes = FALSE)
dev.off()
