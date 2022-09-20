library(glue)
library(here) 
library(sf) # for spatial sf objects 
library(stars) # for spatial stars objects 
library(insol) # for function daylength
library(stringr)
library(rgdal)
library(ggplot2)
library(viridis)
library(SPEI)
library(terra)
library(geotidy)

##=============
##
## Init
##
##=============

nodat <- -9999
EPSG <- 32622
proj.s <- "EPSG:4326"
proj.t <- paste0("EPSG:", EPSG) 
border <- read_sf(here("guyaHelp", "guyane_border", "limites_guyane500.shp"))[1]
border <- st_transform(border[1], crs = EPSG)
write_sf(border, here("guyaHelp", "border_guyane.shp"))
border <- read_sf(here("guyaHelp", "border_guyane.shp"))
bbox <- round(st_bbox(border))
extent <- glue("{bbox[1]} {bbox[2]} {bbox[3]} {bbox[4]}")
Extent <- writeLines(extent, here("guyaHelp", "extent_short.txt"))
Extent <- readLines(here("guyaHelp", "extent_short.txt"))


##=============
##
## Current data reproject and crop for Guyana
##
##=============

for (file in list.files(here("guyaHelp", "chelsa_v2_1", "temp"), pattern = "1.tif", full.names = FALSE)) {
  sourcefile <- here("guyaHelp", "chelsa_v2_1", "temp", file)
  destfile <- here("guyaHelp", "chelsa_v2_1", paste0(str_extract_all(file, "(?<=_)(.*?)(?=_)")[[1]][1],
                                                      "_",
                                                      str_extract_all(file, "(?<=_)(.*?)(?=_)")[[1]][2],
                                                      ".tif"))
  system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff \\
              {sourcefile} {destfile}"))
  write_stars(st_crop(read_stars(destfile), border), destfile, options = c("COMPRESS=LZW","PREDICTOR=2"))
}

##============
##
## Futur data reproject and crop for Guyana
##
##============

dir.create(here("guyaHelp", "chelsa_v2_1", "futur"))
for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){
  for(var in c("tas_", "pr_"))
  {
    dir.create(here("guyaHelp", "chelsa_v2_1", "futur", model))
    for(i in 1:12)
    {
      sourcefile <- here("data_raw", "chelsa_v2_1", "futur", model, "temp", paste0(var, str_pad(i, 2, pad = "0"), ".tif"))
      destfile <- here("guyaHelp", "chelsa_v2_1", "futur", model, paste0(var, str_pad(i, 2, pad = "0"), ".tif"))
      system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
      write_stars(obj = st_crop(read_stars(destfile), border), options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
                  dsn = destfile)
    }
    for (j in c(1, 4, 12)) {
      sourcefile <- here("data_raw", "chelsa_v2_1", "futur", model, "temp", paste0("bio", j, ".tif"))
      destfile <- here("guyaHelp", "chelsa_v2_1", "futur", model,  paste0("bio", j, ".tif"))
      system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
      write_stars(obj = st_crop(read_stars(destfile), border), options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
                  dsn = destfile)
    }
  }
}
bio_1 <- bio_4 <- bio_12 <- 0
for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"))
{
    bio_model_1 <- read_stars(here("guyaHelp", "chelsa_v2_1", "futur", model, paste0("bio", 1, ".tif"))) 
    bio_1 <- bio_1 + bio_model_1
    bio_model_4 <- read_stars(here("guyaHelp", "chelsa_v2_1", "futur", model, paste0("bio", 4, ".tif"))) 
    bio_4 <- bio_4 + bio_model_4
    bio_model_12 <- read_stars(here("guyaHelp", "chelsa_v2_1", "futur", model, paste0("bio", 12, ".tif"))) 
    bio_12 <- bio_12 + bio_model_12
}

bio_ref_1 <- read_stars(here("guyaHelp", "chelsa_v2_1", paste0("bio", str_pad(1, 2, pad = "0"), "_1981-2010.tif")))
bio_mean_1 <- bio_ref_1 - bio_1 / 5
bio_ref_4 <- read_stars(here("guyaHelp", "chelsa_v2_1", paste0("bio", str_pad(4, 2, pad = "0"), "_1981-2010.tif")))
bio_mean_4 <- bio_ref_4 - bio_4 / 5
bio_ref_12 <- read_stars(here("guyaHelp", "chelsa_v2_1", paste0("bio", str_pad(12, 2, pad = "0"), "_1981-2010.tif")))
bio_mean_12 <- bio_ref_12 - bio_12 / 5 

write_stars(bio_mean_1, here("guyaHelp", "bio_mean_1.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))
write_stars(bio_mean_4, here("guyaHelp", "bio_mean_4.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))
write_stars(bio_mean_12, here("guyaHelp", "bio_mean_12.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))

## CWD 
indicator <- read_stars(here("guyaHelp", "bio_mean_1.tif"))
CWD_current <- read_stars(here("guyaHelp", "bio_mean_1.tif"))
CWD_futur <- read_stars(here("guyaHelp", "bio_mean_1.tif"))
names(CWD_current) <- "CWD_current"
names(CWD_futur) <- "CWD_futur"
CWD_current[[1]] <- 0
CWD_futur[[1]] <- 0
t_month_current <- pr_month_current <- rep(0, 12)
t_month_futur <- pr_month_futur <- rep(0, 12)
# coord of no null cell
coord_latlon <- coordinates(spTransform(as_Spatial(st_as_sf(indicator)), CRS("+proj=longlat +datum=WGS84")))
increment <- 1
for (x in 1:st_dimensions(CWD_current)$x$to) # /!\ 35h to run
{
  for (y in 1:st_dimensions(CWD_current)$y$to) 
  {
    # get T_mean for each cell for each month
    if (!is.na(indicator[[1]][x, y]))
    {
      t_month_futur <- pr_month_futur <- rep(0, 12)
      for (month in 1:12) 
      {
        # t_month_current[month] <- read_stars(here("guyaHelp", "chelsa_v2_1", paste0("tas_", str_pad(month, 2, pad = "0"), ".tif")))[[1]][x, y]
        # pr_month_current[month] <- read_stars(here("guyaHelp", "chelsa_v2_1", paste0("pr_", str_pad(month, 2, pad = "0"), ".tif")))[[1]][x, y]
        for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")) {
          t_month_futur[month] <- t_month_futur + read_stars(here("guyaHelp", "chelsa_v2_1", "futur", model, paste0("tas_", str_pad(month, 2, pad = "0"), ".tif")))[[1]][x, y] / 5
          pr_month_futur[month] <- pr_month_futur + read_stars(here("guyaHelp", "chelsa_v2_1", "futur", model, paste0("pr_", str_pad(month, 2, pad = "0"), ".tif")))[[1]][x, y] / 5
        }
      }
      # CWD_current[[1]][x, y] <- sum(pmax(thornthwaite(Tave = t_month_current, lat = coord_latlon[increment, 2]) - pr_month_current, 0))
      CWD_futur[[1]][x, y] <- sum(pmax(thornthwaite(Tave = t_month_futur, lat = coord_latlon[increment, 2]) - pr_month_futur, 0))
      increment <- increment + 1
    }
  }
}
# write_stars(CWD_current, here("guyaHelp", "CWD_current.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))
write_stars(CWD_futur, here("guyaHelp", "CWD_futur.tif"), options = c("COMPRESS=LZW","PREDICTOR=2"))

##============
##
## Plots
##
##============
bio_mean_1 <- read_stars(here("guyaHelp", "bio_mean_1.tif"))
bio_mean_4 <- read_stars(here("guyaHelp", "bio_mean_4.tif"))
bio_mean_12 <- read_stars(here("guyaHelp", "bio_mean_12.tif"))
bio1 <- read_stars(here("guyaHelp", "chelsa_v2_1", "bio01_1981-2010.tif"))
bio4 <- read_stars(here("guyaHelp", "chelsa_v2_1", "bio04_1981-2010.tif"))
bio12 <- read_stars(here("guyaHelp", "chelsa_v2_1", "bio12_1981-2010.tif"))
CWD_current <- st_crop(read_stars(here("guyaHelp", "CWD_current.tif")), border)
CWD_futur <- st_crop(read_stars(here("guyaHelp", "CWD_futur.tif")), border)
CWD_diff <- CWD_current - CWD_futur
dir.create(here("guyaHelp", "plot"))

# Bio 1
ggplot() +
  geom_stars(data = bio1) +
  scale_fill_gradientn(colours = magma(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("Temperature (°C)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "bio1_current.png"))

ggplot() +
  geom_stars(data = bio_mean_1) +
  scale_fill_gradientn(colours = magma(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("Temperature (°C)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "bio1_diff.png"))

# Bio4
ggplot() +
  geom_stars(data = bio4) +
  scale_fill_gradientn(colours = viridis(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("T seasonality \n (°C x 10 SD)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "bio4_current.png"))

ggplot() +
  geom_stars(data = bio_mean_4) +
  scale_fill_gradientn(colours = viridis(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("T seasonality \n (°C x 10 SD)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "bio4_diff.png"))

# Bio 12
ggplot() +
  geom_stars(data = bio12) +
  scale_fill_gradientn(colours = cividis(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("Precipitation \n (mm/year)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "bio12_current.png"))

ggplot() +
  geom_stars(data = bio_mean_12) +
  scale_fill_gradientn(colours = cividis(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("Precipitation (mm/year)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "bio12_diff.png"))

# CWD
ggplot() +
  geom_stars(data = CWD_current) +
  scale_fill_gradientn(colours = plasma(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("Climatic water deficit \n (mm/year)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "cwd_current.png"))

ggplot() +
  geom_stars(data = CWD_diff) +
  scale_fill_gradientn(colours = plasma(9), na.value = "transparent") + 
  coord_fixed() +
  theme_void() +
  ggtitle("Climatic water deficit \n (mm/year)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5))
ggsave(here("guyaHelp", "plot", "cwd_diff.png"))