library(haven)
library(rgee)
library(sf)
library(dplyr)
library(tidyr)
library(splitstackshape)
#library(sjPlot)
library(ggplot2)
#library(margins)

useremail = "xxx@yyy.zz"

ee_Initialize(email=useremail)

# parameters
noise_floor_parametric_space = c(0.15, 0.2, 0.21, 0.22, 0.25, 0.3)
buffer_size_parametric_space = c(100, 1500, 2500)

args <- expand.grid(noise_floor_parametric_space, buffer_size_parametric_space)

######################
# 1 Minigrid database
######################

output_sens <- list()

for (i in 1:nrow(args)){

noise_floor <- args[i,1]
buffer_size_m <- args[i,2]

print(paste(noise_floor, buffer_size_m))
print((i/18) * 100)

nc <- read_dta("minigrid_11.6.2020.dta")

nc$id = 1:nrow(nc)
nc_bk = nc 

nc <- nc[!is.na(nc$longitude),]
nc_bkk = nc 

nc = nc %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% dplyr::mutate(id = rownames(nc)) %>%  dplyr::select(id, geometry) %>% sf_as_ee()

bufferBy = function(size) {
   function(feature) {
    return(feature$buffer(size))   
}}

nc2 = nc$map(bufferBy(buffer_size_m))

replacement = ee$Image(0)

nl19 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2019-01-01', '2020-01-01')$select('avg_rad')$median()$subtract(0.125)
nl19 = nl19$where(nl19$lt(noise_floor), replacement)

nl18 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2018-01-01', '2019-01-01')$select('avg_rad')$median()$subtract(0.125)
nl18 = nl18$where(nl18$lt(noise_floor), replacement)

nl17 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2017-01-01', '2018-01-01')$select('avg_rad')$median()$subtract(0.125)
nl17 = nl17$where(nl17$lt(noise_floor), replacement)

nl16 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2016-01-01', '2017-01-01')$select('avg_rad')$median()
nl16 = nl16$where(nl16$lt(noise_floor), replacement)

nl15 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2015-01-01', '2016-01-01')$select('avg_rad')$median()
nl15 = nl15$where(nl15$lt(noise_floor), replacement)

nl14 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2014-01-01', '2015-01-01')$select('avg_rad')$median()
nl14 = nl14$where(nl14$lt(noise_floor), replacement)

nl13 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2013-01-01', '2014-01-01')$select('avg_rad')$median()
nl13 = nl13$where(nl13$lt(noise_floor), replacement)

nl12 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2012-01-01', '2013-01-01')$select('avg_rad')$median()
nl12 = nl12$where(nl12$lt(noise_floor), replacement)

nl19_minigrid_2 <- (nl19$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl18_minigrid_2 <- (nl18$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl17_minigrid_2 <- (nl17$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl16_minigrid_2 <- (nl16$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl15_minigrid_2 <- (nl15$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl14_minigrid_2 <- (nl14$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl13_minigrid_2 <- (nl13$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl12_minigrid_2 <- (nl12$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum


#

all = cbind(nl19_minigrid_2, nl18_minigrid_2, nl17_minigrid_2, nl16_minigrid_2, nl15_minigrid_2, nl14_minigrid_2, nl13_minigrid_2, nl12_minigrid_2)

all<- as.data.frame(all)

colnames(all) <- c(2019:2012)

nc <- read.csv("minigrid_11.6.2020.csv")

nc$id = 1:nrow(nc)

nc <- nc[!is.na(nc$longitude),]

all = bind_cols(all, nc)

write.csv(all, paste0("mg_database_with_NTL_data_", noise_floor, "_", buffer_size_m,".csv"))

# 

all$operational = as.character(all$operational)

all_melt <- dplyr::select(all, 1:8, 10, 17, 27, 29) %>% tidyr::gather(key = "key", value="value", 1:8)

all$year_instal <- as.numeric(all$year_instal)

summary_mg = all_melt %>% filter(operational=="Operating" | operational=="Under project") %>% group_by(country, key) %>% summarise(share_detected = sum(value>0)/n())

write.csv(summary_mg, paste0("country_level_detection", noise_floor, "_", buffer_size_m,".csv"))

output_sens[[i]]<-summary_mg
}


# plot detection ranges 

output_sens_unwrapped <- do.call(bind_rows, output_sens)

output_sens_unwrapped$noise <- rep(args$Var1, 216)
output_sens_unwrapped$buffer <-rep(args$Var2, 216)

# output_sens_unwrapped <- do.call(rbind, lapply(list.files(pattern = "country_level_detection"), function(X){read.csv(X)}))
# 
# files <- list.files(pattern = "country_level_detection")
# 
# output_sens_unwrapped$noise <- unlist(lapply(1:length(files), function(i){rep(gsub("_", "", substr(files, 24, 27))[i], 23)}))
# output_sens_unwrapped$buffer <- unlist(lapply(1:length(files), function(i){rep(gsub(".csv", "", sub('.*\\_', '', files))[i], 23)}))

ggplot(output_sens_unwrapped %>% group_by(country) %>% filter(share_detected==max(share_detected)) %>% mutate(morethanhalf = ifelse(share_detected>=0.5, 1, 0)))+
  theme_gray()+
  geom_point(aes(x=country, y=share_detected, colour=as.factor(morethanhalf)), size=2)+
  xlab("Country")+
  ylab("% of operational mini-grids \ndetected through NTL")+
  ggtitle("Mini-grid detection (maximised detection by noise floor, buffer size, year)")+
  geom_hline(yintercept = 0.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_discrete(name="Detected")

ggsave("detection_exercise_maximised.png", last_plot(), device = "png", scale=1.5)

ggplot(output_sens_unwrapped %>% group_by(country) %>% filter(key==2018) %>%  filter(share_detected==max(share_detected)) %>% mutate(morethanhalf = ifelse(share_detected>=0.5, 1, 0)))+
  theme_gray()+
  geom_point(aes(x=country, y=share_detected, colour=as.factor(morethanhalf)), size=2)+
  xlab("Country")+
  ylab("% of operational mini-grids \ndetected through NTL")+
  ggtitle("Mini-grid detection (maximised detection by noise floor, buffer size)")+
  geom_hline(yintercept = 0.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_discrete(name="Detected")


######################
# 2 Sites in Nigerian World Bank survey 
######################

nc1 <- read.csv("control_nigeria.csv")

nc1 <- nc1[!is.na(nc1$LON_DD_MOD),]

nc1 = nc1 %>% st_as_sf(., coords = c("LON_DD_MOD", "LAT_DD_MOD"), crs = 4326)  %>% dplyr::mutate(id = rownames(nc1)) %>%  dplyr::select(id, geometry) %>% sf_as_ee()

bufferBy = function(size) {
  function(feature) {
    return(feature$buffer(size))   
  }}

nc1 = nc1$map(bufferBy(1000))

# 
# nc <- read.csv("all_SSA_MGS.csv")
# nc <- nc[!is.na(nc$X),]
# nc = nc %>% st_as_sf(., coords = c("X", "Y"), crs = 4326)
# 
# nc <- subset(nc, nc$Country=="Nigeria")

# 
# library(ggplot2)
# 
# library(rnaturalearth)
# 
# world<-ne_countries(country = "Nigeria", returnclass="sf")
# 
# ggplot()+
#   geom_sf(data=world, colour="lightgrey")+
#   geom_sf(data=nc1, aes(colour="Survey no access"), size=2, alpha=0.75)+
#   geom_sf(data=nc, aes(colour="Database mini-grids"), size=2, alpha=0.75)+
#   theme(legend.position = "bottom")
# 
# ggsave("map_nigeria.png", last_plot(), device = "png", scale=3)

replacement = ee$Image(0)

nl19 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2019-01-01', '2020-01-01')$select('avg_rad')$median()$subtract(0.15)
nl19 = nl19$where(nl19$lt(noise_floor), replacement)

nl18 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2018-01-01', '2019-01-01')$select('avg_rad')$median()$subtract(0.15)
nl18 = nl18$where(nl18$lt(noise_floor), replacement)

nl17 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2017-01-01', '2018-01-01')$select('avg_rad')$median()$subtract(0.15)
nl17 = nl17$where(nl17$lt(noise_floor), replacement)

nl16 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2016-01-01', '2017-01-01')$select('avg_rad')$median()
nl16 = nl16$where(nl16$lt(noise_floor), replacement)

nl15 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2015-01-01', '2016-01-01')$select('avg_rad')$median()
nl15 = nl15$where(nl15$lt(noise_floor), replacement)

nl14 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2014-01-01', '2015-01-01')$select('avg_rad')$median()
nl14 = nl14$where(nl14$lt(noise_floor), replacement)

nl13 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2013-01-01', '2014-01-01')$select('avg_rad')$median()
nl13 = nl13$where(nl13$lt(noise_floor), replacement)

nl12 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2012-01-01', '2013-01-01')$select('avg_rad')$median()
nl12 = nl12$where(nl12$lt(noise_floor), replacement)

nl19_minigrid_2 <- ee_extract(x = nl19, y = nc1, fun = ee$Reducer$sum(), scale=400)
nl18_minigrid_2 <- ee_extract(x = nl18, y = nc1, fun = ee$Reducer$sum(), scale=400)
nl17_minigrid_2 <- ee_extract(x = nl17, y = nc1, fun = ee$Reducer$sum(), scale=400)
nl16_minigrid_2 <- ee_extract(x = nl16, y = nc1, fun = ee$Reducer$sum(), scale=400)
nl15_minigrid_2 <- ee_extract(x = nl15, y = nc1, fun = ee$Reducer$sum(), scale=400)
nl14_minigrid_2 <- (nl14$reduceRegions(reducer = ee$Reducer$sum(), collection=nc1, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum
nl13_minigrid_2 <- (nl13$reduceRegions(reducer = ee$Reducer$sum(), collection=nc1, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum
nl12_minigrid_2 <- (nl12$reduceRegions(reducer = ee$Reducer$sum(), collection=nc1, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

#

all = cbind(nl19_minigrid_2, nl18_minigrid_2, nl17_minigrid_2, nl16_minigrid_2, nl15_minigrid_2, nl14_minigrid_2, nl13_minigrid_2, nl12_minigrid_2)

colnames(all) <- c(2019:2012)

nc1 <- read.csv("control_nigeria.csv")

nc1 = nc1 %>% st_as_sf(., coords = c("LON_DD_MOD", "LAT_DD_MOD"), crs = 4326) 

nc1$X = st_coordinates(nc1$geometry)[,1]
nc1$Y = st_coordinates(nc1$geometry)[,2]

all = bind_cols(all, nc1)
all$geometry=NULL

write.csv(all, paste0("nigeria_with_NTL_data_", noise_floor, "_", buffer_size_m,".csv"))

# 3 Bind Senegalian survey data (from Jorg Peters) with mini-grid database
library(haven)
library(sf)
library(dplyr)

nc2 <- read.csv("senegal_publiclighting.csv") 
nc2 = nc2 %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) 

nc2$lon = st_coordinates(nc2)[,1]
nc2$lat = st_coordinates(nc2)[,2]

nc3 = nc2 %>% st_as_sf(crs=4326) %>%  dplyr::mutate(id = rownames(nc2)) %>%  dplyr::select(id, geometry) %>% sf_as_ee()

bufferBy = function(size) {
  function(feature) {
    return(feature$buffer(size))   
  }}

prova3 = nc3$map(bufferBy(5000))

nl19_minigrid_2 <- (nl19$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl18_minigrid_2 <- (nl18$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl17_minigrid_2 <- (nl17$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl16_minigrid_2 <- (nl16$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl15_minigrid_2 <- (nl15$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl14_minigrid_2 <- (nl14$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl13_minigrid_2 <- (nl13$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl12_minigrid_2 <- (nl12$reduceRegions(reducer = ee$Reducer$sum(), collection=prova3, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

#

all = cbind(nl19_minigrid_2, nl18_minigrid_2, nl17_minigrid_2, nl16_minigrid_2, nl15_minigrid_2, nl14_minigrid_2, nl13_minigrid_2, nl12_minigrid_2)

colnames(all) <- c(2019:2012)

nc <- read.csv("senegal_publiclighting.csv")

nc$id = 1:nrow(nc)
all <- as.data.frame(all)

all = bind_cols(all, nc)

all <- st_as_sf(all, coords = c("longitude", "latitude"), crs = 4326)

nc <- read.csv("minigrid_11.6.2020.csv")
nc <- nc[!is.na(nc$longitude),]
nc = nc %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326)

library(nngeo)
prova <- st_join(all, nc, join = st_nn, maxdist = 5000, left=T) 

write.csv(prova, paste0("senegal_with_NTL_data_", noise_floor, "_", buffer_size_m,".csv"))


######################
######################

# population density
# load libraries
library(tidyverse)
library(sf)
library(ggrepel)
library(mapview)
library(raster)
library(exactextractr)
library(lubridate)
library(reshape2)
library(R.utils)
library(rgee)
library(stars)
library(ncdf4)
library(fasterize)
library(gdistance)
library(wbstats)
library(haven)

nc <- read_dta("D:/OneDrive - IIASA/MCC/minigrid_11.6.2020.dta")
nc$id = 1:nrow(nc)
nc_bk = nc 

nc <- nc[!is.na(nc$longitude),]
nc_bkk = nc 

nc = nc %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% dplyr::mutate(id = rownames(nc)) %>%  dplyr::select(id, geometry)

setwd("F:/wealth_climate")
pop <- raster("GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif")  

nc <- st_transform(nc, 3395) %>% st_buffer(2500) %>% st_transform(4326)

nc$pop <- exact_extract(pop, nc, fun="sum")
#https://docs.google.com/spreadsheets/d/1JWvMr_iXlLD_qC4WKqrPy0ug_WQxAKALKsL4gBFqZvY/edit?invite=CMy-2swI&ts=60ccadc4#gid=0

nc$area = as.numeric(st_area(nc))/1e6

nc$popdens <- nc$pop/nc$area

# travel time city
cities <- read.csv("worldcities.csv")
cities <- st_as_sf(cities, coords = c("lng", "lat"), crs=4326)
points = as.data.frame(st_coordinates(cities))
temp <- dim(points)
n.points <- temp[1]
T.filename <- 'study.area.T.rds'
T.GC.filename <- 'study.area.T.GC.rds'
T.GC <- readRDS(T.GC.filename)
xy.data.frame <- data.frame()
xy.data.frame[1:n.points,1] <- points[,1]
xy.data.frame[1:n.points,2] <- points[,2]
xy.matrix <- as.matrix(xy.data.frame)
t34_new <- accCost(T.GC, xy.matrix)
nc$city_tt <- exact_extract(t34_new, nc, fun="mean")
nc$city_tt <- ifelse(is.infinite(nc$city_tt), median(nc$city_tt), nc$city_tt)

# cropland 
rasterOptions(tmpdir='F:/wealth_climate/')
cropland = raster("cultivated_land.tif")>0.5
nc$cropland <- exact_extract(cropland, nc, fun="sum")

# distance to grid
grid <- read_sf("D:/OneDrive - IIASA/Current papers/Latent demand air cooling/onsset/Processing/existing_grid.shp")

cropland_rr <- projectRaster(cropland, crs="+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

xy <- st_coordinates(grid)

dist <- distanceFromPoints(cropland_rr, xy[,1:2])

nc$grid_dist <- exact_extract(dist, nc, fun="mean")

##
nc$geometry=NULL

nc_m <- merge(nc_bk, nc, by="id", all.x=T)
nc_m <- nc_m[!is.na(nc_m$longitude),]

write_dta(nc_m, "D:/OneDrive - IIASA/MCC/covariates_mg.dta", version = 14)

#write.csv(nc_m, "D:/OneDrive - IIASA/MCC/covariates_mg.csv")

###########

library(haven)
library(rgee)
library(sf)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(sjPlot)
library(ggplot2)
library(margins)

ee_Initialize()

setwd('D:/OneDrive - IIASA/MCC/dynamic minigrid')

noise_floor <- 0.15
buffer_size_m <- 1000

nc <- read_dta("D:/OneDrive - IIASA/MCC/minigrid_11.6.2020.dta")

nc$id = 1:nrow(nc)
nc_bk = nc 

nc <- nc[!is.na(nc$longitude),]

nc = nc %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% dplyr::mutate(id = rownames(nc)) %>%  dplyr::select(id, geometry) %>% sf_as_ee()

bufferBy = function(size) {
  function(feature) {
    return(feature$buffer(size))   
  }}

nc2 = nc$map(bufferBy(buffer_size_m))

replacement = ee$Image(0)

nl19 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2019-01-01', '2020-01-01')$select('avg_rad')$median()$subtract(0.15)
nl19 = nl19$where(nl19$lt(noise_floor), replacement)
nl19 = nl19$where(nl19$gt(100), replacement)

nl18 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2018-01-01', '2019-01-01')$select('avg_rad')$median()$subtract(0.15)
nl18 = nl18$where(nl18$lt(noise_floor), replacement)
nl18 = nl18$where(nl18$gt(100), replacement)

nl17 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2017-01-01', '2018-01-01')$select('avg_rad')$median()$subtract(0.15)
nl17 = nl17$where(nl17$lt(noise_floor), replacement)
nl17 = nl17$where(nl17$gt(100), replacement)

nl16 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2016-01-01', '2017-01-01')$select('avg_rad')$median()
nl16 = nl16$where(nl16$lt(noise_floor), replacement)
nl16 = nl16$where(nl16$gt(100), replacement)

nl15 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2015-01-01', '2016-01-01')$select('avg_rad')$median()
nl15 = nl15$where(nl15$lt(noise_floor), replacement)
nl15 = nl15$where(nl15$gt(100), replacement)

nl14 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2014-01-01', '2015-01-01')$select('avg_rad')$median()
nl14 = nl14$where(nl14$lt(noise_floor), replacement)
nl14 = nl14$where(nl14$gt(100), replacement)

nl13 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2013-01-01', '2014-01-01')$select('avg_rad')$median()
nl13 = nl13$where(nl13$lt(noise_floor), replacement)
nl13 = nl13$where(nl13$gt(100), replacement)

nl12 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2012-01-01', '2013-01-01')$select('avg_rad')$median()
nl12 = nl12$where(nl12$lt(noise_floor), replacement)
nl12 = nl12$where(nl12$gt(100), replacement)

nl19_minigrid_2 <- (nl19$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl18_minigrid_2 <- (nl18$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl17_minigrid_2 <- (nl17$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl16_minigrid_2 <- (nl16$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl15_minigrid_2 <- (nl15$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl14_minigrid_2 <- (nl14$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl13_minigrid_2 <- (nl13$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

nl12_minigrid_2 <- (nl12$reduceRegions(reducer = ee$Reducer$sum(), collection=nc2, scale=400) %>% ee_as_sf() %>% dplyr::select(sum) %>% st_set_geometry(NULL))$sum

#

all = cbind(nl19_minigrid_2, nl18_minigrid_2, nl17_minigrid_2, nl16_minigrid_2, nl15_minigrid_2, nl14_minigrid_2, nl13_minigrid_2, nl12_minigrid_2)

all<- as.data.frame(all)

colnames(all) <- c(2019:2012)

nc <- read.csv("D:/OneDrive - IIASA/MCC/minigrid_11.6.2020.csv")

nc$id = 1:nrow(nc)

nc <- nc[!is.na(nc$longitude),]

all_n = bind_cols(all, nc)

all_n = all_n %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(3395) %>% st_buffer(1000) %>% st_transform(4326) 

library(raster)

ghs_smod <- raster('D:/OneDrive - IIASA/Current papers/Prod_Uses_Agriculture/Repo/onsset/input/ghs_layer_smod_2015.tif')

all_n$ghs_smod = exactextractr::exact_extract(ghs_smod, all_n, "majority")

all_n$geometry=NULL

write.csv(all_n, paste0("mg_database_with_NTL_data_", noise_floor, "_", buffer_size_m,"_with_urbanisation.csv"))


#############
#############

noise_floor <- 0.15
buffer_size_m <- 1000

nc <- read_dta("D:/OneDrive - IIASA/MCC/minigrid_11.6.2020.dta")

nc$id = 1:nrow(nc)
nc_bk = nc 

nc <- nc[!is.na(nc$longitude),]

nc = nc %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% dplyr::mutate(id = rownames(nc)) %>%  dplyr::select(id, geometry) %>% sf_as_ee()

bufferBy = function(size) {
  function(feature) {
    return(feature$buffer(size))   
  }}

nc2 = nc$map(bufferBy(buffer_size_m))

replacement = ee$Image(0)

##

cropland = ee$Image("USGS/GFSAD1000_V1")
cropland = cropland$mask(cropland$gt(0))

allIndex<-function(start_date="",end_date="", delay=32){
  # define range dates
  start_date <- as.Date(start_date, format="%Y-%m-%d")
  end_date <- as.Date(end_date, format="%Y-%m-%d")
  
  range <- seq(start_date,
               end_date,
               by = "month")
  
  out<-lapply(range, function(i){
    ndvi = ee$ImageCollection('LANDSAT/LC08/C01/T1_32DAY_NDVI')$filterDate(as.character(i), as.character(i+delay))$mosaic()$select('NDVI')$mask(cropland$gt(0))
    
    NDVI <- (ndvi$reduceRegions(reducer = ee$Reducer$mean(), collection=nc2, scale=100) %>% ee_as_sf() %>% dplyr::select(mean) %>% st_set_geometry(NULL))$mean
    
    NDVI
    
  })}  

Sys.time()
prova<-allIndex(start_date="2014-01-01",end_date="2020-01-01", delay=32)
Sys.time()

prova_2 <- as.data.frame(do.call(cbind, prova))

range <- seq(as.Date("2014-01-01", format="%Y-%m-%d"),
             as.Date("2020-01-01", format="%Y-%m-%d"),
             by = "month")

range <- substr(as.character(range), 1, 7)

colnames(prova_2) <- range

all_n <- cbind(all_n, prova_2)

write.csv(all_n, paste0("mg_database_with_NTL_data_", noise_floor, "_", buffer_size_m,"_with_urbanisation_and_NDVI.csv"))

###

library(rgee)
library(sf)
library(dplyr)
library(haven)
ee_Initialize()

setwd('D:/OneDrive - IIASA/MCC')

values = c(0.25)

output <- lapply(values, function(variable){
  print(variable)
  data <- read_stata("senegal_publiclighting.dta")
  
  data$sharemg <- data$minigrid / data$households
  #data <- filter(data, sharemg>0.5)
  
  data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)
  
  data = data %>% st_transform(3395) %>% 
    st_buffer(5000) %>% st_transform(4326) %>% dplyr::select(geometry)
  
  replacement = ee$Image(0)
  
  nl19 =  ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$filterDate('2019-01-01', '2020-01-01')$select('avg_rad')$median()$subtract(0.125)
  nl19 = nl19$where(nl19$lt(variable), replacement)
  
  l <- ee_extract(x = nl19, y = data, fun = ee$Reducer$sum(), scale=450) 
  
  data <- read_stata("senegal_publiclighting.dta")
  
  data$sharemg <- data$minigrid / data$households
  #data <- filter(data, sharemg>0.5)
  
  bind <- bind_cols(l, data)
  
  bind$pl = ifelse(bind$publiclighting_working==1, 1, 0)
  
  bind$threshold = variable
  
  bind
  
})

bind = do.call(rbind, output)

t.test(bind$rad[bind$publiclighting_working==1], bind$rad[bind$publiclighting_working==0])

#######

library(tidyverse)
library(sf)
library(rgee)
ee_Initialize()
setwd("D:/OneDrive - IIASA/MCC/NDVI")

#

buffer_size_m <- 5000
period_measurement <- 7 # days before and after beginning of harvesting where to average EVI

# import crop harvest data (maize)
r <- readxl::read_xlsx("201222_Harvest_Season_SSA_v2.xlsx")
r$iso3c <- countrycode::countrycode(r$Location, 'country.name', 'iso3c')

r <- r %>% group_by(Crop, iso3c) %>% summarise(Harvest.median=mean(Harvest.median)) %>% ungroup() 
r <- dplyr::filter(r, Crop=="Maize")


# for each of the locations in the mini-grid database, calculate for every month between xxxx and yyyy the EVI and NDVI indexes

nc <- read.csv("D:/OneDrive - IIASA/MCC/minigrid_11.6.2020.csv", stringsAsFactors = F)
nc$iso3c <- countrycode::countrycode(nc$country, 'country.name', 'iso3c')

nc$id = 1:nrow(nc)
nc_bk = nc 

nc <- nc[!is.na(nc$iso3c),]
nc <- nc[!is.na(nc$longitude),]
nc <-  nc[!is.na(nc$year_instal),]
#nc <-  nc[nc$`Operating status`=="Operating",]

nc <- merge(nc, r, by="iso3c", all.x=T)
nc <- nc[!is.na(nc$Harvest.median),]
nc$Harvest.median <- as.integer(nc$Harvest.median)

outputs <- list()

get_data <- function(i){
  
  print(i)
  print(timestamp())
  
  nc_f <- nc[i,] 
  
  nc2 = nc_f %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% dplyr::mutate(id = rownames(nc_f)) %>%  dplyr::select(id, geometry) %>% sf_as_ee()
  
  bufferBy = function(size) {
    function(feature) {
      return(feature$buffer(size))   
    }}
  
  nc2 = nc2$map(bufferBy(5000))
  
  #
  cropland = ee$Image("USGS/GFSAD1000_V1")
  cropland = cropland$mask(cropland$gt(0))
  
  dataset = ee$ImageCollection('MODIS/MOD09GA_006_EVI')$filter(ee$Filter$date('2010-01-01', '2020-08-31'))$filterBounds(nc2)
  
  dataset = dataset$select('EVI')
  
  # l <- ee$List$sequence(0, 247)
  # 
  # f0 <- function(x){
  #   start = ee$Date('2010-01-01')$advance(x, 'month')
  #   end = start$advance(1, 'month')
  #   return(ee$ImageCollection(dataset)$filterDate(start, end)$mean())}
  # 
  # dataset = l$map(ee_utils_pyfunc(f0))
  # dataset = ee$ImageCollection(dataset)
  
  #filter = ee$Filter$calendarRange(nc_f$Harvest.median,nc_f$Harvest.median, 'day_of_year')
  #dataset = dataset$filter(filter)
  
  # GET cropland
  
  f1 = function(image) {
    return(image$mask(cropland))   
  }
  
  EVI = dataset$map(f1)
  
  # addDate = function(image){
  #   doy = image$date()$format('yyyy-MM-dd')
  #   doyBand = ee$Image(doy)$rename('doy')
  #   doyBand = doyBand$updateMask(image$select('EVI')$mask())
  #   return(image$addBands(doyBand))}
  # 
  # EVI = EVI$map(addDate)
  
  ## get mean values 
  EVI_processed = EVI$map(function(image){
    return(image$reduceRegions(nc2, ee$Reducer$mean(), 463.3127))})$flatten()
  
  EVI_processed = ee_as_sf(EVI_processed) %>% st_set_geometry(NULL) %>% pull(mean)
  
  #Smothing with adaptive smooth Savitzky-Golay filter 
  library(signal)
  Myorder=1
  EVI_processed_filtered=signal::sgolayfilt(EVI_processed, p=Myorder, n=31, ts=1)
  
  # plot.ts(EVI_processed)
  # plot.ts(EVI_processed_filtered)
  
  return(EVI_processed_filtered)
  
}

library(purrr)
get_data2 <- possibly(get_data, otherwise = NA)

for(i in 1106:nrow(nc)) {
  outputs[[i]] <- get_data2(i)
  save(outputs, file = "outputs.Rdata")
}

# process failed
for(i in which(is.na(outputs))) {
  a <- get_data2(i)
  it <- 1
  while (is.na(a) & it<5){
    a <- get_data2(i)
    it <- it + 1
  }
  print(paste0(i, " success"))
  outputs[[i]] <- a
  save(outputs, file = "outputs_n.Rdata")
}

# calculate yearly metrics
f_ext <- function(days, nyears, hmd, obs){
  o = vector()
  for (i in 1:nyears){
    o[i] <- mean(outputs[[obs]][((hmd-days)+i*365):((hmd+days)+i*365)], na.rm = T)
  }
  return(o)
}

outputs_n <- list()

for (i in 1:length(outputs)){
  outputs_n[[i]] <- f_ext(period_measurement, 10, nc$Harvest.median[i], i)
}

save(outputs, file = "outputs_yearly.Rdata")

p <- as.data.frame(do.call(rbind, outputs_n))
nc <- as.data.frame(cbind(nc, p))

# write to csv
write.csv(nc, "mg_with_evi_110121.csv")
