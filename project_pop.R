
a <- raster("baseYr_total_2000.nc4") # MAKE AT 1 km
b <- raster("ssp2_total_2050.nc4")
c <- raster("ssp5_total_2050.nc4")
d <- raster("ssp1_total_2050.tif")
e <- raster("ssp3_total_2050.tif")

rasters <- list.files(pattern="_1km", path="AGEPOP", full.names = T, recursive = T)
rasters <- lapply(rasters, raster)

###############################

r1 <- stack(rasters)

# a <- projectRaster(a, r1[[1]])
# b <- projectRaster(b, r1[[1]])
# c <- projectRaster(c, r1[[1]])

b <- (b/a)
c <- (c/a)
d <- (d/a)
e <- (e/a)

# project current worldpop using gao growth rates

names(r1) <- read_rds("names_r.rds")

r1_g <- stack(r1[[1]])
r1_g[[1]] <- sum(r1[[c(1:4, 11)]]) + sum(r1[[c(19:22, 29)]]) 
r1_g[[2]] <- sum(r1[[c(5:10, 12:15)]])+ sum(r1[[c(23:28, 30:33)]])
r1_g[[3]] <-  sum(r1[[16:18]]) + sum(r1[[34:36]])

gc()

#b <- aggregate(b, fact=30, fun="mean")
b <- projectRaster(b, r1_g)

r1_g_245 <- r1_g * (b) # make all pop groups grow at the same rate


#c <- aggregate(c, fact=30, fun="mean")
c <- projectRaster(c, r1_g)

r1_g_585 <- r1_g * (c) # make all pop groups grow at the same rate

#d <- aggregate(d, fact=30, fun="mean")
d <- projectRaster(d, r1_g)

r1_g_126 <- r1_g * (d)

#e <- aggregate(e, fact=30, fun="mean")
e <- projectRaster(e, r1_g)

r1_g_370 <- r1_g * (e)


# reduce resolution (and computational burden)

r1_g <- aggregate(r1_g, fact=5, fun="sum")
r1_g_245 <- aggregate(r1_g_245, fact=5, fun="sum")
r1_g_585 <- aggregate(r1_g_585, fact=5, fun="sum")
r1_g_126 <- aggregate(r1_g_126, fact=5, fun="sum")
r1_g_370 <- aggregate(r1_g_370, fact=5, fun="sum")

r1_g <- stack(r1_g)
r1_g_245 <- stack(r1_g_245)
r1_g_585 <- stack(r1_g_585)
r1_g_126 <- stack(r1_g_126)
r1_g_370 <- stack(r1_g_370)

# harmonise total pop per country in 2050 with iiasa value (adjustment)

ssp_stats <- read.csv("SspDb_country_data_2013-06-12.csv")

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1", "SSP2", "SSP3", "SSP5"))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")
ssp_stats <- dplyr::select(ssp_stats, SCENARIO, REGION, X2050)

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)

countries_sf <- merge(wrld_simpl, ssp_stats, by.x="ISO3", by.y="REGION")

countries_sf$pop_245 <- exact_extract(sum(r1_g_245), countries_sf, "sum", max_cells_in_memory=1e09)
countries_sf$adj_fact_245 <- countries_sf$pop_245/(countries_sf$X2050 * 1e6)

countries_sf$pop_245[countries_sf$SCENARIO=="SSP5"] <- NA
countries_sf$adj_fact_245[countries_sf$SCENARIO=="SSP5"] <- NA

countries_sf$pop_585 <- exact_extract(sum(r1_g_585), countries_sf, "sum", max_cells_in_memory=1e09)
countries_sf$adj_fact_585 <- countries_sf$pop_585/(countries_sf$X2050 * 1e6)

countries_sf$pop_585[countries_sf$SCENARIO=="SSP2"] <- NA
countries_sf$adj_fact_585[countries_sf$SCENARIO=="SSP2"] <- NA

countries_sf$pop_126 <- exact_extract(sum(r1_g_126), countries_sf, "sum", max_cells_in_memory=1e09)
countries_sf$adj_fact_126 <- countries_sf$pop_126/(countries_sf$X2050 * 1e6)

countries_sf$pop_370 <- exact_extract(sum(r1_g_370), countries_sf, "sum", max_cells_in_memory=1e09)
countries_sf$adj_fact_370 <- countries_sf$pop_370/(countries_sf$X2050 * 1e6)

# df with cell value and country of belonging to apply adjustment factor


wrld_simpl$REGION_n <- as.numeric(wrld_simpl$ISO3)
regions_r <- fasterize::fasterize(wrld_simpl, raster(r1_g_245[[1]]), field="REGION_n")
regions <- dplyr::select(wrld_simpl, ISO3, REGION_n)
regions$geometry<-NULL
regions$reg <- regions$REGION_n
regions$REGION_n <- NULL

df1 <- data.frame(values(r1_g_245), reg = values(regions_r))
df1$id  <- 1:nrow(df1)
df1 <-  plyr::join(df1, regions, by="reg", type = "left")
df1 <- df1[order(df1$id), ]

countries_sf$geometry <- NULL

df1 <- plyr::join(df1, countries_sf %>% dplyr::filter(SCENARIO=="SSP2") , by="ISO3", type = "left")
df1 <- df1[order(df1$id), ]

df1$layer.1 <- df1$layer.1 / df1$adj_fact_245
df1$layer.2 <- df1$layer.2 / df1$adj_fact_245
df1$layer.3 <- df1$layer.3 / df1$adj_fact_245

values(r1_g_245[[1]]) <- df1$layer.1 
values(r1_g_245[[2]]) <- df1$layer.2
values(r1_g_245[[3]]) <- df1$layer.3


# df with cell value and country of belonging to apply adjustment factor


wrld_simpl$REGION_n <- as.numeric(wrld_simpl$ISO3)
regions_r <- fasterize::fasterize(wrld_simpl, raster(r1_g_585[[1]]), field="REGION_n")
regions <- dplyr::select(wrld_simpl, ISO3, REGION_n)
regions$geometry<-NULL

df2 <- data.frame(values(r1_g_585), reg = values(regions_r))
df2$id  <- 1:nrow(df2)
df2 <- merge(df2, regions, by.x="reg", by.y="REGION_n", all.x=T)
df2 <- df2[order(df2$id), ]

countries_sf$geometry <- NULL

df2 <- merge(df2, countries_sf %>% dplyr::filter(SCENARIO=="SSP5") , by.x="ISO3", by.y="ISO3", all.x=T)
df2 <- df2[order(df2$id), ]

df2$layer.1 <- df2$layer.1 / df2$adj_fact_585
df2$layer.2 <- df2$layer.2 / df2$adj_fact_585
df2$layer.3 <- df2$layer.3 / df2$adj_fact_585

values(r1_g_585[[1]]) <- df2$layer.1 
values(r1_g_585[[2]]) <- df2$layer.2
values(r1_g_585[[3]]) <- df2$layer.3

##############


wrld_simpl$REGION_n <- as.numeric(wrld_simpl$ISO3)
regions_r <- fasterize::fasterize(wrld_simpl, raster(r1_g_126[[1]]), field="REGION_n")
regions <- dplyr::select(wrld_simpl, ISO3, REGION_n)
regions$geometry<-NULL

df3 <- data.frame(values(r1_g_126), reg = values(regions_r))
df3$id  <- 1:nrow(df3)
df3 <- merge(df3, regions, by.x="reg", by.y="REGION_n", all.x=T)
df3 <- df3[order(df3$id), ]

countries_sf$geometry <- NULL

df3 <- merge(df3, countries_sf %>% dplyr::filter(SCENARIO=="SSP1") , by.x="ISO3", by.y="ISO3", all.x=T)
df3 <- df3[order(df3$id), ]

df3$layer.1 <- df3$layer.1 / df3$adj_fact_126
df3$layer.2 <- df3$layer.2 / df3$adj_fact_126
df3$layer.3 <- df3$layer.3 / df3$adj_fact_126

values(r1_g_126[[1]]) <- df3$layer.1 
values(r1_g_126[[2]]) <- df3$layer.2
values(r1_g_126[[3]]) <- df3$layer.3

###############


wrld_simpl$REGION_n <- as.numeric(wrld_simpl$ISO3)
regions_r <- fasterize::fasterize(wrld_simpl, raster(r1_g_370[[1]]), field="REGION_n")
regions <- dplyr::select(wrld_simpl, ISO3, REGION_n)
regions$geometry<-NULL

df4 <- data.frame(values(r1_g_370), reg = values(regions_r))
df4$id  <- 1:nrow(df4)
df4 <- merge(df4, regions, by.x="reg", by.y="REGION_n", all.x=T)
df4 <- df4[order(df4$id), ]

countries_sf$geometry <- NULL

df4 <- merge(df4, countries_sf %>% dplyr::filter(SCENARIO=="SSP3") , by.x="ISO3", by.y="ISO3", all.x=T)
df4 <- df4[order(df4$id), ]

df4$layer.1 <- df4$layer.1 / df4$adj_fact_370
df4$layer.2 <- df4$layer.2 / df4$adj_fact_370
df4$layer.3 <- df4$layer.3 / df4$adj_fact_370

values(r1_g_370[[1]]) <- df4$layer.1 
values(r1_g_370[[2]]) <- df4$layer.2
values(r1_g_370[[3]]) <- df4$layer.3

rm(df1, df2, df3, df4)
gc()

# shift age and gender % in each grid cell (keeping total count constant) using iiasa (country-level) age ratios as a share of total pop.

ssp_stats <- read.csv("SspDb_country_data_2013-06-12.csv")

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1", "SSP2", "SSP3", "SSP5"))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats$group <- ifelse(ssp_stats$VARIABLE %in% unique(ssp_stats$VARIABLE)[c(1,2,4, 11, 22, 23, 25, 32)], "0-19", NA)
ssp_stats$group <- ifelse(ssp_stats$VARIABLE %in% unique(ssp_stats$VARIABLE)[c(3, 5:10, 12, 13, 14, 15, 24, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36)], "19-69", ssp_stats$group)
ssp_stats$group <- ifelse(ssp_stats$VARIABLE %in% unique(ssp_stats$VARIABLE)[c(16:21, 3, 37:42, 24)], "69+", ssp_stats$group)


ssp_stats <- group_by(ssp_stats, SCENARIO, REGION, group) %>% dplyr::summarise_if(., is.numeric, "sum")

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::mutate(X2020=X2020/sum(X2020), X2050=X2050/sum(X2050))

ssp_stats <- dplyr::select(ssp_stats, SCENARIO, REGION, group, X2020, X2050)

ssp_stats <- pivot_wider(ssp_stats, names_from = group, values_from = c(X2020, X2050))

ssp_stats$`diff_0-19` <-  ssp_stats$`X2050_0-19` / ssp_stats$`X2020_0-19`
ssp_stats$`diff_19-69` <- ssp_stats$`X2050_19-69` / ssp_stats$`X2020_19-69`
ssp_stats$`diff_69+` <- ssp_stats$`X2050_69+` / ssp_stats$`X2020_69+`

### ADJUST PERCENTAGE GROWTH RATES IN EACH CELL WITH AVERAGE NATIONAL GROWTH RATE FROM IIASA

wrld_simpl$REGION_n <- as.numeric(wrld_simpl$ISO3)
regions_r <- fasterize::fasterize(wrld_simpl, raster(r1_g_245[[1]]), field="REGION_n")
regions <- dplyr::select(wrld_simpl, ISO3, REGION_n)
regions$geometry<-NULL

df1 <- data.frame(values(r1_g_245), reg = values(regions_r))
df1$share1 <- df1$layer.1/(df1$layer.1 + df1$layer.2 + df1$layer.3)
df1$share2 <- df1$layer.2/(df1$layer.1 + df1$layer.2 + df1$layer.3)
df1$share3 <- df1$layer.3/(df1$layer.1 + df1$layer.2 + df1$layer.3)

df1$id  <- 1:nrow(df1)
df1 <- merge(df1, regions, by.x="reg", by.y="REGION_n", all.x=T)
df1 <- df1[order(df1$id), ]

df1 <- merge(df1,ssp_stats %>% dplyr::filter(SCENARIO=="SSP2") , by.x="ISO3", by.y="REGION", all.x=T)
df1 <- df1[order(df1$id), ]

#

df2 <- data.frame(values(r1_g_585), reg = values(regions_r))
df2$share1 <- df2$layer.1/(df2$layer.1 + df2$layer.2 + df2$layer.3)
df2$share2 <- df2$layer.2/(df2$layer.1 + df2$layer.2 + df2$layer.3)
df2$share3 <- df2$layer.3/(df2$layer.1 + df2$layer.2 + df2$layer.3)

df2$id  <- 1:nrow(df2)
df2 <- merge(df2, regions, by.x="reg", by.y="REGION_n", all.x=T)
df2 <- df2[order(df2$id), ]

df2 <- merge(df2,ssp_stats %>% dplyr::filter(SCENARIO=="SSP5") , by.x="ISO3", by.y="REGION", all.x=T)
df2 <- df2[order(df2$id), ]

#

df3 <- data.frame(values(r1_g_126), reg = values(regions_r))
df3$share1 <- df3$layer.1/(df3$layer.1 + df3$layer.2 + df3$layer.3)
df3$share2 <- df3$layer.2/(df3$layer.1 + df3$layer.2 + df3$layer.3)
df3$share3 <- df3$layer.3/(df3$layer.1 + df3$layer.2 + df3$layer.3)

df3$id  <- 1:nrow(df3)
df3 <- merge(df3, regions, by.x="reg", by.y="REGION_n", all.x=T)
df3 <- df3[order(df3$id), ]

df3 <- merge(df3,ssp_stats %>% dplyr::filter(SCENARIO=="SSP1") , by.x="ISO3", by.y="REGION", all.x=T)
df3 <- df3[order(df3$id), ]

#

df4 <- data.frame(values(r1_g_370), reg = values(regions_r))
df4$share1 <- df4$layer.1/(df4$layer.1 + df4$layer.2 + df4$layer.3)
df4$share2 <- df4$layer.2/(df4$layer.1 + df4$layer.2 + df4$layer.3)
df4$share3 <- df4$layer.3/(df4$layer.1 + df4$layer.2 + df4$layer.3)

df4$id  <- 1:nrow(df4)
df4 <- merge(df4, regions, by.x="reg", by.y="REGION_n", all.x=T)
df4 <- df4[order(df4$id), ]

df4 <- merge(df4,ssp_stats %>% dplyr::filter(SCENARIO=="SSP3") , by.x="ISO3", by.y="REGION", all.x=T)
df4 <- df4[order(df4$id), ]

### RESCALE THE ORIGINAL VALUES BY THE NEW SHARES

df1$share1_2050 <- df1$share1 * df1$`diff_0-19`
df1$share2_2050 <- df1$share2 * df1$`diff_19-69`
df1$share3_2050 <- df1$share3 * df1$`diff_69+`

df1$share1_2050_adj  <- df1$share1_2050  / (df1$share3_2050 + df1$share2_2050 + df1$share1_2050)
df1$share2_2050_adj  <- df1$share2_2050  / (df1$share3_2050 + df1$share2_2050 + df1$share1_2050)
df1$share3_2050_adj  <- df1$share3_2050  / (df1$share3_2050 + df1$share2_2050 + df1$share1_2050)

#

df2$share1_2050 <- df2$share1 * df2$`diff_0-19`
df2$share2_2050 <- df2$share2 * df2$`diff_19-69`
df2$share3_2050 <- df2$share3 * df2$`diff_69+`

df2$share1_2050_adj  <- df2$share1_2050  / (df2$share3_2050 + df2$share2_2050 + df2$share1_2050)
df2$share2_2050_adj  <- df2$share2_2050  / (df2$share3_2050 + df2$share2_2050 + df2$share1_2050)
df2$share3_2050_adj  <- df2$share3_2050  / (df2$share3_2050 + df2$share2_2050 + df2$share1_2050)

df3$share1_2050 <- df3$share1 * df3$`diff_0-19`
df3$share2_2050 <- df3$share2 * df3$`diff_19-69`
df3$share3_2050 <- df3$share3 * df3$`diff_69+`

df3$share1_2050_adj  <- df3$share1_2050  / (df3$share3_2050 + df3$share2_2050 + df3$share1_2050)
df3$share2_2050_adj  <- df3$share2_2050  / (df3$share3_2050 + df3$share2_2050 + df3$share1_2050)
df3$share3_2050_adj  <- df3$share3_2050  / (df3$share3_2050 + df3$share2_2050 + df3$share1_2050)

df4$share1_2050 <- df4$share1 * df4$`diff_0-19`
df4$share2_2050 <- df4$share2 * df4$`diff_19-69`
df4$share3_2050 <- df4$share3 * df4$`diff_69+`

df4$share1_2050_adj  <- df4$share1_2050  / (df4$share3_2050 + df4$share2_2050 + df4$share1_2050)
df4$share2_2050_adj  <- df4$share2_2050  / (df4$share3_2050 + df4$share2_2050 + df4$share1_2050)
df4$share3_2050_adj  <- df4$share3_2050  / (df4$share3_2050 + df4$share2_2050 + df4$share1_2050)

### BRING BACK EVERYTHING TO RASTERS

df1$r1_g_245_20501_adj <- (df1$layer.1  + df1$layer.2 + df1$layer.3) * df1$share1_2050_adj
df1$r1_g_245_20502_adj <- (df1$layer.1  + df1$layer.2 + df1$layer.3) * df1$share2_2050_adj
df1$r1_g_245_20503_adj <- (df1$layer.1  + df1$layer.2 + df1$layer.3) * df1$share3_2050_adj

values(r1_g_245[[1]]) <- df1$r1_g_245_20501_adj 
values(r1_g_245[[2]]) <- df1$r1_g_245_20502_adj
values(r1_g_245[[3]]) <- df1$r1_g_245_20503_adj

#

df2$r1_g_585_20501_adj <- (df2$layer.1  + df2$layer.2 + df2$layer.3) * df2$share1_2050_adj
df2$r1_g_585_20502_adj <- (df2$layer.1  + df2$layer.2 + df2$layer.3) * df2$share2_2050_adj
df2$r1_g_585_20503_adj <- (df2$layer.1  + df2$layer.2 + df2$layer.3) * df2$share3_2050_adj

values(r1_g_585[[1]]) <- df2$r1_g_585_20501_adj 
values(r1_g_585[[2]]) <- df2$r1_g_585_20502_adj
values(r1_g_585[[3]]) <- df2$r1_g_585_20503_adj

#

df3$r1_g_126_20501_adj <- (df3$layer.1  + df3$layer.2 + df3$layer.3) * df3$share1_2050_adj
df3$r1_g_126_20502_adj <- (df3$layer.1  + df3$layer.2 + df3$layer.3) * df3$share2_2050_adj
df3$r1_g_126_20503_adj <- (df3$layer.1  + df3$layer.2 + df3$layer.3) * df3$share3_2050_adj

values(r1_g_126[[1]]) <- df3$r1_g_126_20501_adj 
values(r1_g_126[[2]]) <- df3$r1_g_126_20502_adj
values(r1_g_126[[3]]) <- df3$r1_g_126_20503_adj

#

df4$r1_g_370_20501_adj <- (df4$layer.1  + df4$layer.2 + df4$layer.3) * df4$share1_2050_adj
df4$r1_g_370_20502_adj <- (df4$layer.1  + df4$layer.2 + df4$layer.3) * df4$share2_2050_adj
df4$r1_g_370_20503_adj <- (df4$layer.1  + df4$layer.2 + df4$layer.3) * df4$share3_2050_adj

values(r1_g_370[[1]]) <- df4$r1_g_370_20501_adj 
values(r1_g_370[[2]]) <- df4$r1_g_370_20502_adj
values(r1_g_370[[3]]) <- df4$r1_g_370_20503_adj

#
# check consistency

# wrld_simpl <- wrld_simpl %>% filter(ISO3=="ITA")
# wrld_simpl <- fasterize::fasterize(wrld_simpl, r1_g_585[[1]])
# r1_g_585 <- mask(r1_g_585, wrld_simpl)
# 
# gh <- as.data.frame(r1_g_585)
# 
# sum(gh$layer.1, na.rm=T) / (sum(gh$layer.1, na.rm=T) + sum(gh$layer.2, na.rm=T) + sum(gh$layer.3, na.rm=T))
# sum(gh$layer.2, na.rm=T) / (sum(gh$layer.1, na.rm=T) + sum(gh$layer.2, na.rm=T) + sum(gh$layer.3, na.rm=T))
# sum(gh$layer.3, na.rm=T) / (sum(gh$layer.1, na.rm=T) + sum(gh$layer.2, na.rm=T) + sum(gh$layer.3, na.rm=T))
# 
# ssp_stats %>% filter(REGION=="ITA" & SCENARIO == "SSP5") 

##################################################################

terra::writeRaster(rast(r1_g), paste0(getwd(), "/r1_2020.tif"), overwrite=T)
terra::writeRaster(rast(r1_g_245), paste0(getwd(), "/r1_g_245_2050_final.tif"), overwrite=T)
terra::writeRaster(rast(r1_g_585), paste0(getwd(), "/r1_g_585_2050_final.tif"), overwrite=T)
terra::writeRaster(rast(r1_g_126), paste0(getwd(), "/r1_g_126_2050_final.tif"), overwrite=T)
terra::writeRaster(rast(r1_g_370), paste0(getwd(), "/r1_g_370_2050_final.tif"), overwrite=T)
