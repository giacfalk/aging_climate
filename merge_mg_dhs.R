list <- list.files("F:/DHS/IRHR")

list <- list[grep("IRHR", list)]

df <- data.frame(name=list, stringsAsFactors = F)
df$cc <- substr(df$name, 1, 2)
df$wave <- substr(df$name, 7, 8)
df$wave <- gsub("[^0-9.-]", "", df$wave)

library(tidyverse)
#df <- group_by(df, cc) %>% filter(wave==max(wave))

list <- df$name

library(haven)

setwd("F:/DHS")

liston <- list()

list2 <- list.files("F:/DHS", recursive = T, pattern = ".shp")
list2 <- list2[-grep("xml", list2)]

jj <- gsub(".shp", "", sub('.*\\/', '', list2))
jj <- gsub("FL", "s.dta", jj)
jj <- gsub("GE", "IRHR", jj)
jj <- paste0(substr(jj, 1, 7), "s.dta")

list <-list[list %in% jj == TRUE]

 
for(i in 1:length(list)){

  print(i)
  
s <- read_stata(paste0("F:/DHS/IRHR/", list[i]))
s$year <- s$v007

#years <- as.numeric(substr(sub("^.*?/", "", list2), 4,7))

library(sf)

j <- read_sf(list2[match(list[i], jj)])

liston[[i]] <- merge(s, j, by.x="v001", by.y="DHSCLUST")
}

library(nngeo)
mg <- read.csv("D:/OneDrive - IIASA/MCC/NDVI/mg_with_evi_110121.csv") %>% st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>% st_transform(3395)

for(i in 1:length(list)){
  print(i)
  liston[[i]] <- liston[[i]] %>%  filter(duplicated(geometry) == FALSE)
  liston[[i]] <- st_transform((liston[[i]] %>% st_as_sf()), 3395)
}

liston2 <- list()

for(i in 1:length(list)){
  print(i)
  liston2[[i]] <- st_join(mg, (liston[[i]] %>% st_as_sf()), st_nn, k = 1, maxdist = 10000, left=F)
  }

#liston_m <- do.call(rbind, liston2) %>% st_as_sf()

library(dplyr)
for(i in 1:length(list)){
  print(i)
  liston2[[i]] <- liston2[[i]] %>% 
  group_by(DHSCC, hhid) %>%
  filter(n()>1) %>% 
  ungroup()
}

#colnames(liston_m2)[colnames(liston_m2)==id] <- "minigrid_if"

#write_sf(liston_m2, "D:/OneDrive - IIASA/MCC/DHSmerged_dhs_all_couplemissing.gpkg")

for(i in 1:length(list)){
  liston2[[i]]$geometry=NULL
}

library(purrr)
liston3 <- keep(liston2, ~ nrow(.x) > 0)

for(i in 1:length(liston3)){
write.csv(liston3[[i]], paste0(i, ".csv"))
}
