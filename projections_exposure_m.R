# 0) path and libraries
library(raster)
library(sf)
library(tidyverse)
library(rasterVis)
require(rgdal)
require(maptools)
library(pbapply)
library(terra)
library(knitr)
library(kableExtra)
library(modelsummary)
library(openxlsx)
library(xtable)
library(ggforce)
library(maptools)
library(weights)
library(spatstat)
library(rworldmap)
library(scales)
library(patchwork)
library(stars)
library(viridis)
library(devtools)
library(stargazer)
library(readxl)
library(nominatimlite)
library(urbnmapr)

# set the working directory
stub <- "./"
setwd(stub)

# set the raster elaborating temporary directory
rasterOptions(tmpdir = "./")

clip_raster_with_sf <- function(rast, sf){
  
  r_mask <- fasterize::fasterize(sf, rast)
  r_masked <- raster::mask(rast, r_mask)
  
  return(r_masked)
  
}

###################

source("project_pop.R")

r1_g <- stack("r1_2020.tif")
r1_ssp2_2050 <- stack("r1_g_245_2050_final.tif")
r1_ssp5_2050 <- stack("r1_g_585_2050_final.tif")
r1_ssp1_2050 <- stack("r1_g_126_2050_final.tif")
r1_ssp3_2050 <- stack("r1_g_370_2050_final.tif")
#
cdds <- list.files(pattern="cdd24", full.names = T)
cdds <- cdds[-grep("ensemble", cdds)]
cdds <- lapply(cdds, function(X){stack(rast(X))})
#

c1 <- mean(cdds[[5]], na.rm=T)
c2 <- mean(cdds[[2]], na.rm=T)
c3 <- mean(cdds[[4]], na.rm=T)
c4 <- mean(cdds[[1]], na.rm=T)
c5 <- mean(cdds[[3]], na.rm=T)

#
c1 <- projectRaster(c1, r1_g[[1]])
c2 <- projectRaster(c2, r1_g[[1]])
c3 <- projectRaster(c3, r1_g[[1]])
c4 <- projectRaster(c4, r1_g[[1]])
c5 <- projectRaster(c5, r1_g[[1]])

#
c1_view <- c1
c1_view[c1_view==0] <- NA
#
c2_view <- c2
c2_view[c2_view==0] <- NA
#
c3_view <- c3
c3_view[c3_view==0] <- NA
#
c4_view <- c4
c4_view[c4_view==0] <- NA
#
c5_view <- c5
c5_view[c5_view==0] <- NA
#
#plot historical CDDs vs current share of above 60 in every grid cell
#
#
r1_tot <- sum(r1_g, na.rm=T)
r1_ssp2_2050_tot <- sum(r1_ssp2_2050, na.rm=T)
r1_ssp5_2050_tot <- sum(r1_ssp5_2050, na.rm=T)
r1_ssp1_2050_tot <- sum(r1_ssp1_2050, na.rm=T)
r1_ssp3_2050_tot <- sum(r1_ssp3_2050, na.rm=T)
#
#
#
#
r1_view <- (r1_g[[3]]) / r1_tot
r1_view[r1_view==0 | is.infinite(r1_view)] <- NA
#
r1_ssp2_2050[[3]][r1_ssp2_2050[[3]]<10] <- NA
r2_view <- (r1_ssp2_2050[[3]]) / r1_ssp2_2050_tot
r2_view[r2_view==0 | is.infinite(r2_view)] <- NA
#
r1_ssp5_2050[[3]][r1_ssp5_2050[[3]]<10] <- NA
r3_view <- (r1_ssp5_2050[[3]]) / r1_ssp5_2050_tot
r3_view[r3_view==0 | is.infinite(r3_view)] <- NA
#
r1_ssp1_2050[[3]][r1_ssp1_2050[[3]]<10] <- NA
r4_view <- (r1_ssp1_2050[[3]]) / r1_ssp1_2050_tot
r4_view[r4_view==0 | is.infinite(r4_view)] <- NA
#
r1_ssp3_2050[[3]][r1_ssp3_2050[[3]]<10] <- NA
r5_view <- (r1_ssp3_2050[[3]]) / r1_ssp3_2050_tot
r5_view[r5_view==0 | is.infinite(r5_view)] <- NA
#

r1_view[is.na(r3_view)] <- NA

#
#
#
#
tasmaxs <- list.files(pattern="tasmax|tmax", full.names = T)
tasmaxs <- tasmaxs[-grep("ensemble", tasmaxs)]
tasmaxs <- lapply(tasmaxs, function(X){stack(rast(X))})
#
c1_tasmax <- stackApply(tasmaxs[[5]], 1, median, na.rm=T)
c2_tasmax <- stackApply(tasmaxs[[2]], 1, median, na.rm=T)
c3_tasmax <- stackApply(tasmaxs[[4]], 1, median, na.rm=T)
c4_tasmax <- stackApply(tasmaxs[[1]], 1, median, na.rm=T)
c5_tasmax <- stackApply(tasmaxs[[3]], 1, median, na.rm=T)
#
c1_tasmax <- projectRaster(c1_tasmax, r1_g[[1]])
c2_tasmax <- projectRaster(c2_tasmax, r1_g[[1]])
c3_tasmax <- projectRaster(c3_tasmax, r1_g[[1]])
c4_tasmax <- projectRaster(c4_tasmax, r1_g[[1]])
c5_tasmax <- projectRaster(c5_tasmax, r1_g[[1]])
#
c1_tasmax_view <- c1_tasmax
c2_tasmax_view <- c2_tasmax
c3_tasmax_view <- c3_tasmax
c4_tasmax_view <- c4_tasmax
c5_tasmax_view <- c5_tasmax
#
#######
#
#multi-GCM data
#
multigcm_hist <- stack(rast("cdd24_global_hist_all_models.nc"))
multigcm_245 <- stack(rast("cdd24_global_245_all_models.nc"))
multigcm_585 <- stack(rast("cdd24_global_585_all_models.nc"))
multigcm_126 <- stack(rast("cdd24_global_126_all_models.nc"))
multigcm_370 <-  stack(rast("cdd24_global_370_all_models.nc"))
#
multigcm_hist <- projectRaster(multigcm_hist, r1_g[[1]])
multigcm_245 <- projectRaster(multigcm_245, r1_g[[1]])
multigcm_585 <- projectRaster(multigcm_585, r1_g[[1]])
multigcm_126 <- projectRaster(multigcm_126, r1_g[[1]])
multigcm_370 <- projectRaster(multigcm_370, r1_g[[1]])

names(multigcm_hist) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")
names(multigcm_126) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0")
names(multigcm_245) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")
names(multigcm_370) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0")
names(multigcm_585) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")

#
#
#
multigcm_tasmax_hist <- stack(rast("tasmax95_global_hist_all_models.nc"))
multigcm_tasmax_245 <- stack(rast("tasmax95_global_245_all_models.nc"))
multigcm_tasmax_585 <- stack(rast("tasmax95_global_585_all_models.nc"))
multigcm_tasmax_126 <- stack(rast("tasmax95_global_126_all_models.nc"))
multigcm_tasmax_370 <- stack(rast("tasmax95_global_370_all_models.nc"))
#
multigcm_tasmax_hist <- projectRaster(multigcm_tasmax_hist, r1_g[[1]])
multigcm_tasmax_245 <- projectRaster(multigcm_tasmax_245, r1_g[[1]])
multigcm_tasmax_585 <- projectRaster(multigcm_tasmax_585, r1_g[[1]])
multigcm_tasmax_126 <- projectRaster(multigcm_tasmax_126, r1_g[[1]])
multigcm_tasmax_370 <- projectRaster(multigcm_tasmax_370, r1_g[[1]])

names(multigcm_tasmax_hist) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")
names(multigcm_tasmax_126) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
names(multigcm_tasmax_245) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")
names(multigcm_tasmax_370) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
names(multigcm_tasmax_585) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")

#
###
#
multigcm_hotday_hist <- stack(rast("hotdays_global_hist_all_models.nc"))
multigcm_hotday_245 <- stack(rast("hotdays_global_245_all_models.nc"))
multigcm_hotday_585 <- stack(rast("hotdays_global_585_all_models.nc"))
multigcm_hotday_126 <- stack(rast("hotdays_global_126_all_models.nc"))
multigcm_hotday_370 <- stack(rast("hotdays_global_370_all_models.nc"))

#
multigcm_hotday_hist <- projectRaster(multigcm_hotday_hist, r1_g[[1]])
multigcm_hotday_245 <- projectRaster(multigcm_hotday_245, r1_g[[1]])
multigcm_hotday_585 <- projectRaster(multigcm_hotday_585, r1_g[[1]])
multigcm_hotday_126 <- projectRaster(multigcm_hotday_126, r1_g[[1]])
multigcm_hotday_370 <- projectRaster(multigcm_hotday_370, r1_g[[1]])
#
names(multigcm_hotday_hist) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")
names(multigcm_hotday_126) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
names(multigcm_hotday_245) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")
names(multigcm_hotday_370) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
names(multigcm_hotday_585) <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2", "FGOALS-g3"," GFDL-CM4", "GFDL-ESM4", "GISS-E2-1-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM")

#
#
##########
#
hotdays <- list.files(pattern="hotdays", full.names = T)
hotdays <- hotdays[-grep("ensemble", hotdays)]
hotdays <- lapply(hotdays, function(X){stack(rast(X))})
#
c1_hotday <- mean(hotdays[[5]], na.rm=T)
c2_hotday <- mean(hotdays[[2]], na.rm=T)
c3_hotday <- mean(hotdays[[4]], na.rm=T)
c4_hotday <- mean(hotdays[[1]], na.rm=T)
c5_hotday <- mean(hotdays[[3]], na.rm=T)
#
c1_hotday_view <- c1_hotday
c2_hotday_view <- c2_hotday
c3_hotday_view <- c3_hotday
c4_hotday_view <- c4_hotday
c5_hotday_view <- c5_hotday
#
c1_hotday <- projectRaster(c1_hotday, r1_g[[1]])
c2_hotday <- projectRaster(c2_hotday, r1_g[[1]])
c3_hotday <- projectRaster(c3_hotday, r1_g[[1]])
c4_hotday <- projectRaster(c4_hotday, r1_g[[1]])
c5_hotday <- projectRaster(c5_hotday, r1_g[[1]])

save.image("data_ready_v2.Rdata")

#

# load("data_ready_v2.Rdata")

###################
### Figure 1 

data(wrld_simpl)

line = .8
cex = 1
side = 3
adj=-0.085

wrld_simpl <- subset(wrld_simpl, wrld_simpl$ISO3!="ATA")

setwd(stub)

# logarithmic palettes
# https://stackoverflow.com/questions/38627398/r-image-intensities-in-the-log-scale

# r2_view_d <- r2_view - r1_view
# r3_view_d <- r3_view - r1_view
# 
# c2_tasmax_view_d <- c2_tasmax_view - c1_tasmax_view
# c3_tasmax_view_d <- c3_tasmax_view - c1_tasmax_view
# 
# c2_view_d <- c2_view - c1_view
# c3_view_d <- c3_view - c1_view
# 
# c2_hotday_view_d <- c2_hotday_view - c1_hotday_view
# c3_hotday_view_d <- c3_hotday_view - c1_hotday_view

# r1_view_d <- r1_view
# r1_view_d[r1_view_d>0.4] <- 0.4
# 
# r1_view_d <- r1_view_d*100
# 
# r2_view_d[r2_view_d>0.4] <- 0.4
# r2_view_d <- r2_view_d*100
# 
# r3_view_d[r3_view_d>0.4] <- 0.4
# r3_view_d <- r3_view_d*100


# #
# 
# breaks <- c(0, 0.025, 0.05, 0.1, 0.2, 0.4)*100
# col_breaks <- rev(magma(length(breaks)))
# 
# breaks_cdd <- c(30, 350, 700, 1500, 5000)
# col_breaks_cdd <- rev(heat.colors(length(breaks_cdd)))
# 
# breaks_tmax <- c(10, 20, 30, 35, 40, 45, 50)
# col_breaks_tmax <- rev(heat.colors(length(breaks_tmax)))
# 
# #
# 
# breaks_delta <- c(-0.001, 0, 0.025, 0.05, 0.1, 0.2, 0.4)*100
# col_breaks_delta <- rev(magma(length(breaks_delta)))
# 
# breaks_cdd_delta <- c(30, 60, 150, 300, 1000)
# col_breaks_cdd_delta <- rev(heat.colors(length(breaks_cdd_delta)))
# 
# breaks_tmax_delta <- c(0.5, 1, 1.5, 2, 2.5)
# col_breaks_tmax_delta <- rev(heat.colors(length(breaks_tmax_delta)))
# 
# #
# 
# pdf("maps.pdf", height = 5*1.2, width = 8*1.2)
# par(mfrow = c(3, 3), mai=c(.01, .01, .01, .01),
#     oma = c(1,1,0,0) + 0.01,
#     mar = c(0,0,0,0) + 0.01)
# 
# plot(r1_view_d, legend = FALSE, axes=FALSE, box=FALSE ,col = col_breaks, breaks= breaks, ylim=range(-60:90), cex.main=1.1)
# title("Population % aged 69+, current", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# legend("bottom", legend = breaks, fill = col_breaks, border = FALSE, bty = "n", horiz=T)
# #mtext("A", side=side, line=line, cex=cex, adj=adj)
# 
# plot(c1_view, legend = FALSE, axes=FALSE, box=FALSE, col = col_breaks_cdd, breaks= breaks_cdd, ylim=range(-60:90), cex.main=1.1 )
# title("Historical CDDs (1995-2014 avg.)", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# legend("bottom", legend = breaks_cdd, fill = col_breaks_cdd, border = FALSE, bty = "n", horiz=T)
# #mtext("B", side=side, line=line, cex=cex, adj=adj)
# 
# plot(c1_tasmax, legend = FALSE, axes=FALSE, box=FALSE ,col = col_breaks_tmax, breaks= breaks_tmax, ylim=range(-60:90), cex.main=1.1 )
# title("Hist. 95th pctl. Tmax (1995-2014 avg.)", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# legend("bottom", legend = breaks_tmax, fill = col_breaks_tmax, border = FALSE, bty = "n", horiz=T)
# #mtext("C", side=side, line=line, cex=cex, adj=adj)
# 
# plot(r2_view_d, legend = FALSE, axes=FALSE, box=FALSE, col = col_breaks_delta, breaks_delta= breaks_delta, ylim=range(-60:90), cex.main=1.1 )
# title("Pct. pts. change in pop. aged 69+, 2050, SSP245", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# #mtext("D", side=side, line=line, cex=cex, adj=adj)
# 
# plot(c2_view, legend = FALSE, axes=FALSE, box=FALSE, col = col_breaks_cdd_delta, breaks= breaks_cdd_delta, ylim=range(-60:90), cex.main=1.1 )
# title("Change in CDDs (SSP245, 2050)", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# #mtext("E", side=side, line=line, cex=cex, adj=adj)
# 
# plot(c2_tasmax_view, legend = FALSE, axes=FALSE, box=FALSE ,col = col_breaks_tmax_delta, breaks= breaks_tmax_delta, ylim=range(-60:90), cex.main=1.1 )
# title("Change in 95th pctl. Tmax (SSP245, 2050)", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# #mtext("F", side=side, line=line, cex=cex, adj=adj)
# 
# plot(r3_view_d, legend = FALSE, axes=FALSE, box=FALSE, col = col_breaks_delta, breaks= breaks_delta, ylim=range(-60:90), cex.main=1.1 )
# title("Pct. pts. change in pop. aged 69+, 2050, SSP585", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# legend("bottom", legend = breaks_delta, fill = col_breaks_delta, border = FALSE, bty = "n", horiz=T)
# #mtext("G", side=side, line=line, cex=cex, adj=adj)
# 
# plot(c3_view, legend = FALSE, axes=FALSE, box=FALSE, col = col_breaks_cdd_delta, breaks= breaks_cdd_delta, ylim=range(-60:90), cex.main=1.1 )
# title("Change in CDDs (SSP585, 2050)", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# legend("bottom", legend = breaks_cdd_delta, fill = col_breaks_cdd_delta, border = FALSE, bty = "n", horiz=T)
# #mtext("H", side=side, line=line, cex=cex, adj=adj)
# 
# plot(c3_tasmax_view, legend = FALSE, axes=FALSE, box=FALSE ,col = col_breaks_tmax_delta, breaks= breaks_tmax_delta, ylim=range(-60:90), cex.main=1.1 )
# title("Change in 95th pctl. Tmax  (SSP585, 2050)", line = -2, cex.main=1.1)
# par(xpd=TRUE)
# plot(wrld_simpl, add=TRUE, fill=NA, lwd=0.005, border="grey")
# legend("bottom", legend = breaks_tmax_delta, fill = col_breaks_tmax_delta, border = FALSE, bty = "n", horiz=T)
# #mtext("I", side=side, line=line, cex=cex, adj=adj)
# 
# dev.off()

r2_view_d <- r2_view 
r3_view_d <- r3_view 

c2_tasmax_view_d <- c2_tasmax_view
c3_tasmax_view_d <- c3_tasmax_view

c2_view_d <- c2_view
c3_view_d <- c3_view

c2_hotday_view_d <- c2_hotday_view
c3_hotday_view_d <- c3_hotday_view 


#####

f1 <- function(X){
  
  qq <- round(quantile(X, c(0.05, 0.3, 0.6, 0.95)), 2)
  return(qq)
}

f2 <- function(X,Y){
  
  cut(X, breaks=c(-Inf, Y, Inf), labels=c(paste0("< ", Y[1]), paste0(Y[1], " - ", Y[2]), paste0(Y[2], " - ", Y[3]), paste0(Y[3], " - ", Y[4]), paste0("> ", Y[4])))
  
}

f2_pct <- function(X,Y){
  
  cut(X, breaks=c(-Inf, Y, Inf), labels=c(paste0("< ", Y[1]*100), paste0(Y[1]*100, " - ", Y[2]*100), paste0(Y[2]*100, " - ", Y[3]*100), paste0(Y[3]*100, " - ", Y[4]*100), paste0("> ", Y[4]*100)))
  
}


rr = na.omit(values(r3_view_d))

r1_view_d <- st_as_sf(st_as_stars(aggregate(r1_view, 6, "mean")))
r1_view_d$layer_d <- f2_pct(r1_view_d$layer, f1(rr))

r2_view_d <- st_as_sf(st_as_stars(aggregate(r2_view_d, 6, "mean")))
r2_view_d$layer_d <- f2_pct(r2_view_d$layer, f1(rr))

r3_view_d <- st_as_sf(st_as_stars(aggregate(r3_view_d, 6, "mean")))
r3_view_d$layer_d <- f2_pct(r3_view_d$layer, f1(rr))

rr = na.omit(values(c3_tasmax_view_d))

c1_tasmax_view_d <- st_as_sf(st_as_stars(aggregate(c1_tasmax_view, 6, "mean")))
c1_tasmax_view_d$index_1 = ifelse(c1_tasmax_view_d$index_1==0, NA, c1_tasmax_view_d$index_1)
c1_tasmax_view_d$layer_d <- f2(c1_tasmax_view_d$index_1, f1(rr))

c2_tasmax_view_d <- st_as_sf(st_as_stars(aggregate(c2_tasmax_view_d, 6, "mean")))
c2_tasmax_view_d$index_1 = ifelse(c2_tasmax_view_d$index_1==0, NA, c2_tasmax_view_d$index_1)
c2_tasmax_view_d$layer_d <- f2(c2_tasmax_view_d$index_1, f1(rr))

c3_tasmax_view_d <- st_as_sf(st_as_stars(aggregate(c3_tasmax_view_d, 6, "mean")))
c3_tasmax_view_d$index_1 = ifelse(c3_tasmax_view_d$index_1==0, NA, c3_tasmax_view_d$index_1)
c3_tasmax_view_d$layer_d <- f2(c3_tasmax_view_d$index_1, f1(rr))

rr = na.omit(values(c3_view_d))

c1_view_d <- st_as_sf(st_as_stars(aggregate(c1_view, 6, "mean")))
c1_view_d$layer = ifelse(c1_view_d$layer==0, NA, c1_view_d$layer)
c1_view_d$layer_d <- f2(c1_view_d$layer, f1(rr))

c2_view_d <- st_as_sf(st_as_stars(aggregate(c2_view_d, 6, "mean")))
c2_view_d$layer = ifelse(c2_view_d$layer==0, NA, c2_view_d$layer)
c2_view_d$layer_d <- f2(c2_view_d$layer, f1(rr))

c3_view_d <- st_as_sf(st_as_stars(aggregate(c3_view_d, 6, "mean")))
c3_view_d$layer = ifelse(c3_view_d$layer==0, NA, c3_view_d$layer)
c3_view_d$layer_d <- f2(c3_view_d$layer, f1(rr))

rr = na.omit(values(c3_hotday_view_d))

c1_hotday_view_d <- st_as_sf(st_as_stars(aggregate(c1_hotday_view, 6, "mean")))
c1_hotday_view_d$layer = ifelse(c1_hotday_view_d$layer==0, NA, c1_hotday_view_d$layer)
c1_hotday_view_d$layer_d <- f2(c1_hotday_view_d$layer, f1(rr[rr>0]))

c2_hotday_view_d <- st_as_sf(st_as_stars(aggregate(c2_hotday_view_d, 6, "mean")))
c2_hotday_view_d$layer = ifelse(c2_hotday_view_d$layer==0, NA, c2_hotday_view_d$layer)
c2_hotday_view_d$layer_d <- f2(c2_hotday_view_d$layer, f1(rr[rr>0]))

c3_hotday_view_d <- st_as_sf(st_as_stars(aggregate(c3_hotday_view_d, 6, "mean")))
c3_hotday_view_d$layer = ifelse(c3_hotday_view_d$layer==0, NA, c3_hotday_view_d$layer)
c3_hotday_view_d$layer_d <- f2(c3_hotday_view_d$layer, f1(rr[rr>0]))

####

fig1a <-  ggplot() + geom_sf(data = r1_view_d
                             ,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl)
          , fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig1b <-  ggplot() + geom_sf(data = r2_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "bottom", legend.direction = "horizontal",  legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

fig1c <-  ggplot() + geom_sf(data = r3_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.direction = "horizontal", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig2a <-  ggplot() + geom_sf(data = c1_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig2b <-  ggplot() + geom_sf(data = c2_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "bottom",legend.direction = "horizontal",  legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

fig2c <-  ggplot() + geom_sf(data = c3_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.direction = "horizontal", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig3a <-  ggplot() + geom_sf(data = c1_tasmax_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig3b <-  ggplot() + geom_sf(data = c2_tasmax_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "bottom",legend.direction = "horizontal",  legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

fig3c <-  ggplot() + geom_sf(data = c3_tasmax_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.direction = "horizontal", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig4a <-  ggplot() + geom_sf(data = c1_hotday_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")

fig4b <-  ggplot() + geom_sf(data = c2_hotday_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "bottom",legend.direction = "horizontal",  legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

fig4c <-  ggplot() + geom_sf(data = c3_hotday_view_d,
                             aes(fill = layer_d), colour=NA)+
  scale_fill_discrete(type = heat.colors(5, rev = T), name="", na.value = NA)+
  geom_sf(data=st_as_sf(wrld_simpl), fill=NA, size =.005)+
  theme_classic()+
  #coord_sf(ylim=c(-34, 36), xlim=c(-20, 51))+
  theme(legend.position = "none", legend.direction = "horizontal", legend.text=element_text(size=8), axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(size = 10, hjust = 0.5), aspect.ratio=0.55)+

  xlab("")+
  ylab("")


top_row <- (fig1a + fig1b + fig1c) + plot_annotation(tag_levels = list("A", "B", "C"))

row2 <- (fig2a + fig2b + fig2c) + plot_annotation(tag_levels = list(c("D", "E", "F")))

row3 <- (fig3a + fig3b + fig3c) + plot_annotation(tag_levels = list(c("G", "H", "I")))

bottom_row <- (fig4a + fig4b + fig4c) + plot_annotation(tag_levels = list(c("J", "K", "L")))

f1 <- cowplot::plot_grid(top_row, row2, row3, bottom_row, nrow=4)

ggsave("maps.pdf", f1, scale=3.5, width = 4, height = 4)

### bivariate map(s)

source("biv_funcs.R")

rhist <- stack(r1_view, c1)
names(rhist) <- c("Aged 69+ (%)", "CDDs")

r245 <- stack(r2_view, c2)
names(r245) <- c("Aged 69+ (%)", "CDDs")

r585 <- stack(r3_view, c3)
names(r585) <- c("Aged 69+ (%)", "CDDs")

r126 <- stack(r4_view, c4)
names(r126) <- c("Aged 69+ (%)", "CDDs")

r370 <- stack(r5_view, c5)
names(r370) <- c("Aged 69+ (%)", "CDDs")

quantile(na.omit(values(rhist[[1]])), seq(0, 1, 0.25))
quantile(na.omit(values(rhist[[2]])), seq(0, 1, 0.2))

regions <- st_as_sf(rworldmap::countriesLow)
sf_use_s2(F)
regions_globe = dplyr::summarise(regions)
sf_use_s2(T)
regions = filter(regions, REGION !="Antarctica" & !is.na(REGION))
regions$REGION_n <- as.numeric(regions$REGION)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="REGION_n")
regions <- dplyr::select(regions, REGION, REGION_n)

########

for (scen in c("rhist", "r126", "r245", "r370", "r585")){
  
  scenario = ifelse(scen=="rhist", "2020", ifelse(scen=="r126", "2020-2050 (SSP126)", ifelse(scen=="r245", "2020-2050 (SSP245)", ifelse(scen=="r370", "2020-2050 (SSP370)",  "2020-2050 (SSP585)"))))
  r = get(scen)
  r = raster::crop(r, extent(regions_globe))
  
  for (i in 1:2){
    r[[i]] = clip_raster_with_sf(r[[i]], regions_globe)
  }
  
  
  region = "global"
  
  # create the bivariate raster
  bivmap <- bivariate.map(rasterx = r[["Aged.69....."]], rastery = r[["CDDs"]],
                          export.colour.matrix = TRUE, outname = "bivMapCols",
                          colormatrix = col.matrix, nquantiles = nBreaks)
  
  bivMapDF <- as.data.frame(bivmap, xy = TRUE) %>%
    tbl_df() %>%
    dplyr::rename("BivValue" = 3) %>%
    pivot_longer(., names_to = "Variable", values_to = "bivVal", cols = BivValue)
  
  # Make the map using ggplot
  map <- ggplot(data=bivMapDF) +
    ggtitle(paste0("Share of individuals aged 69+ and cooling degree days, ", scenario, ", ", region))+
    geom_raster(data=bivMapDF, aes(x = x, y = y, fill = bivVal)) +
    scale_fill_gradientn(colours = col.matrix, na.value = "transparent")+
    theme_classic()+
    theme(legend.position = "none", axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    guides(fill = guide_colorsteps(barwidth = 1, barheight = 5))+
    xlab("")+
    ylab("")+
    geom_sf(data=regions, fill=NA, size =.3)
  
  ggsave(plot = map,
         filename = paste0("map_biv_cdd_", region, "_", scen, ".pdf"),
         device = "pdf", path = getwd(),
         dpi = 320, scale=1)
  
}

###

source("biv_funcs.R")

rhist <- stack(r1_view, c1_tasmax)
names(rhist) <- c("Aged 69+ (%)", "TMAX95")

r245 <- stack(r2_view, c2_tasmax)
names(r245) <- c("Aged 69+ (%)", "TMAX95")

r585 <- stack(r3_view, c3_tasmax)
names(r585) <- c("Aged 69+ (%)", "TMAX95")

r126 <- stack(r4_view, c4_tasmax)
names(r126) <- c("Aged 69+ (%)", "TMAX95")

r370 <- stack(r5_view, c5_tasmax)
names(r370) <- c("Aged 69+ (%)", "TMAX95")


for (scen in c("rhist", "r126", "r245", "r370", "r585")){
  
  scenario = ifelse(scen=="rhist", "2020", ifelse(scen=="r126", "2020-2050 (SSP126)", ifelse(scen=="r245", "2020-2050 (SSP245)", ifelse(scen=="r370", "2020-2050 (SSP370)",  "2020-2050 (SSP585)"))))
  r = get(scen)
  r = raster::crop(r, extent(regions_globe))
  
  for (i in 1:2){
    r[[i]] = clip_raster_with_sf(r[[i]], regions_globe)
  }
  
  
  region = "global"
  
  # create the bivariate raster
  bivmap <- bivariate.map(rasterx = r[["Aged.69....."]], rastery = r[["TMAX95"]],
                          export.colour.matrix = TRUE, outname = "bivMapCols",
                          colormatrix = col.matrix, nquantiles = nBreaks)
  
  bivMapDF <- as.data.frame(bivmap, xy = TRUE) %>%
    tbl_df() %>%
    dplyr::rename("BivValue" = 3) %>%
    pivot_longer(., names_to = "Variable", values_to = "bivVal", cols = BivValue)
  
  # Make the map using ggplot
  map <- ggplot(data=bivMapDF) +
    ggtitle(paste0("Share of individuals aged 69+ and 95th pctl. maximum temperature, ", scenario, ", ", region))+
    geom_raster(data=bivMapDF, aes(x = x, y = y, fill = bivVal)) +
    scale_fill_gradientn(colours = col.matrix, na.value = "transparent")+
    theme_classic()+
    theme(legend.position = "none", axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    guides(fill = guide_colorsteps(barwidth = 1, barheight = 5))+
    xlab("")+
    ylab("")+
    geom_sf(data=regions, fill=NA, size =.3)
  
  ggsave(plot = map,
         filename = paste0("map_biv_tmax_", region, "_", scen, ".pdf"),
         device = "pdf", path = getwd(),
         dpi = 320, scale=1)
  
}

source("biv_funcs.R")

rhist <- stack(r1_view, c1_hotday)
names(rhist) <- c("Aged 69+ (%)", "# hot days")

r245 <- stack(r2_view, c2_hotday)
names(r245) <- c("Aged 69+ (%)", "# hot days")

r585 <- stack(r3_view, c3_hotday)
names(r585) <- c("Aged 69+ (%)", "# hot days")

r126 <- stack(r4_view, c4_hotday)
names(r126) <- c("Aged 69+ (%)", "# hot days")

r370 <- stack(r5_view, c5_hotday)
names(r370) <- c("Aged 69+ (%)", "# hot days")

for (scen in c("rhist", "r126", "r245", "r370", "r585")){
  
  scenario = ifelse(scen=="rhist", "2020", ifelse(scen=="r126", "2020-2050 (SSP126)", ifelse(scen=="r245", "2020-2050 (SSP245)", ifelse(scen=="r370", "2020-2050 (SSP370)",  "2020-2050 (SSP585)"))))
  r = get(scen)
  r = raster::crop(r, extent(regions_globe))
  
  for (i in 1:2){
    r[[i]] = clip_raster_with_sf(r[[i]], regions_globe)
  }
  
  
  region = "global"
  
  # create the bivariate raster
  bivmap <- bivariate.map(rasterx = r[["Aged.69....."]], rastery = r[["X..hot.days"]],
                          export.colour.matrix = TRUE, outname = "bivMapCols",
                          colormatrix = col.matrix, nquantiles = nBreaks)
  
  bivMapDF <- as.data.frame(bivmap, xy = TRUE) %>%
    tbl_df() %>%
    dplyr::rename("BivValue" = 3) %>%
    pivot_longer(., names_to = "Variable", values_to = "bivVal", cols = BivValue)
  
  # Make the map using ggplot
  map <- ggplot(data=bivMapDF) +
    ggtitle(paste0("Share of individuals aged 69+ and frequency of hot days, ", scenario, " , ", region))+
    geom_raster(data=bivMapDF, aes(x = x, y = y, fill = bivVal)) +
    scale_fill_gradientn(colours = col.matrix, na.value = "transparent")+
    theme_classic()+
    theme(legend.position = "none", axis.line=element_blank(),axis.text.x=element_blank(), legend.title=element_text(size=8, face = "bold"), panel.border = element_rect(colour = "transparent", fill=NA, size=0),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    guides(fill = guide_colorsteps(barwidth = 1, barheight = 5))+
    xlab("")+
    ylab("")+
    geom_sf(data=regions, fill=NA, size =.3)
  
  ggsave(plot = map,
         filename = paste0("map_biv_hotdays_", region, "_", scen, ".pdf"),
         device = "pdf", path = getwd(),
         dpi = 320, scale=1)
  
}

###################
### Figure 2

asinh_trans <- function() {
  # Source: https://stackoverflow.com/questions/14504869/histogram-with-negative-logarithmic-scale-in-r
  scales::trans_new(name = "asinh", transform = function(x) asinh(x), inverse = function(x) sinh(x))
}

# calculate CDDs threshold


df1 <- data.frame(exp = values(c1), tmax=values(c1_tasmax))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

summary(df1$exp[df1$tmax>37.5])

# > summary(df1$exp[df1$tmax>37.5])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.8592  717.0837 1214.8964 1204.9524 1635.2951 2900.5663 

###############
###############

# generate region raster

regions <- st_as_sf(rworldmap::countriesLow)
regions$REGION_n <- as.numeric(regions$REGION)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="REGION_n")
regions <- dplyr::select(regions, REGION, REGION_n)
regions$geometry<-NULL


# global figure

df1 <- data.frame(values(sum(r1_g)),  exp = values(c1), values(multigcm_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050)),  exp = values(c2), values(multigcm_245))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050)), exp = values(c3), values(multigcm_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050)), exp = values(c4), values(multigcm_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050)), exp = values(c5), values(multigcm_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

# df1 <- df1 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..)/sum(values.sum.r1_g.., na.rm=T))
# df2 <- df2 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..)/sum(values.sum.r1_g.., na.rm=T))
# df3 <- df3 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..)/sum(values.sum.r1_g.., na.rm=T))

# df1 <- df1 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..))
# df2 <- df2 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..))
# df3 <- df3 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..))


df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.sum.r1_g.. = cumsum(values.sum.r1_g..))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_a <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.sum.r1_g../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("Population, total, billion") + xlab("Yearly CDD exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + scale_x_continuous(breaks = seq(0, 5000, 1000)) + theme(legend.position = "none", legend.direction = "horizontal") + geom_vline(xintercept = 1200, linetype="dashed", colour="purple")

##

#

df1 <- data.frame(values(sum(r1_g)),  exp = values(c1), reg = values(regions_r), values(multigcm_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050)),  exp = values(c2), reg = values(regions_r), values(multigcm_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050)), exp = values(c3), reg = values(regions_r), values(multigcm_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050)), exp = values(c4), reg = values(regions_r), values(multigcm_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050)), exp = values(c5), reg = values(regions_r), values(multigcm_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_b <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("Population-weighted CDD exposure") + ylab("Region") + theme(legend.position = "none", legend.direction = "horizontal") + geom_vline(xintercept = 1200, linetype="dashed", colour="purple")

#

df1 <- data.frame(values(r1_g[[3]]),  exp = values(c1), values(multigcm_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050[[3]]),  exp = values(c2), values(multigcm_245))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(r1_ssp5_2050[[3]]), exp = values(c3), values(multigcm_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(r1_ssp1_2050[[3]]), exp = values(c4), values(multigcm_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(r1_ssp3_2050[[3]]), exp = values(c5), values(multigcm_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

#

# exposure change

sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>1200]) / sum(df$values.r1_g..3...[df$scen=="Historical"])

sum(df$values.r1_g..3...[df$scen=="SSP1(26)" & df$exp>1200]) / sum(df$values.r1_g..3...[df$scen=="SSP1(26)"])
sum(df$values.r1_g..3...[df$scen=="SSP2(45)" & df$exp>1200]) / sum(df$values.r1_g..3...[df$scen=="SSP2(45)"])
sum(df$values.r1_g..3...[df$scen=="SSP3(70)" & df$exp>1200]) / sum(df$values.r1_g..3...[df$scen=="SSP3(70)"])
sum(df$values.r1_g..3...[df$scen=="SSP5(85)" & df$exp>1200]) / sum(df$values.r1_g..3...[df$scen=="SSP5(85)"])

sum(df$values.r1_g..3...[df$scen=="SSP1(26)" & df$exp>1200])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>1200])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP2(45)" & df$exp>1200])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>1200])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP3(70)" & df$exp>1200])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>1200])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP5(85)" & df$exp>1200])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>1200])/1e9

#####

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.r1_g..3... = cumsum(values.r1_g..3...))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_c <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.r1_g..3.../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("69+ population, billion") + xlab("Yearly CDD exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + scale_x_continuous(breaks = seq(0, 5000, 1000)) + theme(legend.position = "none", legend.direction = "horizontal") + geom_vline(xintercept = 1200, linetype="dashed", colour="purple")

#, trans="asinh"

##

#

df1 <- data.frame(values(r1_g[[3]]),  exp = values(c1), reg = values(regions_r), values(multigcm_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050[[3]]),  exp = values(c2), reg = values(regions_r), values(multigcm_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(r1_ssp5_2050[[3]]), exp = values(c3), reg = values(regions_r), values(multigcm_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(r1_ssp1_2050[[3]]), exp = values(c4), reg = values(regions_r), values(multigcm_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(r1_ssp3_2050[[3]]), exp = values(c5), reg = values(regions_r), values(multigcm_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL

# exposure level change, by region
pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="SSP1(26)"], values.r1_g..3...[scen=="SSP1(26)"]))) - pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="Historical"], values.r1_g..3...[scen=="Historical"])))

pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="SSP2(45)"], values.r1_g..3...[scen=="SSP2(45)"]))) - pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="Historical"], values.r1_g..3...[scen=="Historical"])))

pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="SSP3(70)"], values.r1_g..3...[scen=="SSP3(70)"]))) - pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="Historical"], values.r1_g..3...[scen=="Historical"])))

pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="SSP5(85)"], values.r1_g..3...[scen=="SSP5(85)"]))) - pull(df %>% group_by(reg) %>% dplyr::summarise(a = weighted.mean(exp[scen=="Historical"], values.r1_g..3...[scen=="Historical"])))

###

df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_d <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("69+ population-weighted CDD exposure") + ylab("Region") + theme(legend.position = "bottom", legend.direction = "horizontal") + geom_vline(xintercept = 1200, linetype="dashed", colour="purple")

#

df1 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c1), values(multigcm_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050[[1:2]])),  exp = values(c2), values(multigcm_245))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050[[1:2]])), exp = values(c3), values(multigcm_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050[[1:2]])),  exp = values(c4), values(multigcm_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050[[1:2]])), exp = values(c5), values(multigcm_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.sum.r1_g..1.2.... = cumsum(values.sum.r1_g..1.2....))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_e <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.sum.r1_g..1.2..../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("<69 population, billion") + xlab("Yearly CDD exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + scale_x_continuous(breaks = seq(0, 5000, 1000)) + theme(legend.position = "none", legend.direction = "horizontal") + geom_vline(xintercept = 1200, linetype="dashed", colour="purple")

##

#

df1 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c1), reg = values(regions_r), values(multigcm_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c2), reg = values(regions_r), values(multigcm_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_g[[1:2]])), exp = values(c3), reg = values(regions_r), values(multigcm_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_g[[1:2]])), exp = values(c4), reg = values(regions_r), values(multigcm_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_g[[1:2]])), exp = values(c5), reg = values(regions_r), values(multigcm_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_f <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("69+ population weighted CDD exposure") + ylab("Region") + theme(legend.position = "bottom", legend.direction = "horizontal") + geom_vline(xintercept = 1200, linetype="dashed", colour="purple")

##


fig_1 <- fig_1_panel_a + fig_1_panel_b + fig_1_panel_c + fig_1_panel_d + fig_1_panel_e + fig_1_panel_f + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect", nrow = 3) & theme(legend.position = 'bottom') 

ggsave("fig_2_CDD.pdf", fig_1, height=9.5, width = 10)

#############

# global figure

df1 <- data.frame(values(sum(r1_g)),  exp = values(c1_tasmax), values(multigcm_tasmax_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050)),  exp = values(c2_tasmax), values(multigcm_tasmax_245))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050)), exp = values(c3_tasmax), values(multigcm_tasmax_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050)), exp = values(c4_tasmax), values(multigcm_tasmax_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050)), exp = values(c5_tasmax), values(multigcm_tasmax_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

# df1 <- df1 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..)/sum(values.sum.r1_g.., na.rm=T))
# df2 <- df2 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..)/sum(values.sum.r1_g.., na.rm=T))
# df3 <- df3 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..)/sum(values.sum.r1_g.., na.rm=T))

# df1 <- df1 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..))
# df2 <- df2 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..))
# df3 <- df3 %>% arrange(exp) %>% mutate(f1=cumsum(values.sum.r1_g..))


df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.sum.r1_g.. = cumsum(values.sum.r1_g..))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_g <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.sum.r1_g../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("Population, total, billion") + xlab("Yearly TMAX95 exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + theme(legend.position = "none", legend.direction = "horizontal")+ geom_vline(xintercept = 37.5, linetype="dashed", colour="purple") + coord_cartesian(xlim = c(20, 60))

df1 <- data.frame(values(sum(r1_g)),  exp = values(c1_tasmax), reg = values(regions_r), values(multigcm_tasmax_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050)),  exp = values(c2_tasmax), reg = values(regions_r), values(multigcm_tasmax_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050)), exp = values(c3_tasmax), reg = values(regions_r), values(multigcm_tasmax_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050)), exp = values(c4_tasmax), reg = values(regions_r), values(multigcm_tasmax_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050)), exp = values(c5_tasmax), reg = values(regions_r), values(multigcm_tasmax_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_h <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("Population-weighted TMAX95 exposure") + ylab("Region") + theme(legend.position = "none", legend.direction = "horizontal")+ geom_vline(xintercept = 37.5, linetype="dashed", colour="purple") + coord_cartesian(xlim = c(20, 60))

#####


df1 <- data.frame(values(r1_g[[3]]),  exp = values(c1_tasmax), values(multigcm_tasmax_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050[[3]]),  exp = values(c2_tasmax), values(multigcm_tasmax_245))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(r1_ssp5_2050[[3]]), exp = values(c3_tasmax), values(multigcm_tasmax_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(r1_ssp1_2050[[3]]), exp = values(c4_tasmax), values(multigcm_tasmax_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(r1_ssp3_2050[[3]]), exp = values(c5_tasmax), values(multigcm_tasmax_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

# exposure change

sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>37.5]) / sum(df$values.r1_g..3...[df$scen=="Historical"])

sum(df$values.r1_g..3...[df$scen=="SSP1(26)" & df$exp>37.5]) / sum(df$values.r1_g..3...[df$scen=="SSP1(26)"])
sum(df$values.r1_g..3...[df$scen=="SSP2(45)" & df$exp>37.5]) / sum(df$values.r1_g..3...[df$scen=="SSP2(45)"])
sum(df$values.r1_g..3...[df$scen=="SSP3(70)" & df$exp>37.5]) / sum(df$values.r1_g..3...[df$scen=="SSP3(70)"])
sum(df$values.r1_g..3...[df$scen=="SSP5(85)" & df$exp>37.5]) / sum(df$values.r1_g..3...[df$scen=="SSP5(85)"])

sum(df$values.r1_g..3...[df$scen=="SSP1(26)" & df$exp>37.5])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>37.5])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP2(45)" & df$exp>37.5])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>37.5])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP3(70)" & df$exp>37.5])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>37.5])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP5(85)" & df$exp>37.5])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>37.5])/1e9

#

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.r1_g..3... = cumsum(values.r1_g..3...))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_i <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.r1_g..3.../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("69+ population, billion") + xlab("Yearly TMAX95 exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + theme(legend.position = "none", legend.direction = "horizontal") + geom_vline(xintercept = 37.5, linetype="dashed", colour="purple") + coord_cartesian(xlim = c(20, 60))

##

#

df1 <- data.frame(values(r1_g[[3]]),  exp = values(c1_tasmax), reg = values(regions_r), values(multigcm_tasmax_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050[[3]]),  exp = values(c2_tasmax), reg = values(regions_r), values(multigcm_tasmax_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(r1_ssp5_2050[[3]]), exp = values(c3_tasmax), reg = values(regions_r), values(multigcm_tasmax_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(r1_ssp1_2050[[3]]), exp = values(c4_tasmax), reg = values(regions_r), values(multigcm_tasmax_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(r1_ssp3_2050[[3]]), exp = values(c5_tasmax), reg = values(regions_r), values(multigcm_tasmax_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_j <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("69+ population-weighted TMAX95 exposure") + ylab("Region") + theme(legend.position = "bottom", legend.direction = "horizontal")+ coord_cartesian(xlim = c(20, 60)) + geom_vline(xintercept = 37.5, linetype="dashed", colour="purple")

#

df1 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c1_tasmax), values(multigcm_tasmax_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050[[1:2]])),  exp = values(c2_tasmax), values(multigcm_tasmax_245))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050[[1:2]])), exp = values(c3_tasmax), values(multigcm_tasmax_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050[[1:2]])), exp = values(c4_tasmax), values(multigcm_tasmax_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050[[1:2]])), exp = values(c5_tasmax), values(multigcm_tasmax_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.sum.r1_g..1.2.... = cumsum(values.sum.r1_g..1.2....))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_k <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.sum.r1_g..1.2..../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("<69 population, billion") + xlab("Yearly TMAX95 exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + theme(legend.position = "none", legend.direction = "horizontal")+ geom_vline(xintercept = 37.5, linetype="dashed", colour="purple") + coord_cartesian(xlim = c(20, 60))


df1 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c1_tasmax), reg = values(regions_r), values(multigcm_tasmax_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c2_tasmax), reg = values(regions_r), values(multigcm_tasmax_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_g[[1:2]])), exp = values(c3_tasmax), reg = values(regions_r), values(multigcm_tasmax_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]


df4 <- data.frame(values(sum(r1_g[[1:2]])), exp = values(c4_tasmax), reg = values(regions_r), values(multigcm_tasmax_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_g[[1:2]])), exp = values(c5_tasmax), reg = values(regions_r), values(multigcm_tasmax_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_l <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("<69 population-weighted TMAX95 exposure") + ylab("Region") + theme(legend.position = "bottom", legend.direction = "horizontal")+ geom_vline(xintercept = 37.5, linetype="dashed", colour="purple")+ coord_cartesian(xlim = c(20, 60))


fig_1 <- fig_1_panel_g + fig_1_panel_h + fig_1_panel_i + fig_1_panel_j  + fig_1_panel_k + fig_1_panel_l + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect", nrow = 3) & theme(legend.position = 'bottom') 

ggsave("fig_2_TMAX.pdf", fig_1, height=9.5, width = 10)


###################
### Figure 2

# generate region raster

regions <- st_as_sf(rworldmap::countriesLow)
regions$REGION_n <- as.numeric(regions$REGION)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="REGION_n")
regions <- dplyr::select(regions, REGION, REGION_n)
regions$geometry<-NULL


# global figure

df1 <- data.frame(values(sum(r1_g)),  exp = values(c1_hotday), values(multigcm_hotday_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050)),  exp = values(c2_hotday), values(multigcm_hotday_245[[1:13]]))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050)), exp = values(c3_hotday), values(multigcm_hotday_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]


df4 <- data.frame(values(sum(r1_ssp1_2050)), exp = values(c4_hotday), values(multigcm_hotday_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050)), exp = values(c5_hotday), values(multigcm_hotday_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)


df$id <- NULL

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.sum.r1_g.. = cumsum(values.sum.r1_g..))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_m <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.sum.r1_g../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("Population, total, billion") + xlab("Yearly #HD / yr. exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + scale_x_continuous(breaks = seq(0, 200, 50)) + theme(legend.position = "none", legend.direction = "horizontal")

##

#

df1 <- data.frame(values(sum(r1_g)),  exp = values(c1_hotday), reg = values(regions_r), values(multigcm_hotday_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050)),  exp = values(c2_hotday), reg = values(regions_r), values(multigcm_hotday_245[[1:13]]), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050)), exp = values(c3_hotday), reg = values(regions_r), values(multigcm_hotday_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050)), exp = values(c4_hotday), reg = values(regions_r), values(multigcm_hotday_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050)), exp = values(c5_hotday), reg = values(regions_r), values(multigcm_hotday_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_n <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("Population-weighted #HD / yr. exposure") + ylab("Region") + theme(legend.position = "none", legend.direction = "horizontal")+ coord_cartesian(xlim = c(0, 200))

#

df1 <- data.frame(values(r1_g[[3]]),  exp = values(c1_hotday), values(multigcm_hotday_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050[[3]]),  exp = values(c2_hotday), values(multigcm_hotday_245[[1:13]]))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(r1_ssp5_2050[[3]]), exp = values(c3_hotday), values(multigcm_hotday_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(r1_ssp1_2050[[3]]), exp = values(c4_hotday), values(multigcm_hotday_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(r1_ssp3_2050[[3]]), exp = values(c5_hotday), values(multigcm_hotday_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

# exposure change

sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>30]) / sum(df$values.r1_g..3...[df$scen=="Historical"])

sum(df$values.r1_g..3...[df$scen=="SSP1(26)" & df$exp>30]) / sum(df$values.r1_g..3...[df$scen=="SSP1(26)"])
sum(df$values.r1_g..3...[df$scen=="SSP2(45)" & df$exp>30]) / sum(df$values.r1_g..3...[df$scen=="SSP2(45)"])
sum(df$values.r1_g..3...[df$scen=="SSP3(70)" & df$exp>30]) / sum(df$values.r1_g..3...[df$scen=="SSP3(70)"])
sum(df$values.r1_g..3...[df$scen=="SSP5(85)" & df$exp>30]) / sum(df$values.r1_g..3...[df$scen=="SSP5(85)"])

sum(df$values.r1_g..3...[df$scen=="SSP1(26)" & df$exp>30])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>30])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP2(45)" & df$exp>30])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>30])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP3(70)" & df$exp>30])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>30])/1e9
sum(df$values.r1_g..3...[df$scen=="SSP5(85)" & df$exp>30])/1e9 - sum(df$values.r1_g..3...[df$scen=="Historical" & df$exp>30])/1e9

#

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.r1_g..3... = cumsum(values.r1_g..3...))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_o <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.r1_g..3.../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("69+ population, billion") + xlab("Yearly #HD / yr. exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + scale_x_continuous(breaks = seq(0, 200, 50)) + theme(legend.position = "none", legend.direction = "horizontal")

##

#

df1 <- data.frame(values(r1_g[[3]]),  exp = values(c1_hotday), reg = values(regions_r), values(multigcm_hotday_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050[[3]]),  exp = values(c2_hotday), reg = values(regions_r), values(multigcm_hotday_245[[1:13]]), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(r1_ssp5_2050[[3]]), exp = values(c3_hotday), reg = values(regions_r), values(multigcm_hotday_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(r1_ssp1_2050[[3]]), exp = values(c4_hotday), reg = values(regions_r), values(multigcm_hotday_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(r1_ssp3_2050[[3]]), exp = values(c5_hotday), reg = values(regions_r), values(multigcm_hotday_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)


df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_p <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("69+ population-weighted #HD / yr. exposure") + ylab("Region") + theme(legend.position = "bottom", legend.direction = "horizontal")+ coord_cartesian(xlim = c(0, 200))

#

df1 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c1_hotday), values(multigcm_hotday_hist))
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050[[1:2]])),  exp = values(c2_hotday), values(multigcm_hotday_245[[1:13]]))
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050[[1:2]])), exp = values(c3_hotday), values(multigcm_hotday_585))
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050[[1:2]])), exp = values(c4_hotday), values(multigcm_hotday_126))
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050[[1:2]])), exp = values(c5_hotday), values(multigcm_hotday_370))
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$id <- NULL

df <- reshape2::melt(df, c(1,17))

df <- df %>%
  group_by(scen, variable) %>% 
  arrange(value) %>%
  dplyr::mutate(values.sum.r1_g..1.2.... = cumsum(values.sum.r1_g..1.2....))

df$type <- ifelse(df$variable=="exp", "Median", "GCMs range")

fig_1_panel_q <- ggplot(df %>% sample_frac(0.001), aes(x=value, y=values.sum.r1_g..1.2..../1e9, colour=scen, group=interaction(scen,variable), alpha=type)) + theme_classic() + geom_step(show.legend = FALSE) + ylab("<69 population, billion") + xlab("Yearly #HD / yr. exposure (CMIP6)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + scale_x_continuous(breaks = seq(0, 200, 50)) + theme(legend.position = "none", legend.direction = "horizontal")

##

#

df1 <- data.frame(values(sum(r1_g[[1:2]])),  exp = values(c1_hotday), reg = values(regions_r), values(multigcm_hotday_hist), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(sum(r1_ssp2_2050[[1:2]])),  exp = values(c2_hotday), reg = values(regions_r), values(multigcm_hotday_245[[1:13]]), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1] <- colnames(df1)[1]

df3 <- data.frame(values(sum(r1_ssp5_2050[[1:2]])), exp = values(c3_hotday), reg = values(regions_r), values(multigcm_hotday_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1] <- colnames(df1)[1]

df4 <- data.frame(values(sum(r1_ssp1_2050[[1:2]])), exp = values(c4_hotday), reg = values(regions_r), values(multigcm_hotday_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1] <- colnames(df1)[1]

df5 <- data.frame(values(sum(r1_ssp3_2050[[1:2]])), exp = values(c5_hotday), reg = values(regions_r), values(multigcm_hotday_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1] <- colnames(df1)[1]

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL
df2 <- reshape2::melt(df, c(1, 3, 18:19))
colnames(df2)[1] <- "pop"

df2_m <- df2 %>% filter(variable=="exp") %>%  group_by(scen, reg) %>% dplyr::summarise(med_gcm = weighted.median(value, pop, na.rm=T))

fig_1_panel_r <- ggplot()+ theme_classic() + geom_boxplot(data=df2 %>% slice_sample(prop=.01), aes(y=reg, x=value, fill=scen, weight=pop), outlier.shape=NA) + scale_fill_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") + geom_point(data=df2_m, aes(y=reg, x=med_gcm, colour=scen), show.legend = F, position=position_dodge(0.75), shape=23, size=2.5) + scale_color_manual(values=rep("cyan", 5)) + xlab("<69 population-weighted #HD / yr. exposure") + ylab("Region") + theme(legend.position = "bottom", legend.direction = "horizontal")+ coord_cartesian(xlim = c(0, 200))

fig_1 <- fig_1_panel_m + fig_1_panel_n + fig_1_panel_o + fig_1_panel_p + fig_1_panel_q + fig_1_panel_r + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect", nrow = 3) & theme(legend.position = 'bottom') 

ggsave("fig_2_hotdays.pdf", fig_1, height=9.5, width = 10)

####



fig_1 <- fig_1_panel_c + fig_1_panel_d + fig_1_panel_o + fig_1_panel_p + fig_1_panel_i + fig_1_panel_j  + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect", nrow = 3) & theme(legend.position = 'bottom') 

ggsave("fig_2_new.pdf", fig_1, height=9.5, width = 10)
ggsave("fig_2_new.png", fig_1, height=9.5, width = 10)

# https://www.thelancet.com/journals/lanplh/article/PIIS2542-5196(21)00079-6/fulltext#:~:text=Short%20exposures%20to%20temperatures%20above,%2C%20fish%2C%20and%20agricultural%20crops.[ Short exposures to temperatures above 35?C with high humidity, or above 40?C with low humidity, can be lethal]


###################
### Figure 3 

df <- data.frame(values(r1_g), exp = values(c1), reg = values(regions_r), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% arrange(exp) %>% mutate(f1=cumsum(layer.1+layer.2)/sum(layer.1+layer.2), f3=cumsum(layer)/sum(layer))

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("r1"))
df2 <- reshape2::melt(df, 1:3)

df2$variable <- ifelse(df2$variable=="f1", "Age group <69", as.character(df2$variable))
df2$variable <- ifelse(df2$variable=="f3", "Age group 69+", as.character(df2$variable))

g <- ggplot(df2, aes(x=exp, y=value, group=reg, colour=reg)) + theme_classic() + geom_line() + ylab("Share of sub-population (%)") + xlab("Yearly CDD exposure (CMIP6 1995 - 2014)")+ggsci::scale_colour_npg(name="Region") + facet_wrap(vars(variable)) + scale_x_continuous(breaks = seq(0, 2500, 500)) + scale_y_continuous(breaks = seq(0, 1, .1), labels = scales::label_percent())

ggsave("fig_3.pdf", g, width = 10)

# 
# h <- ggplot(df2, aes(x=exp, y=value, group=variable, colour=variable)) + geom_line() + ylab("Share of sub-population (%)") + xlab("Yearly CDD exposure (CMIP6 1995 - 2014)")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Region") + facet_wrap(vars(reg))
# 
# 
# h


###################
### Figure 4

names(r1_g)[3] <- "layer.3"

df1 <- data.frame(values(r1_g),  exp = values(c1), reg = values(regions_r), xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050),  exp = values(c2), reg = values(regions_r), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1:3] <- colnames(df1)[1:3]

df3 <- data.frame(values(r1_ssp5_2050), exp = values(c3), reg = values(regions_r), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1:3] <- colnames(df1)[1:3]

df4 <- data.frame(values(r1_ssp1_2050), exp = values(c4), reg = values(regions_r), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1:3] <- colnames(df1)[1:3]

df5 <- data.frame(values(r1_ssp3_2050), exp = values(c5), reg = values(regions_r), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1:3] <- colnames(df1)[1:3]

df1 <- df1 %>% dplyr::group_by(reg) %>% arrange(exp) %>% mutate(f1=cumsum(layer.1+layer.2)/sum(layer.1+layer.2), f3=cumsum(layer.3)/sum(layer.3))

df2 <- df2 %>% dplyr::group_by(reg) %>% arrange(exp) %>% mutate(f1=cumsum(layer.1+layer.2)/sum(layer.1+layer.2), f3=cumsum(layer.3)/sum(layer.3))

df3 <- df3 %>% dplyr::group_by(reg) %>% arrange(exp) %>% mutate(f1=cumsum(layer.1+layer.2)/sum(layer.1+layer.2), f3=cumsum(layer.3)/sum(layer.3))

df4 <- df4 %>% dplyr::group_by(reg) %>% arrange(exp) %>% mutate(f1=cumsum(layer.1+layer.2)/sum(layer.1+layer.2), f3=cumsum(layer.3)/sum(layer.3))

df5 <- df5 %>% dplyr::group_by(reg) %>% arrange(exp) %>% mutate(f1=cumsum(layer.1+layer.2)/sum(layer.1+layer.2), f3=cumsum(layer.3)/sum(layer.3))

df1$scen <- "Historical"
df2$scen <- "SSP2(45)"
df3$scen <- "SSP5(85)"
df4$scen <- "SSP1(26)"
df5$scen <- "SSP3(70)"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL

df <- dplyr::select(df, -starts_with("layer"))
df2 <- reshape2::melt(df, c(1:3, 6))

df2$variable <- ifelse(df2$variable=="f1", "Age group <69", as.character(df2$variable))
df2$variable <- ifelse(df2$variable=="f3", "Age group 69+", as.character(df2$variable))

g2 <- ggplot(df2, aes(x=exp+1, y=value, group=scen, colour=scen)) +
  theme_classic() + geom_line() + ylab("Share of sub-population (%)") +
  xlab("Yearly CDD exposure")+scale_color_manual(values=c("#302f2b", "#fcec8d", "#ffbe0a", "#ff9830", "#b51209"),name="Scenario") +
  facet_wrap(vars(reg, variable), ncol = 6, scales = "free") +
  scale_y_continuous(breaks = seq(0, 1, .25), labels = scales::label_percent())  +
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave("g2.pdf", g2, width = 11, height = 6, scale=0.85)

###################
### Figure 5 


regions <- st_as_sf(rworldmap::countriesLow)
regions$REGION_n <- as.numeric(regions$REGION)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="REGION_n")
regions <- dplyr::select(regions, REGION, REGION_n)
regions$geometry<-NULL

df1 <- data.frame(values(r1_g),  exp = values(c1), reg = values(regions_r), values(multigcm_hist),  xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050),  exp = values(c2), reg = values(regions_r), values(multigcm_245[[1:13]]), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1:3] <- colnames(df1)[1:3]

df3 <- data.frame(values(r1_ssp5_2050), exp = values(c3), reg = values(regions_r), values(multigcm_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1:3] <- colnames(df1)[1:3]

df4 <- data.frame(values(r1_ssp1_2050), exp = values(c4), reg = values(regions_r), values(multigcm_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1:3] <- colnames(df1)[1:3]

df5 <- data.frame(values(r1_ssp3_2050), exp = values(c5), reg = values(regions_r), values(multigcm_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1:3] <- colnames(df1)[1:3]

df1$scen <- "Baseline (historical climate, 2020 population)"
df2$scen <- "SSP2(45), year 2050"
df3$scen <- "SSP5(85), year 2050"
df4$scen <- "SSP1(26), year 2050"
df5$scen <- "SSP3(70), year 2050"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL

####

df2 <- reshape2::melt(df, c(4:21))

df2$variable <- ifelse(df2$variable=="layer.1" | df2$variable=="layer.2", "Age group <69", as.character(df2$variable))
df2$variable <- ifelse(df2$variable=="layer", "Age group 69+", as.character(df2$variable))

colnames(df2)[19:20]  <- c("pop_group", "pop_value")

#gdata::keep(df2, sure=T)
gc()

# extract numbers

stat <- df2 %>% group_by(reg, scen, pop_group) %>% dplyr::summarise(exp_r = weighted.mean(exp, pop_value, na.rm=T))

stat <- stat %>% group_by(reg) %>% dplyr::summarise(min_exp_r=min(exp_r), max_exp_r=max(exp_r))

write.csv(stat, "stat_cdd_bins.csv")

## run a statistical test

stat <- (df2 %>% group_by(reg, scen) %>% dplyr::summarise(t= as.numeric(wtd.t.test(exp[pop_group=="Age group <69"], exp[pop_group=="Age group 69+"], weight=pop_value[pop_group=="Age group <69"], weighty=pop_value[pop_group=="Age group 69+"])$coefficients[1]), p= as.numeric(wtd.t.test(exp[pop_group=="Age group <69"], exp[pop_group=="Age group 69+"], weight=pop_value[pop_group=="Age group <69"], weighty=pop_value[pop_group=="Age group 69+"])$coefficients[3])))

stat$sig <- ifelse(stat$p < 0.1 & p>=0.05, "*", ifelse(stat$p < 0.05 & p>=0.01, "**", "***"))

write.csv(stat, "stat_cdd_bins_test.csv")

r <- read.csv("stat_cdd_bins_test.csv")
r$X <- NULL

stargazer::stargazer(r, summary=F, out = "ttest_tab.tex", rownames = FALSE)

####

rm(df, df1, df4, df5)
gc()

df3 <- reshape2::melt(df2, c(2, 17:20))

rm(df2)
gc()

df3$exp_q <- NA
df3$exp_q[df3$value < 30] <- "0-30"
df3$exp_q[df3$value >= 30 & df3$value < 350] <-  "30-350" 
df3$exp_q[df3$value >= 350 & df3$value < 1000] <- "350-1000" 
df3$exp_q[df3$value >= 1000 & df3$value < 2500] <- "1000-2500"
df3$exp_q[df3$value > 2500] <- ">2500"

df3$value = NULL

df2_s <- dplyr::group_by(df3, scen, exp_q, reg, variable, pop_group) %>% dplyr::summarise(value=sum(pop_value, na.rm=T))

rm(df3)
gc()

df2_s$exp_q <- factor(df2_s$exp_q, levels = c("0-30", "30-350", "350-1000", "1000-2500", ">2500"))

###

df2_s2 <- df2_s %>% filter(variable=="exp")

df2_s2$value <- df2_s2$value/1e6

df2_s2 <- df2_s2 %>% 
  mutate_if(is.numeric, round, 0)

datasummary(value * reg * scen ~ pop_group * exp_q * (Mean),  fmt=0, data = df2_s2, output = "table_si2.tex")


###

df2_s <-df2_s %>% ungroup() %>% complete(exp_q, scen, reg, variable, pop_group, fill = list(value=0))

df2_s <- dplyr::group_by(df2_s, exp_q, reg, variable, pop_group) %>% dplyr::mutate(value=value - value[scen == "Baseline (historical climate, 2020 population)"])

write.csv(df2_s, "table_si2.csv")

df2_s <- filter(df2_s, scen != "Baseline (historical climate, 2020 population)")

#

df2_s$exp_q = factor(df2_s$exp_q, ordered = T)
df2_s$exp_q = as.numeric(df2_s$exp_q)

g3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP5(85), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP5(85), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("Population change from 2020 (million people), range across GCMs") + xlab("CDD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") + scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_world <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP5(85), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP5(85), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#

g3_si <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP2(45), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP2(45), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("Population change from 2020 (million people), range across GCMs") + xlab("CDD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") + scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_si_world <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP2(45), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP2(45), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)
  
##

g3_si2 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP1(26), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP1(26), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("Population change from 2020 (million people), range across GCMs") + xlab("CDD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") + scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_si_world2 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP1(26), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP1(26), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#

g3_si3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP3(70), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP3(70), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("Population change from 2020 (million people), range across GCMs") + xlab("CDD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") + scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_si_world3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP3(70), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP3(70), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("0", "30+", "350+", "1000+", "2500+")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#############################################

df1 <- data.frame(values(r1_g),  exp = values(c1_tasmax), reg = values(regions_r), values(multigcm_tasmax_hist),  xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050),  exp = values(c2_tasmax), reg = values(regions_r), values(multigcm_tasmax_245), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1:3] <- colnames(df1)[1:3]

df3 <- data.frame(values(r1_ssp5_2050), exp = values(c3_tasmax), reg = values(regions_r), values(multigcm_tasmax_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1:3] <- colnames(df1)[1:3]

df4 <- data.frame(values(r1_ssp1_2050), exp = values(c4_tasmax), reg = values(regions_r), values(multigcm_tasmax_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1:3] <- colnames(df1)[1:3]

df5 <- data.frame(values(r1_ssp3_2050), exp = values(c5_tasmax), reg = values(regions_r), values(multigcm_tasmax_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1:3] <- colnames(df1)[1:3]

df1$scen <- "Baseline (historical climate, 2020 population)"
df2$scen <- "SSP2(45), year 2050"
df3$scen <- "SSP5(85), year 2050"
df4$scen <- "SSP1(26), year 2050"
df5$scen <- "SSP3(70), year 2050"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL

####

df2 <- reshape2::melt(df, c(4:21))

df2$variable <- ifelse(df2$variable=="layer.1" | df2$variable=="layer.2", "Age group <69", as.character(df2$variable))
df2$variable <- ifelse(df2$variable=="layer", "Age group 69+", as.character(df2$variable))

colnames(df2)[19:20]  <- c("pop_group", "pop_value")

df3 <- reshape2::melt(df2, c(2, 17:20))

rm(df, df1, df2, df4, df5)
gc()

df3$value_q <- NA
df3$value_q[df3$value < 24] <- "<24"
df3$value_q[df3$value >= 24 & df3$value < 30] <-  "24-30" 
df3$value_q[df3$value >= 30 & df3$value < 35] <- "30-35"
df3$value_q[df3$value >= 35 & df3$value < 37.5] <- "35-37.5"
df3$value_q[df3$value > 37.5] <- ">37.5"

df3$value = NULL

df2_s <- dplyr::group_by(df3, scen, value_q, reg, variable, pop_group) %>% dplyr::summarise(value=sum(pop_value, na.rm=T))

rm(df3)
gc()

df2_s$exp_q <- factor(df2_s$value_q, levels = c("<24", "24-30", "30-35", "35-37.5", ">37.5"))

#

df2_s <-df2_s %>% ungroup() %>% complete(exp_q, scen, reg, variable, pop_group, fill = list(value=0))

df2_s <- dplyr::group_by(df2_s, exp_q, reg, variable, pop_group) %>% dplyr::mutate(value=value - value[scen == "Baseline (historical climate, 2020 population)"])

df2_s <- filter(df2_s, scen != "Baseline (historical climate, 2020 population)")

##

df2_s$exp_q = factor(df2_s$exp_q, ordered = T)
df2_s$exp_q = as.numeric(df2_s$exp_q)

g3_b <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP5(85), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP5(85), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("Population change from 2020 (million people), range across GCMs") + xlab("TMAX95 exposure (°C)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group")+ scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_b_world <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP5(85), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP5(85), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#

g3_b_si <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) +  geom_smooth(data=df2_s %>% filter(scen=="SSP2(45), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP2(45), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("TMAX95 exposure (°C)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") +  scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_b_si_world <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP2(45), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP2(45), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

##

g3_b_si_2 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) +  geom_smooth(data=df2_s %>% filter(scen=="SSP1(26), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP1(26), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("TMAX95 exposure (°C)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") + scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_b_si_world_2 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP1(26), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP1(26), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#

g3_b_si_3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) +  geom_smooth(data=df2_s %>% filter(scen=="SSP3(70), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP3(70), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("TMAX95 exposure (°C)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group") + scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_b_si_world_3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP3(70), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP3(70), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<24", "24-30", "30-35", "35-37.5", ">37.5")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

# g3_b <- ggplot(df2_s) + theme_classic()+ geom_hline(yintercept = 0, linetype="dashed", alpha=0.75, lwd=0.15) + geom_boxplot(aes(x=exp_q, y=value/1e6, fill=pop_group), outlier.shape=NA, lwd=0.15) + ylab("") + xlab("95th percentile maximum temperature (C degrees) exposure") + facet_wrap(vars(reg, scen), scales = "free", nrow=6) + scale_fill_manual(values=c("#ffbe0a", "#b51209"),name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal")

#

df1 <- data.frame(values(r1_g),  exp = values(c1_hotday), reg = values(regions_r), values(multigcm_hotday_hist),  xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050),  exp = values(c2_hotday), reg = values(regions_r), values(multigcm_hotday_245[[1:13]]), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1:3] <- colnames(df1)[1:3]

df3 <- data.frame(values(r1_ssp5_2050), exp = values(c3_hotday), reg = values(regions_r), values(multigcm_hotday_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1:3] <- colnames(df1)[1:3]

df4 <- data.frame(values(r1_ssp1_2050), exp = values(c4_hotday), reg = values(regions_r), values(multigcm_hotday_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1:3] <- colnames(df1)[1:3]

df5 <- data.frame(values(r1_ssp3_2050), exp = values(c5_hotday), reg = values(regions_r), values(multigcm_hotday_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1:3] <- colnames(df1)[1:3]

df1$scen <- "Baseline (historical climate, 2020 population)"
df2$scen <- "SSP2(45), year 2050"
df3$scen <- "SSP5(85), year 2050"
df4$scen <- "SSP1(26), year 2050"
df5$scen <- "SSP3(70), year 2050"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL

####

df2 <- reshape2::melt(df, c(4:21))

df2$variable <- ifelse(df2$variable=="layer.1" | df2$variable=="layer.2", "Age group <69", as.character(df2$variable))
df2$variable <- ifelse(df2$variable=="layer", "Age group 69+", as.character(df2$variable))

colnames(df2)[19:20]  <- c("pop_group", "pop_value")

df3 <- reshape2::melt(df2, c(2, 17:20))

rm(df, df1, df2, df4, df5)
gc()

df3$value_q <- NA
df3$value_q[df3$value < 5] <- "<5"
df3$value_q[df3$value >= 5 & df3$value < 10] <-  "5-10" 
df3$value_q[df3$value >= 10 & df3$value < 25] <- "10-25"
df3$value_q[df3$value >= 25 & df3$value < 50] <- "25-50"
df3$value_q[df3$value > 50] <- ">50"

df3$value = NULL

df2_s <- dplyr::group_by(df3, scen, value_q, reg, variable, pop_group) %>% dplyr::summarise(value=sum(pop_value, na.rm=T))

rm(df3)
gc()

df2_s$exp_q <- factor(df2_s$value_q, levels = c("<5", "5-10", "10-25", "25-50", ">50"))

#

df2_s <-df2_s %>% ungroup() %>% complete(exp_q, scen, reg, variable, pop_group, fill = list(value=0))

df2_s <- dplyr::group_by(df2_s, exp_q, reg, variable, pop_group) %>% dplyr::mutate(value=value - value[scen == "Baseline (historical climate, 2020 population)"])

df2_s <- filter(df2_s, scen != "Baseline (historical climate, 2020 population)")

##

df2_s$exp_q = factor(df2_s$exp_q, ordered = T)
df2_s$exp_q = as.numeric(df2_s$exp_q)

g3_c <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) +  geom_smooth(data=df2_s %>% filter(scen=="SSP5(85), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP5(85), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("#HD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group")+ scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_c_world <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP5(85), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP5(85), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#

g3_c_si <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) +  geom_smooth(data=df2_s %>% filter(scen=="SSP2(45), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP2(45), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("#HD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group")+ scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_c_si_world <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP2(45), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP2(45), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

##

g3_c_si_2 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) +  geom_smooth(data=df2_s %>% filter(scen=="SSP1(26), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP1(26), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("#HD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group")+ scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_c_si_world_2 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP1(26), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP1(26), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

#

g3_c_si_3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP3(70), year 2050"), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP3(70), year 2050" & variable=="exp"), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("#HD exposure (#)") + facet_wrap(vars(reg), scales = "free", nrow=6) + scale_colour_discrete(name="Pop. group")+ scale_fill_discrete(name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

g3_c_si_world_3 <- ggplot() + theme_classic(base_size = 13) + geom_hline(yintercept = 0, alpha=0.75, lwd=0.15) + geom_smooth(data=df2_s %>% filter(scen=="SSP3(70), year 2050") %>% dplyr::group_by(variable, exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, fill=pop_group), colour="transparent", linewidth=0.15, level = 0.99) + geom_point(data=df2_s %>% filter(scen=="SSP3(70), year 2050" & variable=="exp") %>% dplyr::group_by(exp_q, pop_group) %>% dplyr::summarise(value=sum(value, na.rm=T)), aes(x=exp_q, y=value/1e6, colour=pop_group), size=1.5) + ylab("") + xlab("") + scale_fill_discrete(name="Pop. group") + scale_colour_discrete(name="Pop. group")  + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(breaks=c(1:5), labels=c("<5", "5-10", "10-25", "25-50", ">50")) + theme(strip.background = element_blank(), strip.text = element_blank(), aspect.ratio=.5)

# g3_b <- ggplot(df2_s) + theme_classic()+ geom_hline(yintercept = 0, linetype="dashed", alpha=0.75, lwd=0.15) + geom_boxplot(aes(x=exp_q, y=value/1e6, fill=pop_group), outlier.shape=NA, lwd=0.15) + ylab("") + xlab("95th percentile maximum temperature (C degrees) exposure") + facet_wrap(vars(reg, scen), scales = "free", nrow=6) + scale_fill_manual(values=c("#ffbe0a", "#b51209"),name="Pop. group") + theme(legend.position = "bottom", legend.direction = "horizontal")



g3_tot <- (g3_world + g3_c_world + g3_b_world + g3 + g3_c + g3_b)  + plot_layout(ncol=3, guides = "collect") & theme(plot.margin = margin(.1, .01, .05, .01, "cm"), legend.position = "bottom", legend.direction = "horizontal")

ggsave("g3.pdf", g3_tot, width = 8, height = 8, scale=1.35)
ggsave("g3.png", g3_tot, width = 8, height = 8, scale=1.35)

g3_tot_si <- (g3_si_world + g3_c_si_world + g3_b_si_world + g3_si + g3_c_si + g3_b_si)  + plot_layout(ncol=3, guides = "collect") & theme(plot.margin = margin(.1, .01, .05, .01, "cm"), legend.position = "bottom", legend.direction = "horizontal")

ggsave("g3_si_245.pdf", g3_tot_si, width = 8, height = 8, scale=1.35)
ggsave("g3_si_245.png", g3_tot, width = 8, height = 8, scale=1.35)

g3_tot_si <- (g3_si_world2 + g3_c_si_world_2 + g3_b_si_world_2 + g3_si2 + g3_c_si_2 + g3_b_si_2)  + plot_layout(ncol=3, guides = "collect") & theme(plot.margin = margin(.1, .01, .05, .01, "cm"), legend.position = "bottom", legend.direction = "horizontal")

ggsave("g3_si_126.pdf", g3_tot_si, width = 8, height = 8, scale=1.35)
ggsave("g3_si_126.png", g3_tot_si, width = 8, height = 8, scale=1.35)

g3_tot_si <- (g3_si_world3 + g3_c_si_world_3 + g3_b_si_world_3 + g3_si3 + g3_c_si_3 + g3_b_si_3)  + plot_layout(ncol=3, guides = "collect") & theme(plot.margin = margin(.1, .01, .05, .01, "cm"), legend.position = "bottom", legend.direction = "horizontal")

ggsave("g3_si_370.pdf", g3_tot_si, width = 8, height = 8, scale=1.35)
ggsave("g3_si_370.png", g3_tot_si, width = 8, height = 8, scale=1.35)

###

df1 <- data.frame(values(r1_g),  exp = values(c1_tasmax), reg = values(regions_r), values(multigcm_tasmax_hist),  xy=T)
df1$id <- 1:nrow(df1)
df1 <- na.omit(df1)

df2 <- data.frame(values(r1_ssp2_2050),  exp = values(c2_tasmax), reg = values(regions_r), values(multigcm_tasmax_245[[1:13]]), xy=T)
df2$id <- 1:nrow(df2)
df2 <- na.omit(df2)
colnames(df2)[1:3] <- colnames(df1)[1:3]

df3 <- data.frame(values(r1_ssp5_2050), exp = values(c3_tasmax), reg = values(regions_r), values(multigcm_tasmax_585), xy=T)
df3$id <- 1:nrow(df3)
df3 <- na.omit(df3)
colnames(df3)[1:3] <- colnames(df1)[1:3]

df4 <- data.frame(values(r1_ssp1_2050), exp = values(c4_tasmax), reg = values(regions_r), values(multigcm_tasmax_126), xy=T)
df4$id <- 1:nrow(df4)
df4 <- na.omit(df4)
colnames(df4)[1:3] <- colnames(df1)[1:3]

df5 <- data.frame(values(r1_ssp3_2050), exp = values(c5_tasmax), reg = values(regions_r), values(multigcm_tasmax_370), xy=T)
df5$id <- 1:nrow(df5)
df5 <- na.omit(df5)
colnames(df5)[1:3] <- colnames(df1)[1:3]

df1$scen <- "Baseline (historical climate, 2020 population)"
df2$scen <- "SSP2(45), year 2050"
df3$scen <- "SSP5(85), year 2050"
df4$scen <- "SSP1(26), year 2050"
df5$scen <- "SSP3(70), year 2050"

df <- bind_rows(df1, df2, df3, df4, df5)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df$id <- NULL

####

df2 <- reshape2::melt(df, c(4:21))

df2$variable <- ifelse(df2$variable=="layer.1" | df2$variable=="layer.2", "Age group <69", as.character(df2$variable))
df2$variable <- ifelse(df2$variable=="layer", "Age group 69+", as.character(df2$variable))

colnames(df2)[19:20]  <- c("pop_group", "pop_value")

df3 <- reshape2::melt(df2, c(2, 17:20))

rm(df, df1, df2, df4, df5)
gc()

df3$value_q <- df3$value

df3$value = NULL

  df2_s <- df3

rm(df3)
gc()

# calculate share of people living under TAXMAX > 37.5 in 2050

a<- df2_s %>%  ungroup() %>% filter(scen=="Baseline (historical climate, 2020 population)" & (value_q  >= 37.5) & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T)) / df2_s %>%  ungroup() %>% filter(scen=="Baseline (historical climate, 2020 population)" & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T))

b <- df2_s %>%  ungroup() %>% filter(scen=="SSP2(45), year 2050" & (value_q  >= 37.5) & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T)) / df2_s %>%  ungroup() %>% filter(scen=="SSP2(45), year 2050" & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T))

c<- df2_s %>%  ungroup() %>% filter(scen=="SSP5(85), year 2050" & (value_q  >= 37.5) & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T)) / df2_s %>%  ungroup() %>% filter(scen=="SSP5(85), year 2050" & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T))

d <- (df2_s %>%  ungroup() %>% filter(scen=="SSP5(85), year 2050" & (value_q  >= 37.5) & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T))/ 1e9 - df2_s %>%  ungroup() %>% filter(scen=="Baseline (historical climate, 2020 population)" & (value_q  >= 37.5) & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T))/ 1e9)  

e <- (df2_s %>%  ungroup() %>% filter(scen=="SSP5(85), year 2050" & (value_q  >= 37.5) & pop_group=="Age group 69+" & variable=="exp") %>% dplyr::summarise(value=sum(value, na.rm=T))/ 1e9) 

write.csv(bind_rows(a, b, c, d,e), "count_pops_tmax.csv")

rm(df, df1, df2, df3, df4, df5)
gc()

#############

# decomposition analysis

#1) # elderly * CDDs, by region

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_tot <- df

#2)  # elderly * CDDs, by region (project only total pop. growth)


ssp_stats <- read.csv(paste0(stub, "/SspDb_country_data_2013-06-12.csv"))

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1","SSP2","SSP3", "SSP5"))

ssp_stats <- ssp_stats %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")


ssp_stats$gr_2010s <- (ssp_stats$X2020 - ssp_stats$X2010) / ssp_stats$X2010
ssp_stats$gr_2020s <- (ssp_stats$X2030 - ssp_stats$X2020) / ssp_stats$X2020
ssp_stats$gr_2030s <- (ssp_stats$X2040 - ssp_stats$X2030) / ssp_stats$X2030
ssp_stats$gr_2040s <- (ssp_stats$X2050 - ssp_stats$X2040) / ssp_stats$X2040

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)
r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

raster_of_iso <- fasterize::fasterize(wrld_simpl, r, field="UN", "first")

wrld_simpl_c <- dplyr::select(wrld_simpl, UN, ISO3)

ssp_stats <- merge(ssp_stats, wrld_simpl_c, by.x="REGION", by.y="ISO3")

ssp_stats <- st_as_sf(ssp_stats)

grs <- list()
it <- 0

for (scenario in sort(unique(ssp_stats$SCENARIO))){
  
  it <- it+1
  
  gr_2020s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario),c1, field="gr_2020s")
  gr_2030s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1, field="gr_2030s", "first")
  gr_2040s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1, field="gr_2040s", "first")
  
  grs[[it]] <- stack(gr_2020s, gr_2030s, gr_2040s)
  
}


r1_ssp1_2050_tot_pop <- r1_tot * (1+grs[[1]][[1]]) * (1+grs[[1]][[2]]) * (1+grs[[1]][[3]])
r1_ssp2_2050_tot_pop <- r1_tot * (1+grs[[2]][[1]]) * (1+grs[[2]][[2]]) * (1+grs[[2]][[3]])
r1_ssp3_2050_tot_pop <- r1_tot * (1+grs[[3]][[1]]) * (1+grs[[3]][[2]]) * (1+grs[[3]][[3]])
r1_ssp5_2050_tot_pop <- r1_tot * (1+grs[[4]][[1]]) * (1+grs[[4]][[2]]) * (1+grs[[4]][[3]])

df <- data.frame(values(r1_ssp2_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp2_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_245 <- df

#

df <- data.frame(values(r1_ssp5_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp5_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_585 <- df

#

df <- data.frame(values(r1_ssp1_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp1_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_126 <- df

#

df <- data.frame(values(r1_ssp3_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp3_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_370 <- df


#3)  # elderly * CDDs, by region (project only % of pop above 69)

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) /  rowSums(values(r1_ssp2_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) /  rowSums(values(r1_ssp5_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) /  rowSums(values(r1_ssp1_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_126 <- df

#


df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) /  rowSums(values(r1_ssp3_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_370 <- df

#4)  # elderly * CDDs, by region (project only climate change)

df <- data.frame(rowSums(values(r1_g)), exp = values(c2), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c3), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c4), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_126 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c5), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_370 <- df

#

# barplot of percentages of 1=100, + growth due to 2,3,4

df2_tot <- group_by(df2_tot, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_245 <- group_by(df2_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_585 <- group_by(df2_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_126 <- group_by(df2_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_370 <- group_by(df2_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_245 <- group_by(df2_pop_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_585 <- group_by(df2_pop_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_126 <- group_by(df2_pop_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_370 <- group_by(df2_pop_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_share_245 <- group_by(df2_pop_share_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_585 <- group_by(df2_pop_share_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_126 <- group_by(df2_pop_share_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_370 <- group_by(df2_pop_share_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

#

df_decompose <- bind_cols(df2_tot, df2_245, df2_585, df2_126, df2_370, df2_pop_245, df2_pop_585, df2_pop_126, df2_pop_370, df2_pop_share_245, df2_pop_share_585, df2_pop_share_126, df2_pop_share_370)

df_decompose$value...4 <- df_decompose$value...4 - df_decompose$value...2
df_decompose$value...6 <- df_decompose$value...6 - df_decompose$value...2
df_decompose$value...8 <- df_decompose$value...8 - df_decompose$value...2
df_decompose$value...10 <- df_decompose$value...10 - df_decompose$value...2
df_decompose$value...12 <- df_decompose$value...12 - df_decompose$value...2
df_decompose$value...14 <- df_decompose$value...14 - df_decompose$value...2
df_decompose$value...16 <- df_decompose$value...16 - df_decompose$value...2
df_decompose$value...18 <- df_decompose$value...18 - df_decompose$value...2
df_decompose$value...20 <- df_decompose$value...20 - df_decompose$value...2
df_decompose$value...22 <- df_decompose$value...22 - df_decompose$value...2
df_decompose$value...24 <- df_decompose$value...24 - df_decompose$value...2
df_decompose$value...26 <- df_decompose$value...26 - df_decompose$value...2

df_decompose <- reshape2::melt(df_decompose, c(1,3,5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25))

df_decompose$scenario <- NA
df_decompose$scenario[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$scenario[df_decompose$variable=="value...4"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...6"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...8"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...10"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...12"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...14"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...16"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...18"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...20"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...22"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...24"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...26"] <- "SSP370"

df_decompose$type <- NA
df_decompose$type[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$type[df_decompose$variable=="value...4"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...6"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...8"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...10"] <- "Climate change"

df_decompose$type[df_decompose$variable=="value...12"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...14"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...16"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...18"] <- "Pop. growth"

df_decompose$type[df_decompose$variable=="value...20"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...22"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...24"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...26"] <- "Aging"

df_decompose <- dplyr::select(df_decompose, -c(2:14))

colnames(df_decompose)[1] <- "reg"

#

df_decompose_a = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP245")
df_decompose_b = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP585")
df_decompose_c = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP126")
df_decompose_d = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP370")

df_decompose = df_decompose %>% filter(scenario!="Historical")
df_decompose = bind_rows(df_decompose, df_decompose_a, df_decompose_b, df_decompose_c, df_decompose_d)

df_decompose$type <- factor(df_decompose$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

deco_a <- ggplot(df_decompose)+
  theme_classic()+
  geom_col(aes(x=scenario, y=value/1e9, fill=type), show.legend = F)+
  facet_wrap(vars(reg), scales = "free")+
  ylab(expression(paste("Billion",phantom(x), PDD)))+
  xlab("Scenario")+
  scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))

# df_decompose_r <- df_decompose %>% group_by(scenario, reg...1) %>% dplyr::mutate(value=value/sum(value, na.rm=T))
# 
# df_decompose_r <- filter(df_decompose_r, scenario!="Historical")
# 
# df_decompose_r$type <- factor(df_decompose_r$type , levels = c("Climate change", "Pop. growth", "Aging"))
# 
# deco_b <- ggplot(df_decompose_r)+
#   theme_classic()+
#   geom_col(aes(x=scenario, y=value, fill=type), show.legend = F)+
#   facet_wrap(vars(reg...1), scales = "free")+
#   ylab(expression(paste("Percentage growth in",phantom(x), PDD)))+
#   xlab("Scenario")+
#   scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4]))+
# scale_y_continuous(labels = scales::label_percent()) + 
#   theme(legend.position = "none")


df_decompose_r <- df_decompose %>% group_by(scenario, type) %>% dplyr::summarise(value=sum(value, na.rm=T))

df_decompose_r$type <- factor(df_decompose_r$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

df_decompose_r$reg <- "Global"

deco_b <- ggplot(df_decompose_r)+
  theme_classic()+
  geom_col(aes(x=scenario, y=value/1e9, fill=type), show.legend = F)+
  facet_wrap(vars(reg), scales = "free")+
  ylab(expression(paste("Billion",phantom(x), PDD)))+
  xlab("Scenario")+
  scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))

###########################################
###########################################

#1) # elderly * CDDs, by region

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_tot <- df

#2)  # elderly * CDDs, by region (project only total pop. growth)


ssp_stats <- read.csv(paste0(stub, "/SspDb_country_data_2013-06-12.csv"))

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1","SSP2","SSP3", "SSP5"))

ssp_stats <- ssp_stats %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")


ssp_stats$gr_2010s <- (ssp_stats$X2020 - ssp_stats$X2010) / ssp_stats$X2010
ssp_stats$gr_2020s <- (ssp_stats$X2030 - ssp_stats$X2020) / ssp_stats$X2020
ssp_stats$gr_2030s <- (ssp_stats$X2040 - ssp_stats$X2030) / ssp_stats$X2030
ssp_stats$gr_2040s <- (ssp_stats$X2050 - ssp_stats$X2040) / ssp_stats$X2040

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)
r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

raster_of_iso <- fasterize::fasterize(wrld_simpl, r, field="UN", "first")

wrld_simpl_c <- dplyr::select(wrld_simpl, UN, ISO3)

ssp_stats <- merge(ssp_stats, wrld_simpl_c, by.x="REGION", by.y="ISO3")

ssp_stats <- st_as_sf(ssp_stats)

grs <- list()
it <- 0

for (scenario in sort(unique(ssp_stats$SCENARIO))){
  
  it <- it+1
  
  gr_2020s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario),c1_tasmax, field="gr_2020s")
  gr_2030s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_tasmax, field="gr_2030s", "first")
  gr_2040s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_tasmax, field="gr_2040s", "first")
  
  grs[[it]] <- stack(gr_2020s, gr_2030s, gr_2040s)
  
}


r1_ssp1_2050_tot_pop <- r1_tot * (1+grs[[1]][[1]]) * (1+grs[[1]][[2]]) * (1+grs[[1]][[3]])
r1_ssp2_2050_tot_pop <- r1_tot * (1+grs[[2]][[1]]) * (1+grs[[2]][[2]]) * (1+grs[[2]][[3]])
r1_ssp3_2050_tot_pop <- r1_tot * (1+grs[[3]][[1]]) * (1+grs[[3]][[2]]) * (1+grs[[3]][[3]])
r1_ssp5_2050_tot_pop <- r1_tot * (1+grs[[4]][[1]]) * (1+grs[[4]][[2]]) * (1+grs[[4]][[3]])

df <- data.frame(values(r1_ssp2_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp2_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_245 <- df

#

df <- data.frame(values(r1_ssp5_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp5_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_585 <- df

#

df <- data.frame(values(r1_ssp1_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp1_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_126 <- df

#

df <- data.frame(values(r1_ssp3_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp3_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_370 <- df


#3)  # elderly * CDDs, by region (project only % of pop above 69)

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) /  rowSums(values(r1_ssp2_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) /  rowSums(values(r1_ssp5_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) /  rowSums(values(r1_ssp1_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_126 <- df

#


df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) /  rowSums(values(r1_ssp3_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_370 <- df

#4)  # elderly * CDDs, by region (project only climate change)

df <- data.frame(rowSums(values(r1_g)), exp = values(c2_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c3_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c4_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_126 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c5_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_370 <- df

#

# barplot of percentages of 1=100, + growth due to 2,3,4

df2_tot <- group_by(df2_tot, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_245 <- group_by(df2_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_585 <- group_by(df2_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_126 <- group_by(df2_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_370 <- group_by(df2_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_245 <- group_by(df2_pop_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_585 <- group_by(df2_pop_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_126 <- group_by(df2_pop_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_370 <- group_by(df2_pop_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_share_245 <- group_by(df2_pop_share_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_585 <- group_by(df2_pop_share_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_126 <- group_by(df2_pop_share_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_370 <- group_by(df2_pop_share_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

#

df_decompose <- bind_cols(df2_tot, df2_245, df2_585, df2_126, df2_370, df2_pop_245, df2_pop_585, df2_pop_126, df2_pop_370, df2_pop_share_245, df2_pop_share_585, df2_pop_share_126, df2_pop_share_370)

df_decompose$value...4 <- df_decompose$value...4 - df_decompose$value...2
df_decompose$value...6 <- df_decompose$value...6 - df_decompose$value...2
df_decompose$value...8 <- df_decompose$value...8 - df_decompose$value...2
df_decompose$value...10 <- df_decompose$value...10 - df_decompose$value...2
df_decompose$value...12 <- df_decompose$value...12 - df_decompose$value...2
df_decompose$value...14 <- df_decompose$value...14 - df_decompose$value...2
df_decompose$value...16 <- df_decompose$value...16 - df_decompose$value...2
df_decompose$value...18 <- df_decompose$value...18 - df_decompose$value...2
df_decompose$value...20 <- df_decompose$value...20 - df_decompose$value...2
df_decompose$value...22 <- df_decompose$value...22 - df_decompose$value...2
df_decompose$value...24 <- df_decompose$value...24 - df_decompose$value...2
df_decompose$value...26 <- df_decompose$value...26 - df_decompose$value...2

df_decompose <- reshape2::melt(df_decompose, c(1,3,5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25))

df_decompose$scenario <- NA
df_decompose$scenario[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$scenario[df_decompose$variable=="value...4"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...6"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...8"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...10"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...12"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...14"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...16"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...18"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...20"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...22"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...24"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...26"] <- "SSP370"

df_decompose$type <- NA
df_decompose$type[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$type[df_decompose$variable=="value...4"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...6"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...8"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...10"] <- "Climate change"

df_decompose$type[df_decompose$variable=="value...12"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...14"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...16"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...18"] <- "Pop. growth"

df_decompose$type[df_decompose$variable=="value...20"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...22"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...24"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...26"] <- "Aging"

df_decompose <- dplyr::select(df_decompose, -c(2:14))

colnames(df_decompose)[1] <- "reg"

#

df_decompose_a = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP245")
df_decompose_b = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP585")
df_decompose_c = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP126")
df_decompose_d = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP370")

df_decompose = df_decompose %>% filter(scenario!="Historical")
df_decompose = bind_rows(df_decompose, df_decompose_a, df_decompose_b, df_decompose_c, df_decompose_d)

df_decompose$type <- factor(df_decompose$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

deco_c <- ggplot(df_decompose)+
  theme_classic()+
  geom_col(aes(x=scenario, y=value/1e9, fill=type), show.legend = F)+
  facet_wrap(vars(reg), scales = "free")+
  ylab(expression(paste("Billion",phantom(x), PD95)))+
  xlab("Scenario")+
  scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))

# df_decompose_r <- df_decompose %>% group_by(scenario, reg...1) %>% dplyr::mutate(value=value/sum(value, na.rm=T))
# 
# df_decompose_r <- filter(df_decompose_r, scenario!="Historical")
# 
# df_decompose_r$type <- factor(df_decompose_r$type , levels = c("Climate change", "Pop. growth", "Aging"))
# 
# deco_b <- ggplot(df_decompose_r)+
#   theme_classic()+
#   geom_col(aes(x=scenario, y=value, fill=type), show.legend = F)+
#   facet_wrap(vars(reg...1), scales = "free")+
#   ylab(expression(paste("Percentage growth in",phantom(x), PDD)))+
#   xlab("Scenario")+
#   scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4]))+
# scale_y_continuous(labels = scales::label_percent()) + 
#   theme(legend.position = "none")


df_decompose_r <- df_decompose %>% group_by(scenario, type) %>% dplyr::summarise(value=sum(value, na.rm=T))

df_decompose_r$type <- factor(df_decompose_r$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

df_decompose_r$reg <- "Global"

deco_d <- ggplot(df_decompose_r)+
  theme_classic()+
  geom_col(aes(x=scenario, y=value/1e9, fill=type))+
  facet_wrap(vars(reg), scales = "free")+
  ylab(expression(paste("Billion",phantom(x), PD95)))+
  xlab("Scenario")+
  scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))

###########################################
###########################################

#1) # elderly * CDDs, by region

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_tot <- df

#2)  # elderly * CDDs, by region (project only total pop. growth)


ssp_stats <- read.csv(paste0(stub, "/SspDb_country_data_2013-06-12.csv"))

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1","SSP2","SSP3", "SSP5"))

ssp_stats <- ssp_stats %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")


ssp_stats$gr_2010s <- (ssp_stats$X2020 - ssp_stats$X2010) / ssp_stats$X2010
ssp_stats$gr_2020s <- (ssp_stats$X2030 - ssp_stats$X2020) / ssp_stats$X2020
ssp_stats$gr_2030s <- (ssp_stats$X2040 - ssp_stats$X2030) / ssp_stats$X2030
ssp_stats$gr_2040s <- (ssp_stats$X2050 - ssp_stats$X2040) / ssp_stats$X2040

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)
r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

raster_of_iso <- fasterize::fasterize(wrld_simpl, r, field="UN", "first")

wrld_simpl_c <- dplyr::select(wrld_simpl, UN, ISO3)

ssp_stats <- merge(ssp_stats, wrld_simpl_c, by.x="REGION", by.y="ISO3")

ssp_stats <- st_as_sf(ssp_stats)

grs <- list()
it <- 0

for (scenario in sort(unique(ssp_stats$SCENARIO))){
  
  it <- it+1
  
  gr_2020s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario),c1_hotday, field="gr_2020s")
  gr_2030s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_hotday, field="gr_2030s", "first")
  gr_2040s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_hotday, field="gr_2040s", "first")
  
  grs[[it]] <- stack(gr_2020s, gr_2030s, gr_2040s)
  
}


r1_ssp1_2050_tot_pop <- r1_tot * (1+grs[[1]][[1]]) * (1+grs[[1]][[2]]) * (1+grs[[1]][[3]])
r1_ssp2_2050_tot_pop <- r1_tot * (1+grs[[2]][[1]]) * (1+grs[[2]][[2]]) * (1+grs[[2]][[3]])
r1_ssp3_2050_tot_pop <- r1_tot * (1+grs[[3]][[1]]) * (1+grs[[3]][[2]]) * (1+grs[[3]][[3]])
r1_ssp5_2050_tot_pop <- r1_tot * (1+grs[[4]][[1]]) * (1+grs[[4]][[2]]) * (1+grs[[4]][[3]])

df <- data.frame(values(r1_ssp2_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp2_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_245 <- df

#

df <- data.frame(values(r1_ssp5_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp5_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_585 <- df

#

df <- data.frame(values(r1_ssp1_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp1_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_126 <- df

#

df <- data.frame(values(r1_ssp3_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp3_2050_tot_pop.*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"


df <- dplyr::select(df, -starts_with("values"))
df2_pop_370 <- df


#3)  # elderly * CDDs, by region (project only % of pop above 69)

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) /  rowSums(values(r1_ssp2_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) /  rowSums(values(r1_ssp5_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) /  rowSums(values(r1_ssp1_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_126 <- df

#


df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) /  rowSums(values(r1_ssp3_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_pop_share_370 <- df

#4)  # elderly * CDDs, by region (project only climate change)

df <- data.frame(rowSums(values(r1_g)), exp = values(c2_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c3_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c4_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_126 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c5_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg[df$reg==1] <- "Africa"
df$reg[df$reg==2] <- NA
df$reg[df$reg==3] <- "Asia"
df$reg[df$reg==4] <- "Oceania"
df$reg[df$reg==5] <- "Europe"
df$reg[df$reg==6] <- "North America"
df$reg[df$reg==7] <- "South America"

df2_370 <- df

#

# barplot of percentages of 1=100, + growth due to 2,3,4

df2_tot <- group_by(df2_tot, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_245 <- group_by(df2_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_585 <- group_by(df2_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_126 <- group_by(df2_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_370 <- group_by(df2_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_245 <- group_by(df2_pop_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_585 <- group_by(df2_pop_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_126 <- group_by(df2_pop_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_370 <- group_by(df2_pop_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_share_245 <- group_by(df2_pop_share_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_585 <- group_by(df2_pop_share_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_126 <- group_by(df2_pop_share_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_370 <- group_by(df2_pop_share_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

#

df_decompose <- bind_cols(df2_tot, df2_245, df2_585, df2_126, df2_370, df2_pop_245, df2_pop_585, df2_pop_126, df2_pop_370, df2_pop_share_245, df2_pop_share_585, df2_pop_share_126, df2_pop_share_370)

df_decompose$value...4 <- df_decompose$value...4 - df_decompose$value...2
df_decompose$value...6 <- df_decompose$value...6 - df_decompose$value...2
df_decompose$value...8 <- df_decompose$value...8 - df_decompose$value...2
df_decompose$value...10 <- df_decompose$value...10 - df_decompose$value...2
df_decompose$value...12 <- df_decompose$value...12 - df_decompose$value...2
df_decompose$value...14 <- df_decompose$value...14 - df_decompose$value...2
df_decompose$value...16 <- df_decompose$value...16 - df_decompose$value...2
df_decompose$value...18 <- df_decompose$value...18 - df_decompose$value...2
df_decompose$value...20 <- df_decompose$value...20 - df_decompose$value...2
df_decompose$value...22 <- df_decompose$value...22 - df_decompose$value...2
df_decompose$value...24 <- df_decompose$value...24 - df_decompose$value...2
df_decompose$value...26 <- df_decompose$value...26 - df_decompose$value...2

df_decompose <- reshape2::melt(df_decompose, c(1,3,5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25))

df_decompose$scenario <- NA
df_decompose$scenario[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$scenario[df_decompose$variable=="value...4"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...6"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...8"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...10"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...12"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...14"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...16"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...18"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...20"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...22"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...24"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...26"] <- "SSP370"

df_decompose$type <- NA
df_decompose$type[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$type[df_decompose$variable=="value...4"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...6"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...8"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...10"] <- "Climate change"

df_decompose$type[df_decompose$variable=="value...12"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...14"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...16"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...18"] <- "Pop. growth"

df_decompose$type[df_decompose$variable=="value...20"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...22"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...24"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...26"] <- "Aging"

df_decompose <- dplyr::select(df_decompose, -c(2:14))

colnames(df_decompose)[1] <- "reg"

#

df_decompose_a = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP245")
df_decompose_b = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP585")
df_decompose_c = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP126")
df_decompose_d = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP370")

df_decompose = df_decompose %>% filter(scenario!="Historical")
df_decompose = bind_rows(df_decompose, df_decompose_a, df_decompose_b, df_decompose_c, df_decompose_d)

df_decompose$type <- factor(df_decompose$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

deco_e <- ggplot(df_decompose)+
  theme_classic()+
  geom_col(aes(x=scenario, y=value/1e9, fill=type), show.legend = F)+
  facet_wrap(vars(reg), scales = "free")+
  ylab(expression(paste("Billion",phantom(x), PHD)))+
  xlab("Scenario")+
  scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))

# df_decompose_r <- df_decompose %>% group_by(scenario, reg...1) %>% dplyr::mutate(value=value/sum(value, na.rm=T))
# 
# df_decompose_r <- filter(df_decompose_r, scenario!="Historical")
# 
# df_decompose_r$type <- factor(df_decompose_r$type , levels = c("Climate change", "Pop. growth", "Aging"))
# 
# deco_b <- ggplot(df_decompose_r)+
#   theme_classic()+
#   geom_col(aes(x=scenario, y=value, fill=type), show.legend = F)+
#   facet_wrap(vars(reg...1), scales = "free")+
#   ylab(expression(paste("Percentage growth in",phantom(x), PDD)))+
#   xlab("Scenario")+
#   scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4]))+
# scale_y_continuous(labels = scales::label_percent()) + 
#   theme(legend.position = "none")


df_decompose_r <- df_decompose %>% group_by(scenario, type) %>% dplyr::summarise(value=sum(value, na.rm=T))

df_decompose_r$type <- factor(df_decompose_r$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

df_decompose_r$reg <- "Global"

deco_f <- ggplot(df_decompose_r)+
  theme_classic()+
  geom_col(aes(x=scenario, y=value/1e9, fill=type))+
  facet_wrap(vars(reg), scales = "free")+
  ylab(expression(paste("Billion",phantom(x), PHD)))+
  xlab("Scenario")+
  scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))



deco <- deco_a + deco_b + deco_e + deco_f + deco_c + deco_d + plot_annotation(tag_levels = 'A') + plot_layout(ncol=2, widths = c(1,0.5,1,0.5,1,0.5), guides = "collect") & theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


###########################

ggsave("decomposition.pdf", scale=1.5, width = 6, height = 7)
ggsave("decomposition.png", scale=1.5, width = 6, height = 7)

###############

# decomposition analysis AT THE COUNTRY LEVEL


regions <- st_as_sf(rworldmap::countriesLow)
regions <- filter(regions, regions$POP_EST>500000)
regions$ISO3_n <- as.numeric(regions$ISO3)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="ISO3")
regions <- dplyr::select(regions, ISO3, ISO3_n)
regions$geometry<-NULL


#1) # elderly * CDDs, by region

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_tot <- df

#2)  # elderly * CDDs, by region (project only total pop. growth)


ssp_stats <- read.csv(paste0(stub, "/SspDb_country_data_2013-06-12.csv"))

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1","SSP2","SSP3", "SSP5"))

ssp_stats <- ssp_stats %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")


ssp_stats$gr_2010s <- (ssp_stats$X2020 - ssp_stats$X2010) / ssp_stats$X2010
ssp_stats$gr_2020s <- (ssp_stats$X2030 - ssp_stats$X2020) / ssp_stats$X2020
ssp_stats$gr_2030s <- (ssp_stats$X2040 - ssp_stats$X2030) / ssp_stats$X2030
ssp_stats$gr_2040s <- (ssp_stats$X2050 - ssp_stats$X2040) / ssp_stats$X2040

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)
r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

raster_of_iso <- fasterize::fasterize(wrld_simpl, r, field="UN", "first")

wrld_simpl_c <- dplyr::select(wrld_simpl, UN, ISO3)

ssp_stats <- merge(ssp_stats, wrld_simpl_c, by.x="REGION", by.y="ISO3")

ssp_stats <- st_as_sf(ssp_stats)

grs <- list()
it <- 0

for (scenario in sort(unique(ssp_stats$SCENARIO))){
  
  it <- it+1
  
  gr_2020s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario),c1, field="gr_2020s")
  gr_2030s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1, field="gr_2030s", "first")
  gr_2040s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1, field="gr_2040s", "first")
  
  grs[[it]] <- stack(gr_2020s, gr_2030s, gr_2040s)
  
}


r1_ssp1_2050_tot_pop <- r1_tot * (1+grs[[1]][[1]]) * (1+grs[[1]][[2]]) * (1+grs[[1]][[3]])
r1_ssp2_2050_tot_pop <- r1_tot * (1+grs[[2]][[1]]) * (1+grs[[2]][[2]]) * (1+grs[[2]][[3]])
r1_ssp3_2050_tot_pop <- r1_tot * (1+grs[[3]][[1]]) * (1+grs[[3]][[2]]) * (1+grs[[3]][[3]])
r1_ssp5_2050_tot_pop <- r1_tot * (1+grs[[4]][[1]]) * (1+grs[[4]][[2]]) * (1+grs[[4]][[3]])

df <- data.frame(values(r1_ssp2_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp2_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_245 <- df

#

df <- data.frame(values(r1_ssp5_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp5_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_585 <- df

#

df <- data.frame(values(r1_ssp1_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp1_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_126 <- df

#

df <- data.frame(values(r1_ssp3_2050_tot_pop), exp = values(c1), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp3_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_370 <- df


#3)  # elderly * CDDs, by region (project only % of pop above 69)

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) /  rowSums(values(r1_ssp2_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) /  rowSums(values(r1_ssp5_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) /  rowSums(values(r1_ssp1_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_126 <- df

#


df <- data.frame(rowSums(values(r1_g)), exp = values(c1), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) /  rowSums(values(r1_ssp3_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_370 <- df

#4)  # elderly * CDDs, by region (project only climate change)

df <- data.frame(rowSums(values(r1_g)), exp = values(c2), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c3), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c4), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_126 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c5), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_370 <- df

#

# barplot of percentages of 1=100, + growth due to 2,3,4

df2_tot <- group_by(df2_tot, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_245 <- group_by(df2_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_585 <- group_by(df2_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_126 <- group_by(df2_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_370 <- group_by(df2_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_245 <- group_by(df2_pop_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_585 <- group_by(df2_pop_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_126 <- group_by(df2_pop_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_370 <- group_by(df2_pop_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_share_245 <- group_by(df2_pop_share_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_585 <- group_by(df2_pop_share_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_126 <- group_by(df2_pop_share_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_370 <- group_by(df2_pop_share_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

#

df2_tot <- filter(df2_tot, df2_tot$reg %in% df2_pop_126$reg)
df2_245 <- filter(df2_245, df2_245$reg %in% df2_pop_126$reg)
df2_585 <- filter(df2_585, df2_585$reg %in% df2_pop_126$reg)
df2_370 <- filter(df2_370, df2_370$reg %in% df2_pop_126$reg)
df2_126 <- filter(df2_126, df2_126$reg %in% df2_pop_126$reg)

df_decompose <- bind_cols(df2_tot, df2_245, df2_585, df2_126, df2_370, df2_pop_245, df2_pop_585, df2_pop_126, df2_pop_370, df2_pop_share_245, df2_pop_share_585, df2_pop_share_126, df2_pop_share_370)

df_decompose$value...4 <- df_decompose$value...4 - df_decompose$value...2
df_decompose$value...6 <- df_decompose$value...6 - df_decompose$value...2
df_decompose$value...8 <- df_decompose$value...8 - df_decompose$value...2
df_decompose$value...10 <- df_decompose$value...10 - df_decompose$value...2
df_decompose$value...12 <- df_decompose$value...12 - df_decompose$value...2
df_decompose$value...14 <- df_decompose$value...14 - df_decompose$value...2
df_decompose$value...16 <- df_decompose$value...16 - df_decompose$value...2
df_decompose$value...18 <- df_decompose$value...18 - df_decompose$value...2
df_decompose$value...20 <- df_decompose$value...20 - df_decompose$value...2
df_decompose$value...22 <- df_decompose$value...22 - df_decompose$value...2
df_decompose$value...24 <- df_decompose$value...24 - df_decompose$value...2
df_decompose$value...26 <- df_decompose$value...26 - df_decompose$value...2

df_decompose <- reshape2::melt(df_decompose, c(1,3,5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25))

df_decompose$scenario <- NA
df_decompose$scenario[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$scenario[df_decompose$variable=="value...4"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...6"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...8"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...10"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...12"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...14"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...16"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...18"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...20"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...22"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...24"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...26"] <- "SSP370"

df_decompose$type <- NA
df_decompose$type[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$type[df_decompose$variable=="value...4"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...6"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...8"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...10"] <- "Climate change"

df_decompose$type[df_decompose$variable=="value...12"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...14"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...16"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...18"] <- "Pop. growth"

df_decompose$type[df_decompose$variable=="value...20"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...22"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...24"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...26"] <- "Aging"

df_decompose <- dplyr::select(df_decompose, -c(2:14))

colnames(df_decompose)[1] <- "reg"

#

df_decompose_a = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP245")
df_decompose_b = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP585")
df_decompose_c = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP126")
df_decompose_d = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP370")

df_decompose = df_decompose %>% filter(scenario!="Historical")
df_decompose = bind_rows(df_decompose, df_decompose_a, df_decompose_b, df_decompose_c, df_decompose_d)

df_decompose$type <- factor(df_decompose$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

df_decompose_r <- df_decompose

for(i in 1:5){
  ggplot(df_decompose_r)+
    theme_classic()+
    geom_col(aes(x=scenario, y=value/1e9, fill=type))+
    facet_wrap_paginate(vars(reg), scales = "free", ncol = 6, nrow = 6, page = i)+
    ylab(expression(paste("Growth in",phantom(x),phantom(x), PDD, ", (billion)")))+
    xlab("Scenario")+
    scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))+
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0("decomposition_country_", i, ".pdf"), scale=2, width = 8, height = 6)
}

####


regions <- st_as_sf(rworldmap::countriesLow)
regions <- filter(regions, regions$POP_EST>500000)
regions$ISO3_n <- as.numeric(regions$ISO3)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="ISO3")
regions <- dplyr::select(regions, ISO3, ISO3_n)
regions$geometry<-NULL


#1) # elderly * CDDs, by region

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_tot <- df

#2)  # elderly * CDDs, by region (project only total pop. growth)


ssp_stats <- read.csv(paste0(stub, "/SspDb_country_data_2013-06-12.csv"))

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1","SSP2","SSP3", "SSP5"))

ssp_stats <- ssp_stats %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")


ssp_stats$gr_2010s <- (ssp_stats$X2020 - ssp_stats$X2010) / ssp_stats$X2010
ssp_stats$gr_2020s <- (ssp_stats$X2030 - ssp_stats$X2020) / ssp_stats$X2020
ssp_stats$gr_2030s <- (ssp_stats$X2040 - ssp_stats$X2030) / ssp_stats$X2030
ssp_stats$gr_2040s <- (ssp_stats$X2050 - ssp_stats$X2040) / ssp_stats$X2040

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)
r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

raster_of_iso <- fasterize::fasterize(wrld_simpl, r, field="UN", "first")

wrld_simpl_c <- dplyr::select(wrld_simpl, UN, ISO3)

ssp_stats <- merge(ssp_stats, wrld_simpl_c, by.x="REGION", by.y="ISO3")

ssp_stats <- st_as_sf(ssp_stats)

grs <- list()
it <- 0

for (scenario in sort(unique(ssp_stats$SCENARIO))){
  
  it <- it+1
  
  gr_2020s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario),c1_hotday, field="gr_2020s")
  gr_2030s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_hotday, field="gr_2030s", "first")
  gr_2040s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_hotday, field="gr_2040s", "first")
  
  grs[[it]] <- stack(gr_2020s, gr_2030s, gr_2040s)
  
}


r1_ssp1_2050_tot_pop <- r1_tot * (1+grs[[1]][[1]]) * (1+grs[[1]][[2]]) * (1+grs[[1]][[3]])
r1_ssp2_2050_tot_pop <- r1_tot * (1+grs[[2]][[1]]) * (1+grs[[2]][[2]]) * (1+grs[[2]][[3]])
r1_ssp3_2050_tot_pop <- r1_tot * (1+grs[[3]][[1]]) * (1+grs[[3]][[2]]) * (1+grs[[3]][[3]])
r1_ssp5_2050_tot_pop <- r1_tot * (1+grs[[4]][[1]]) * (1+grs[[4]][[2]]) * (1+grs[[4]][[3]])

df <- data.frame(values(r1_ssp2_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp2_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_245 <- df

#

df <- data.frame(values(r1_ssp5_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp5_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_585 <- df

#

df <- data.frame(values(r1_ssp1_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp1_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_126 <- df

#

df <- data.frame(values(r1_ssp3_2050_tot_pop), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp3_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_370 <- df


#3)  # elderly * CDDs, by region (project only % of pop above 69)

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) /  rowSums(values(r1_ssp2_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) /  rowSums(values(r1_ssp5_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) /  rowSums(values(r1_ssp1_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_126 <- df

#


df <- data.frame(rowSums(values(r1_g)), exp = values(c1_hotday), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) /  rowSums(values(r1_ssp3_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_370 <- df

#4)  # elderly * CDDs, by region (project only climate change)

df <- data.frame(rowSums(values(r1_g)), exp = values(c2_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c3_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c4_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_126 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c5_hotday), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_370 <- df

#

# barplot of percentages of 1=100, + growth due to 2,3,4

df2_tot <- group_by(df2_tot, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_245 <- group_by(df2_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_585 <- group_by(df2_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_126 <- group_by(df2_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_370 <- group_by(df2_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_245 <- group_by(df2_pop_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_585 <- group_by(df2_pop_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_126 <- group_by(df2_pop_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_370 <- group_by(df2_pop_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_share_245 <- group_by(df2_pop_share_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_585 <- group_by(df2_pop_share_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_126 <- group_by(df2_pop_share_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_370 <- group_by(df2_pop_share_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

#

df2_tot <- filter(df2_tot, df2_tot$reg %in% df2_pop_126$reg)
df2_245 <- filter(df2_245, df2_245$reg %in% df2_pop_126$reg)
df2_585 <- filter(df2_585, df2_585$reg %in% df2_pop_126$reg)
df2_370 <- filter(df2_370, df2_370$reg %in% df2_pop_126$reg)
df2_126 <- filter(df2_126, df2_126$reg %in% df2_pop_126$reg)

df_decompose <- bind_cols(df2_tot, df2_245, df2_585, df2_126, df2_370, df2_pop_245, df2_pop_585, df2_pop_126, df2_pop_370, df2_pop_share_245, df2_pop_share_585, df2_pop_share_126, df2_pop_share_370)

df_decompose$value...4 <- df_decompose$value...4 - df_decompose$value...2
df_decompose$value...6 <- df_decompose$value...6 - df_decompose$value...2
df_decompose$value...8 <- df_decompose$value...8 - df_decompose$value...2
df_decompose$value...10 <- df_decompose$value...10 - df_decompose$value...2
df_decompose$value...12 <- df_decompose$value...12 - df_decompose$value...2
df_decompose$value...14 <- df_decompose$value...14 - df_decompose$value...2
df_decompose$value...16 <- df_decompose$value...16 - df_decompose$value...2
df_decompose$value...18 <- df_decompose$value...18 - df_decompose$value...2
df_decompose$value...20 <- df_decompose$value...20 - df_decompose$value...2
df_decompose$value...22 <- df_decompose$value...22 - df_decompose$value...2
df_decompose$value...24 <- df_decompose$value...24 - df_decompose$value...2
df_decompose$value...26 <- df_decompose$value...26 - df_decompose$value...2

df_decompose <- reshape2::melt(df_decompose, c(1,3,5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25))

df_decompose$scenario <- NA
df_decompose$scenario[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$scenario[df_decompose$variable=="value...4"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...6"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...8"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...10"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...12"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...14"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...16"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...18"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...20"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...22"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...24"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...26"] <- "SSP370"

df_decompose$type <- NA
df_decompose$type[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$type[df_decompose$variable=="value...4"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...6"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...8"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...10"] <- "Climate change"

df_decompose$type[df_decompose$variable=="value...12"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...14"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...16"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...18"] <- "Pop. growth"

df_decompose$type[df_decompose$variable=="value...20"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...22"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...24"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...26"] <- "Aging"

df_decompose <- dplyr::select(df_decompose, -c(2:14))

colnames(df_decompose)[1] <- "reg"

#

df_decompose_a = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP245")
df_decompose_b = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP585")
df_decompose_c = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP126")
df_decompose_d = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP370")

df_decompose = df_decompose %>% filter(scenario!="Historical")
df_decompose = bind_rows(df_decompose, df_decompose_a, df_decompose_b, df_decompose_c, df_decompose_d)

df_decompose$type <- factor(df_decompose$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

df_decompose_r <- df_decompose

for(i in 1:5){
  ggplot(df_decompose_r)+
    theme_classic()+
    geom_col(aes(x=scenario, y=value/1e6, fill=type))+
    facet_wrap_paginate(vars(reg), scales = "free", ncol = 6, nrow = 6, page = i)+
    ylab(expression(paste("Growth in",phantom(x),phantom(x), PHD, ", (million)")))+
    xlab("Scenario")+
    scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))+
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0("decomposition_country_", i, "_hotdays.pdf"), scale=2, width = 8, height = 6)
}

###

####


regions <- st_as_sf(rworldmap::countriesLow)
regions <- filter(regions, regions$POP_EST>500000)
regions$ISO3_n <- as.numeric(regions$ISO3)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="ISO3")
regions <- dplyr::select(regions, ISO3, ISO3_n)
regions$geometry<-NULL


#1) # elderly * CDDs, by region

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_tot <- df

#2)  # elderly * CDDs, by region (project only total pop. growth)


ssp_stats <- read.csv(paste0(stub, "/SspDb_country_data_2013-06-12.csv"))

ssp_stats$SCENARIO <- substr(ssp_stats$SCENARIO, 1, 4)

ssp_stats <- filter(ssp_stats, SCENARIO %in% c("SSP1","SSP2","SSP3", "SSP5"))

ssp_stats <- ssp_stats %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

ssp_stats <- dplyr::select(ssp_stats, 1,2,3,4,18, 20, 22, 24, 26)

ssp_stats <-filter(ssp_stats, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))

ssp_stats <- group_by(ssp_stats, SCENARIO, REGION) %>% dplyr::summarise_if(., is.numeric, "sum")


ssp_stats$gr_2010s <- (ssp_stats$X2020 - ssp_stats$X2010) / ssp_stats$X2010
ssp_stats$gr_2020s <- (ssp_stats$X2030 - ssp_stats$X2020) / ssp_stats$X2020
ssp_stats$gr_2030s <- (ssp_stats$X2040 - ssp_stats$X2030) / ssp_stats$X2030
ssp_stats$gr_2040s <- (ssp_stats$X2050 - ssp_stats$X2040) / ssp_stats$X2040

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)
r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

raster_of_iso <- fasterize::fasterize(wrld_simpl, r, field="UN", "first")

wrld_simpl_c <- dplyr::select(wrld_simpl, UN, ISO3)

ssp_stats <- merge(ssp_stats, wrld_simpl_c, by.x="REGION", by.y="ISO3")

ssp_stats <- st_as_sf(ssp_stats)

grs <- list()
it <- 0

for (scenario in sort(unique(ssp_stats$SCENARIO))){
  
  it <- it+1
  
  gr_2020s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario),c1_tasmax, field="gr_2020s")
  gr_2030s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_tasmax, field="gr_2030s", "first")
  gr_2040s <- fasterize::fasterize(ssp_stats %>% filter(SCENARIO== scenario), c1_tasmax, field="gr_2040s", "first")
  
  grs[[it]] <- stack(gr_2020s, gr_2030s, gr_2040s)
  
}


r1_ssp1_2050_tot_pop <- r1_tot * (1+grs[[1]][[1]]) * (1+grs[[1]][[2]]) * (1+grs[[1]][[3]])
r1_ssp2_2050_tot_pop <- r1_tot * (1+grs[[2]][[1]]) * (1+grs[[2]][[2]]) * (1+grs[[2]][[3]])
r1_ssp3_2050_tot_pop <- r1_tot * (1+grs[[3]][[1]]) * (1+grs[[3]][[2]]) * (1+grs[[3]][[3]])
r1_ssp5_2050_tot_pop <- r1_tot * (1+grs[[4]][[1]]) * (1+grs[[4]][[2]]) * (1+grs[[4]][[3]])

df <- data.frame(values(r1_ssp2_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp2_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_245 <- df

#

df <- data.frame(values(r1_ssp5_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp5_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_585 <- df

#

df <- data.frame(values(r1_ssp1_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp1_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_126 <- df

#

df <- data.frame(values(r1_ssp3_2050_tot_pop), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=exp*values.r1_ssp3_2050_tot_pop.*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df <- dplyr::select(df, -starts_with("values"))
df2_pop_370 <- df


#3)  # elderly * CDDs, by region (project only % of pop above 69)

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) /  rowSums(values(r1_ssp2_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) /  rowSums(values(r1_ssp5_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) /  rowSums(values(r1_ssp1_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_126 <- df

#


df <- data.frame(rowSums(values(r1_g)), exp = values(c1_tasmax), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) /  rowSums(values(r1_ssp3_2050)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_pop_share_370 <- df

#4)  # elderly * CDDs, by region (project only climate change)

df <- data.frame(rowSums(values(r1_g)), exp = values(c2_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_245 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c3_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_585 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c4_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_126 <- df

#

df <- data.frame(rowSums(values(r1_g)), exp = values(c5_tasmax), reg = values(regions_r), share=values(r1_g[[3]]) /  rowSums(values(r1_g)), xy=T)
df <- na.omit(df)

df <- df %>% dplyr::group_by(reg) %>% mutate(tot=rowSums.values.r1_g..*exp*share)

df$reg <- regions$ISO3[match(df$reg, regions$ISO3_n)]

df2_370 <- df

#

# barplot of percentages of 1=100, + growth due to 2,3,4

df2_tot <- group_by(df2_tot, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_245 <- group_by(df2_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_585 <- group_by(df2_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_126 <- group_by(df2_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_370 <- group_by(df2_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_245 <- group_by(df2_pop_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_585 <- group_by(df2_pop_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_126 <- group_by(df2_pop_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_370 <- group_by(df2_pop_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

df2_pop_share_245 <- group_by(df2_pop_share_245, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_585 <- group_by(df2_pop_share_585, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_126 <- group_by(df2_pop_share_126, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))
df2_pop_share_370 <- group_by(df2_pop_share_370, reg) %>% dplyr::summarise(value=sum(tot, na.rm=T))

#

df2_tot <- filter(df2_tot, df2_tot$reg %in% df2_pop_126$reg)
df2_245 <- filter(df2_245, df2_245$reg %in% df2_pop_126$reg)
df2_585 <- filter(df2_585, df2_585$reg %in% df2_pop_126$reg)
df2_370 <- filter(df2_370, df2_370$reg %in% df2_pop_126$reg)
df2_126 <- filter(df2_126, df2_126$reg %in% df2_pop_126$reg)

df_decompose <- bind_cols(df2_tot, df2_245, df2_585, df2_126, df2_370, df2_pop_245, df2_pop_585, df2_pop_126, df2_pop_370, df2_pop_share_245, df2_pop_share_585, df2_pop_share_126, df2_pop_share_370)

df_decompose$value...4 <- df_decompose$value...4 - df_decompose$value...2
df_decompose$value...6 <- df_decompose$value...6 - df_decompose$value...2
df_decompose$value...8 <- df_decompose$value...8 - df_decompose$value...2
df_decompose$value...10 <- df_decompose$value...10 - df_decompose$value...2
df_decompose$value...12 <- df_decompose$value...12 - df_decompose$value...2
df_decompose$value...14 <- df_decompose$value...14 - df_decompose$value...2
df_decompose$value...16 <- df_decompose$value...16 - df_decompose$value...2
df_decompose$value...18 <- df_decompose$value...18 - df_decompose$value...2
df_decompose$value...20 <- df_decompose$value...20 - df_decompose$value...2
df_decompose$value...22 <- df_decompose$value...22 - df_decompose$value...2
df_decompose$value...24 <- df_decompose$value...24 - df_decompose$value...2
df_decompose$value...26 <- df_decompose$value...26 - df_decompose$value...2

df_decompose <- reshape2::melt(df_decompose, c(1,3,5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25))

df_decompose$scenario <- NA
df_decompose$scenario[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$scenario[df_decompose$variable=="value...4"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...6"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...8"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...10"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...12"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...14"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...16"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...18"] <- "SSP370"

df_decompose$scenario[df_decompose$variable=="value...20"] <- "SSP245"
df_decompose$scenario[df_decompose$variable=="value...22"] <- "SSP585"
df_decompose$scenario[df_decompose$variable=="value...24"] <- "SSP126"
df_decompose$scenario[df_decompose$variable=="value...26"] <- "SSP370"

df_decompose$type <- NA
df_decompose$type[df_decompose$variable=="value...2"] <- "Historical"
df_decompose$type[df_decompose$variable=="value...4"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...6"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...8"] <- "Climate change"
df_decompose$type[df_decompose$variable=="value...10"] <- "Climate change"

df_decompose$type[df_decompose$variable=="value...12"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...14"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...16"] <- "Pop. growth"
df_decompose$type[df_decompose$variable=="value...18"] <- "Pop. growth"

df_decompose$type[df_decompose$variable=="value...20"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...22"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...24"] <- "Aging"
df_decompose$type[df_decompose$variable=="value...26"] <- "Aging"

df_decompose <- dplyr::select(df_decompose, -c(2:14))

colnames(df_decompose)[1] <- "reg"

#

df_decompose_a = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP245")
df_decompose_b = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP585")
df_decompose_c = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP126")
df_decompose_d = df_decompose %>% filter(scenario=="Historical") %>%  mutate(scenario = "SSP370")

df_decompose = df_decompose %>% filter(scenario!="Historical")
df_decompose = bind_rows(df_decompose, df_decompose_a, df_decompose_b, df_decompose_c, df_decompose_d)

df_decompose$type <- factor(df_decompose$type , levels = c("Climate change", "Pop. growth", "Aging", "Historical"))

df_decompose_r <- df_decompose

for(i in 1:5){
  ggplot(df_decompose_r)+
    theme_classic()+
    geom_col(aes(x=scenario, y=value/1e6, fill=type))+
    facet_wrap_paginate(vars(reg), scales = "free", ncol = 6, nrow = 6, page = i)+
    ylab(expression(paste("Growth in",phantom(x),phantom(x), PD, ", (million)")))+
    xlab("Scenario")+
    scale_fill_manual(name="Driver", values = c(ggsci::pal_npg("nrc", alpha =1)(4)[2:4], ggsci::pal_npg("nrc", alpha =1)(9)[1]))+
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0("decomposition_country_", i, "_tmax.pdf"), scale=2, width = 8, height = 6)
}

###############

# Table 1


regions <- st_as_sf(rworldmap::countriesLow)
regions$REGION_n <- as.numeric(regions$REGION)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="REGION_n")
regions <- dplyr::select(regions, REGION, REGION_n)
regions$geometry<-NULL


#region // >69% // <69% //// SHARE of tot pop /// tot POP /// pop_wgtd CDD //// pop_wgtd TMAS

df <- data.frame(totpop=values(r1_tot), cdd = values(multigcm_hist), tasmax = values(multigcm_tasmax_hist), hotdays=values(multigcm_hotday_hist), reg = values(regions_r), share=values(r1_g[[3]]) / values(r1_tot), pdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hist), tmaxdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_tasmax_hist),  hotdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hotday_hist),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "Historical"

df1 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_tot), cdd = values(multigcm_hist), tasmax = values(multigcm_tasmax_hist), hotdays=values(multigcm_hotday_hist), reg = 8, share=values(r1_g[[3]]) / values(r1_tot), pdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hist), tmaxdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_tasmax_hist),  hotdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hotday_hist),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "Historical"

df1_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

###

df <- data.frame(totpop=values(r1_ssp2_2050_tot), cdd = values(multigcm_245), tasmax = values(multigcm_tasmax_245), hotdays=values(multigcm_hotday_245[[1:13]]), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot), pdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_245), tmaxdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_tasmax_245),  hotdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_hotday_245[[1:13]]),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP2, 2050"

df2 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

df <- data.frame(totpop=values(r1_ssp2_2050_tot), cdd = values(multigcm_245), tasmax = values(multigcm_tasmax_245), hotdays=values(multigcm_hotday_245[[1:13]]), reg = 8, share=values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot), pdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_245), tmaxdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_tasmax_245),  hotdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_hotday_245[[1:13]]),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP2, 2050"

df2_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

###

df <- data.frame(totpop=values(r1_ssp5_2050_tot), cdd = values(multigcm_585), tasmax = values(multigcm_tasmax_585), hotdays=values(multigcm_hotday_585), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot), pdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_585), tmaxdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_tasmax_585),  hotdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_hotday_585),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP5, 2050"

df3 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_ssp5_2050_tot), cdd = values(multigcm_585), tasmax = values(multigcm_tasmax_585), hotdays=values(multigcm_hotday_585), reg = 8, share=values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot), pdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_585), tmaxdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_tasmax_585),  hotdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_hotday_585),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP5, 2050"

df3_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()


###


df <- data.frame(totpop=values(r1_ssp1_2050_tot), cdd = values(multigcm_126), tasmax = values(multigcm_tasmax_126), hotdays=values(multigcm_hotday_126), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot), pdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_126), tmaxdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_tasmax_126),  hotdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_hotday_126),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP1, 2050"

df4 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_ssp1_2050_tot), cdd = values(multigcm_126), tasmax = values(multigcm_tasmax_126), hotdays=values(multigcm_hotday_126), reg = 8, share=values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot), pdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_126), tmaxdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_tasmax_126),  hotdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_hotday_126),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP1, 2050"

df4_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

###


df <- data.frame(totpop=values(r1_ssp3_2050_tot), cdd = values(multigcm_370), tasmax = values(multigcm_tasmax_370), hotdays=values(multigcm_hotday_370), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot), pdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_370), tmaxdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_tasmax_370),  hotdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_hotday_370),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP3, 2050"

df5 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_ssp3_2050_tot), cdd = values(multigcm_370), tasmax = values(multigcm_tasmax_370), hotdays=values(multigcm_hotday_370), reg = 8, share=values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot), pdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_370), tmaxdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_tasmax_370),  hotdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_hotday_370),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP3, 2050"

df5_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

###

df <- bind_rows(df1, df1_glo, df2, df2_glo, df3, df3_glo, df4, df4_glo, df5, df5_glo)

rm(df1, df1_glo, df2, df2_glo, df3, df3_glo, df4, df4_glo, df5, df5_glo)
gc()

########

emptycol <- function(x) " "

df <- reshape2::melt(df, c(1, 4:22))
df <- reshape2::melt(df, c(1, 20:22))
colnames(df) <- c("reg", "scenario", "variable", "value", "clim_var", "clim_value")

df$reg <- as.factor(df$reg)

df$value[df$variable=="totpops"] <- df$value[df$variable=="totpops"] / 1e6
df$value[df$variable=="share"] <- df$value[df$variable=="share"] * 100
df$clim_value[grepl("tmaxdd", df$clim_var)] <- df$clim_value[grepl("tmaxdd", df$clim_var)] / 1e6
df$clim_value[grepl("hotdd", df$clim_var)] <- df$clim_value[grepl("hotdd", df$clim_var)] / 1e6
df$clim_value[grepl("pdd", df$clim_var)] <- df$clim_value[grepl("pdd", df$clim_var)] / 1e9

df <- df %>% 
  mutate_if(is.numeric, round, 1)

df$variable <- as.character(df$variable)
df$clim_var <- as.character(df$clim_var)

df$variable[df$variable=="totpops"] <- "Pop. total (million)"
df$variable[df$variable=="share"] <- "% of pop. >69 y.o."

df$clim_var[df$clim_var=="q5_cdd"] <- "Pop. weighted CDDs, 5th PCTL"
df$clim_var[df$clim_var=="q50_cdd"] <- "Pop. weighted CDDs, 50th PCTL"
df$clim_var[df$clim_var=="q95_cdd"] <- "Pop. weighted CDDs, 95th PCTL"

df$clim_var[df$clim_var=="q5_hotdays"] <- "Pop. weighted # hot days, 5th PCTL"
df$clim_var[df$clim_var=="q50_hotdays"] <- "Pop. weighted # hot days, 50th PCTL"
df$clim_var[df$clim_var=="q95_hotdays"] <- "Pop. weighted # hot days, 95th PCTL"

df$clim_var[df$clim_var=="q5_tasmax"] <- "Pop. weighted TMAX95th, 5th PCTL"
df$clim_var[df$clim_var=="q50_tasmax"] <- "Pop. weighted TMAX95th, 50th PCTL"
df$clim_var[df$clim_var=="q95_tasmax"] <- "Pop. weighted TMAX95th, 95th PCTL"

df$clim_var[df$clim_var=="q5_pdd"] <- "Pop. weighted PDDs, 5th PCTL"
df$clim_var[df$clim_var=="q50_pdd"] <- "Pop. weighted PDDs, 50th PCTL"
df$clim_var[df$clim_var=="q95_pdd"] <- "Pop. weighted PDDs, 95th PCTL"

df$clim_var[df$clim_var=="q5_tmaxdd"] <- "Pop. weighted PD95th, 5th PCTL"
df$clim_var[df$clim_var=="q50_tmaxdd"] <- "Pop. weighted PD95th, 50th PCTL"
df$clim_var[df$clim_var=="q95_tmaxdd"] <- "Pop. weighted PD95th, 95th PCTL"

df$clim_var[df$clim_var=="q5_hotdd"] <- "Pop. weighted PHD, 5th PCTL"
df$clim_var[df$clim_var=="q50_hotdd"] <- "Pop. weighted PHD, 50th PCTL"
df$clim_var[df$clim_var=="q95_hotdd"] <- "Pop. weighted PHD, 95th PCTL"


df$reg <- ifelse(!is.na(regions$REGION[match(df$reg, as.numeric(regions$REGION))]), as.character(regions$REGION[match(df$reg, as.numeric(regions$REGION))]), " Global total")
df$reg <- as.factor(df$reg)

df$cli_var_type <- ifelse(grepl("Pop. weighted CDDs,", df$clim_var), "CDDs", ifelse(grepl("Pop. weighted TMAX95th,", df$clim_var), "TMAX95th", ifelse(grepl("Pop. weighted PDDs,", df$clim_var), "PDDs", ifelse(grepl("Pop. weighted PHD,", df$clim_var), "PHDs", ifelse(grepl("Pop. weighted PD95th,", df$clim_var), "PD95th", "Hot days")))))

df_dci <- df %>% group_by(reg, scenario, variable, cli_var_type) %>% mutate(clim_value = paste0(as.character(round(clim_value[2], 1)), " (", as.character(round(clim_value[1], 1)), " - ", as.character(round(clim_value[3], 1)), ")"))

df_dci$cli_var_type <- NULL

df_dci$value <- as.character(df_dci$value)

nrr = nrow(df_dci)

df_dci[c(1:nrr)+nrr,] <- df_dci[1:nrr,]
df_dci[c(1:nrr)+nrr,c(3:4)] <- df_dci[c(1:nrr)+nrr,c(5:6)]
df_dci[,c(5:6)] <- NULL

df_dci <- filter(df_dci, variable!="Pop. weighted CDDs, 5th PCTL" & variable!="Pop. weighted CDDs, 95th PCTL", variable!="Pop. weighted TMAX95th, 5th PCTL" & variable!="Pop. weighted TMAX95th, 95th PCTL" & variable!="Pop. weighted PDDs, 5th PCTL" & variable!="Pop. weighted PDDs, 95th PCTL", variable!="Pop. weighted PD95th, 5th PCTL" & variable!="Pop. weighted PD95th, 95th PCTL" , variable!="Pop. weighted # hot days, 5th PCTL" & variable!="Pop. weighted # hot days, 95th PCTL", variable!="Pop. weighted PHD, 5th PCTL" & variable!="Pop. weighted PHD, 95th PCTL")


df_dci <- distinct(df_dci)
df_dci_hist <- spread(df_dci, key = variable, value = value)
df_dci_hist <- arrange(df_dci_hist, scenario, reg)

df_dci_hist <- df_dci_hist[,c(1, 2, 4, 3, 6, 8, 10, 7, 5, 9)]  

colnames(df_dci_hist)[1] <- "Region"
colnames(df_dci_hist)[2] <- "Scenario"

colnames(df_dci_hist) <- gsub(", 50th PCTL", "", colnames(df_dci_hist))

colnames(df_dci_hist)[6] <- paste0(colnames(df_dci_hist)[6], " (billion)")
colnames(df_dci_hist)[c(8, 10)] <- paste0(colnames(df_dci_hist)[c(8, 10)], " (million)")

tt <- xtable(df_dci_hist)
align(tt) <- rep("c", 11)
digits(tt) <- 1

write_rds(df_dci_hist, "df_dci_region.Rds")

sink("tab1_region.tex")
print(tt, include.rownames=FALSE, booktabs = TRUE, floating.environment = "sidewaystable")
sink()

openxlsx::write.xlsx(dplyr::arrange(df_dci_hist, Region), "age_byreg.xlsx")

####
# summary table by country


regions <- st_as_sf(rworldmap::countriesLow)
regions$REGION_n <- as.numeric(regions$SOV_A3)
regions_r <- fasterize::fasterize(regions, r1_g[[1]], field="REGION_n")
regions <- dplyr::select(regions, SOV_A3, REGION_n)
regions$geometry<-NULL

df <- data.frame(totpop=values(r1_tot), cdd = values(multigcm_hist), tasmax = values(multigcm_tasmax_hist), hotdays=values(multigcm_hotday_hist), reg = values(regions_r), share=values(r1_g[[3]]) / values(r1_tot), pdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hist), tmaxdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_tasmax_hist),  hotdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hotday_hist),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "Historical"

df1 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_tot), cdd = values(multigcm_hist), tasmax = values(multigcm_tasmax_hist), hotdays=values(multigcm_hotday_hist), reg = 8, share=values(r1_g[[3]]) / values(r1_tot), pdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hist), tmaxdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_tasmax_hist),  hotdd=values(r1_tot) * (values(r1_g[[3]]) / values(r1_tot)) * values(multigcm_hotday_hist),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "Historical"

df1_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

###

df <- data.frame(totpop=values(r1_ssp2_2050_tot), cdd = values(multigcm_245), tasmax = values(multigcm_tasmax_245), hotdays=values(multigcm_hotday_245[[1:13]]), reg = values(regions_r), share=values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot), pdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_245), tmaxdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_tasmax_245),  hotdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_hotday_245[[1:13]]),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP2, 2050"

df2 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

df <- data.frame(totpop=values(r1_ssp2_2050_tot), cdd = values(multigcm_245), tasmax = values(multigcm_tasmax_245), hotdays=values(multigcm_hotday_245[[1:13]]), reg = 8, share=values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot), pdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_245), tmaxdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_tasmax_245),  hotdd=values(r1_ssp2_2050_tot) * (values(r1_ssp2_2050[[3]]) / values(r1_ssp2_2050_tot)) * values(multigcm_hotday_245[[1:13]]),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP2, 2050"

df2_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

###

df <- data.frame(totpop=values(r1_ssp5_2050_tot), cdd = values(multigcm_585), tasmax = values(multigcm_tasmax_585), hotdays=values(multigcm_hotday_585), reg = values(regions_r), share=values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot), pdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_585), tmaxdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_tasmax_585),  hotdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_hotday_585),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP5, 2050"

df3 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_ssp5_2050_tot), cdd = values(multigcm_585), tasmax = values(multigcm_tasmax_585), hotdays=values(multigcm_hotday_585), reg = 8, share=values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot), pdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_585), tmaxdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_tasmax_585),  hotdd=values(r1_ssp5_2050_tot) * (values(r1_ssp5_2050[[3]]) / values(r1_ssp5_2050_tot)) * values(multigcm_hotday_585),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP5, 2050"

df3_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()


###


df <- data.frame(totpop=values(r1_ssp1_2050_tot), cdd = values(multigcm_126), tasmax = values(multigcm_tasmax_126), hotdays=values(multigcm_hotday_126), reg = values(regions_r), share=values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot), pdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_126), tmaxdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_tasmax_126),  hotdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_hotday_126),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP1, 2050"

df4 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_ssp1_2050_tot), cdd = values(multigcm_126), tasmax = values(multigcm_tasmax_126), hotdays=values(multigcm_hotday_126), reg = 8, share=values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot), pdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_126), tmaxdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_tasmax_126),  hotdd=values(r1_ssp1_2050_tot) * (values(r1_ssp1_2050[[3]]) / values(r1_ssp1_2050_tot)) * values(multigcm_hotday_126),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP1, 2050"

df4_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

###


df <- data.frame(totpop=values(r1_ssp3_2050_tot), cdd = values(multigcm_370), tasmax = values(multigcm_tasmax_370), hotdays=values(multigcm_hotday_370), reg = values(regions_r), share=values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot), pdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_370), tmaxdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_tasmax_370),  hotdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_hotday_370),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP3, 2050"

df5 <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

#

df <- data.frame(totpop=values(r1_ssp3_2050_tot), cdd = values(multigcm_370), tasmax = values(multigcm_tasmax_370), hotdays=values(multigcm_hotday_370), reg = 8, share=values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot), pdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_370), tmaxdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_tasmax_370),  hotdd=values(r1_ssp3_2050_tot) * (values(r1_ssp3_2050[[3]]) / values(r1_ssp3_2050_tot)) * values(multigcm_hotday_370),xy=T)

df <- na.omit(df)

df_a <- df %>% group_by(reg) %>% dplyr::summarise(totpops=sum(totpop, na.rm = T), share = weighted.mean(share, totpop, na.rm=T))

df_b <- df %>% group_by(reg) %>% dplyr::select(starts_with("cdd"), starts_with("tasmax"), starts_with("hotdays"), totpop, share) %>% dplyr::summarise_all(funs(weighted.mean(., totpop*share, na.rm=T))) %>% ungroup()

df_b_1 <-  dplyr::select(df_b, reg, starts_with("cdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_cdd = as.numeric(quantile(value, 0.25)), q50_cdd =  as.numeric(quantile(value, 0.5)), q95_cdd =  as.numeric(quantile(value, 0.75)))
df_b_2 <-  dplyr::select(df_b, reg, starts_with("tasmax")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tasmax = as.numeric(quantile(value, 0.25)), q50_tasmax =  as.numeric(quantile(value, 0.5)), q95_tasmax =  as.numeric(quantile(value, 0.75)))
df_b_3 <-  dplyr::select(df_b, reg, starts_with("hotdays")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdays = as.numeric(quantile(value, 0.25)), q50_hotdays =  as.numeric(quantile(value, 0.5)), q95_hotdays =  as.numeric(quantile(value, 0.75)))

df_b_1$reg <- NULL
df_b_2$reg <- NULL
df_b_3$reg <- NULL

df_b <- bind_cols(df_b_1, df_b_2, df_b_3)

df_c <- df %>% group_by(reg) %>% dplyr::select(starts_with("pdd"), starts_with("tmaxdd"), starts_with("hotdd")) %>% dplyr::summarise_all(funs(sum(., na.rm=T))) %>% ungroup()

df_c_1 <-  dplyr::select(df_c, reg, starts_with("pdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_pdd = as.numeric(quantile(value, 0.25)), q50_pdd =  as.numeric(quantile(value, 0.5)), q95_pdd =  as.numeric(quantile(value, 0.75)))
df_c_2 <-  dplyr::select(df_c, reg, starts_with("tmaxdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_tmaxdd = as.numeric(quantile(value, 0.25)), q50_tmaxdd =  as.numeric(quantile(value, 0.5)), q95_tmaxdd =  as.numeric(quantile(value, 0.75)))
df_c_3 <-  dplyr::select(df_c, reg, starts_with("hotdd")) %>% reshape2::melt(1) %>% group_by(reg) %>% dplyr::summarise(q5_hotdd = as.numeric(quantile(value, 0.25)), q50_hotdd =  as.numeric(quantile(value, 0.5)), q95_hotdd =  as.numeric(quantile(value, 0.75)))

df_c_1$reg <- NULL
df_c_2$reg <- NULL
df_c_3$reg <- NULL

df_c <- bind_cols(df_c_1, df_c_2, df_c_3)

df <- bind_cols(df_a, df_b, df_c)

df$scenario <- "SSP3, 2050"

df5_glo <- df

rm(df, df_a, df_b, df_c, df_a_1, df_a_2, df_a_3, df_b_1, df_b_2, df_b_3, df_c_1, df_c_2, df_c_3)
gc()

###

df <- bind_rows(df1, df1_glo, df2, df2_glo, df3, df3_glo, df4, df4_glo, df5, df5_glo)

rm(df1, df1_glo, df2, df2_glo, df3, df3_glo, df4, df4_glo, df5, df5_glo)
gc()

########

emptycol <- function(x) " "

df <- reshape2::melt(df, c(1, 4:22))
df <- reshape2::melt(df, c(1, 20:22))
colnames(df) <- c("reg", "scenario", "variable", "value", "clim_var", "clim_value")

df$reg <- as.factor(df$reg)

df$value[df$variable=="totpops"] <- df$value[df$variable=="totpops"] / 1e6
df$value[df$variable=="share"] <- df$value[df$variable=="share"] * 100
df$clim_value[grepl("tmaxdd", df$clim_var)] <- df$clim_value[grepl("tmaxdd", df$clim_var)] / 1e6
df$clim_value[grepl("hotdd", df$clim_var)] <- df$clim_value[grepl("hotdd", df$clim_var)] / 1e6
df$clim_value[grepl("pdd", df$clim_var)] <- df$clim_value[grepl("pdd", df$clim_var)] / 1e9

df <- df %>% 
  mutate_if(is.numeric, round, 1)

df$variable <- as.character(df$variable)
df$clim_var <- as.character(df$clim_var)

df$variable[df$variable=="totpops"] <- "Pop. total (million)"
df$variable[df$variable=="share"] <- "% of pop. >69 y.o."

df$clim_var[df$clim_var=="q5_cdd"] <- "Pop. weighted CDDs, 5th PCTL"
df$clim_var[df$clim_var=="q50_cdd"] <- "Pop. weighted CDDs, 50th PCTL"
df$clim_var[df$clim_var=="q95_cdd"] <- "Pop. weighted CDDs, 95th PCTL"

df$clim_var[df$clim_var=="q5_hotdays"] <- "Pop. weighted # hot days, 5th PCTL"
df$clim_var[df$clim_var=="q50_hotdays"] <- "Pop. weighted # hot days, 50th PCTL"
df$clim_var[df$clim_var=="q95_hotdays"] <- "Pop. weighted # hot days, 95th PCTL"

df$clim_var[df$clim_var=="q5_tasmax"] <- "Pop. weighted TMAX95th, 5th PCTL"
df$clim_var[df$clim_var=="q50_tasmax"] <- "Pop. weighted TMAX95th, 50th PCTL"
df$clim_var[df$clim_var=="q95_tasmax"] <- "Pop. weighted TMAX95th, 95th PCTL"

df$clim_var[df$clim_var=="q5_pdd"] <- "Pop. weighted PDDs, 5th PCTL"
df$clim_var[df$clim_var=="q50_pdd"] <- "Pop. weighted PDDs, 50th PCTL"
df$clim_var[df$clim_var=="q95_pdd"] <- "Pop. weighted PDDs, 95th PCTL"

df$clim_var[df$clim_var=="q5_tmaxdd"] <- "Pop. weighted PD95th, 5th PCTL"
df$clim_var[df$clim_var=="q50_tmaxdd"] <- "Pop. weighted PD95th, 50th PCTL"
df$clim_var[df$clim_var=="q95_tmaxdd"] <- "Pop. weighted PD95th, 95th PCTL"

df$clim_var[df$clim_var=="q5_hotdd"] <- "Pop. weighted PHD, 5th PCTL"
df$clim_var[df$clim_var=="q50_hotdd"] <- "Pop. weighted PHD, 50th PCTL"
df$clim_var[df$clim_var=="q95_hotdd"] <- "Pop. weighted PHD, 95th PCTL"

df$reg <- ifelse(!is.na(regions$SOV_A3[match(df$reg, regions$REGION_n)]), as.character(regions$SOV_A3[match(df$reg, regions$REGION_n)]), " Global total")
df$reg <- as.factor(df$reg)

df$cli_var_type <- ifelse(grepl("Pop. weighted CDDs,", df$clim_var), "CDDs", ifelse(grepl("Pop. weighted TMAX95th,", df$clim_var), "TMAX95th", ifelse(grepl("Pop. weighted PDDs,", df$clim_var), "PDDs", ifelse(grepl("Pop. weighted PHD,", df$clim_var), "PHDs", ifelse(grepl("Pop. weighted PD95th,", df$clim_var), "PD95th", "Hot days")))))

df_dci <- df %>% group_by(reg, scenario, variable, cli_var_type) %>% mutate(clim_value = paste0(as.character(round(clim_value[2], 1)), " (", as.character(round(clim_value[1], 1)), " - ", as.character(round(clim_value[3], 1)), ")"))

df_dci$cli_var_type <- NULL

df_dci$value <- as.character(df_dci$value)

nrr = nrow(df_dci)

df_dci[c(1:nrr)+nrr,] <- df_dci[1:nrr,]
df_dci[c(1:nrr)+nrr,c(3:4)] <- df_dci[c(1:nrr)+nrr,c(5:6)]
df_dci[,c(5:6)] <- NULL

df_dci <- filter(df_dci, variable!="Pop. weighted CDDs, 5th PCTL" & variable!="Pop. weighted CDDs, 95th PCTL", variable!="Pop. weighted TMAX95th, 5th PCTL" & variable!="Pop. weighted TMAX95th, 95th PCTL" & variable!="Pop. weighted PDDs, 5th PCTL" & variable!="Pop. weighted PDDs, 95th PCTL", variable!="Pop. weighted PD95th, 5th PCTL" & variable!="Pop. weighted PD95th, 95th PCTL" , variable!="Pop. weighted # hot days, 5th PCTL" & variable!="Pop. weighted # hot days, 95th PCTL", variable!="Pop. weighted PHD, 5th PCTL" & variable!="Pop. weighted PHD, 95th PCTL")


df_dci <- distinct(df_dci)
df_dci_hist <- spread(df_dci, key = variable, value = value)
df_dci_hist <- arrange(df_dci_hist, scenario, reg)

df_dci_hist <- df_dci_hist[,c(1, 2, 4, 3, 6, 8, 10, 7, 5, 9)]  

colnames(df_dci_hist)[1] <- "Country"
colnames(df_dci_hist)[2] <- "Scenario"

colnames(df_dci_hist) <- gsub(", 50th PCTL", "", colnames(df_dci_hist))

colnames(df_dci_hist)[6] <- paste0(colnames(df_dci_hist)[6], " (billion)")
colnames(df_dci_hist)[c(8, 10)] <- paste0(colnames(df_dci_hist)[c(8, 10)], " (million)")

tt <- xtable(df_dci_hist)
align(tt) <- rep("c", 11)
digits(tt) <- 1

write_rds(df_dci_hist, "df_dci_country.Rds")

tt <- xtable(df_dci_hist)
align(tt) <- rep("c", 11)
digits(tt) <- 1

sink("tab1_country.tex")
print(tt, tabular.environment="longtable", floating=FALSE, booktabs = TRUE)
sink()
