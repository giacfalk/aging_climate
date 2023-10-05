devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)
devtools::install_github("https://github.com/UrbanInstitute/urbnmapr")
library(urbnmapr)

# 

r1_g <- stack(paste0(getwd(), "/r1_2020.tif"))
r1_g_245 <- stack(paste0(getwd(), "/r1_g_245_2050_final.tif"))
r1_g_585 <- stack(paste0(getwd(), "/r1_g_585_2050_final.tif"))
r1_g_126 <- stack(paste0(getwd(), "/r1_g_126_2050_final.tif"))
r1_g_370 <- stack(paste0(getwd(), "/r1_g_370_2050_final.tif"))

# validate for US

r1_ssp5_2050 <- r1_g_585

r1_ssp2_2050 <- r1_g_245 

r1_ssp3_2050 <- r1_g_370

r1_ssp1_2050 <- r1_g_126 

#

load(paste0(getwd(), "/proj_res.RData"))

counties_sf <- get_urbn_map("counties", sf = TRUE)

proj_res$STATE_COUNTY <- paste(proj_res$STATE, proj_res$COUNTY, sep="_")

proj_res2 <- filter(proj_res, grepl(" County", COUNTY))

proj_res2 <- filter(proj_res2, YEAR==2050)

# aggregations for comparability

proj_res2$AGEGR <- ifelse(proj_res2$AGEGR=="0--19", "age_gr_0_19", ifelse(proj_res2$AGEGR=="20--29" | proj_res2$AGEGR=="30--54" | proj_res2$AGEGR=="55--69" , "age_gr_20_69", "age_gr_69+"))

proj_res2 <- group_by(proj_res2, COUNTY, STATE, YEAR, SCEN, AGEGR) %>% dplyr::summarise(TPOP=sum(TPOP), PROP=sum(PROP))

#

proj_res2_ssp5 <- proj_res2 %>% filter(SCEN=="ssp5")

proj_res2_ssp5 <- proj_res2_ssp5 %>%
  pivot_wider(names_from = AGEGR, values_from = c(TPOP, PROP))

proj_res2_ssp2 <- proj_res2 %>% filter(SCEN=="ssp2")

proj_res2_ssp2 <- proj_res2_ssp2 %>%
  pivot_wider(names_from = AGEGR, values_from = c(TPOP, PROP))

#

proj_res3_ssp5 <- merge(proj_res2_ssp5, counties_sf, by.x=c("COUNTY", "STATE"),  by.y=c("county_name", "state_name"))
proj_res3_ssp5 <- st_as_sf(proj_res3_ssp5)
proj_res3_ssp5 <- st_transform(proj_res3_ssp5, 4326)

proj_res3_ssp2 <- merge(proj_res2_ssp2, counties_sf, by.x=c("COUNTY", "STATE"),  by.y=c("county_name", "state_name"))
proj_res3_ssp2 <- st_as_sf(proj_res3_ssp2)
proj_res3_ssp2 <- st_transform(proj_res3_ssp2, 4326)

#

proj_res2_ssp3 <- proj_res2 %>% filter(SCEN=="ssp3")

proj_res2_ssp3 <- proj_res2_ssp3 %>%
  pivot_wider(names_from = AGEGR, values_from = c(TPOP, PROP))

proj_res2_ssp1 <- proj_res2 %>% filter(SCEN=="ssp1")

proj_res2_ssp1 <- proj_res2_ssp1 %>%
  pivot_wider(names_from = AGEGR, values_from = c(TPOP, PROP))

#

proj_res3_ssp3 <- merge(proj_res2_ssp3, counties_sf, by.x=c("COUNTY", "STATE"),  by.y=c("county_name", "state_name"))
proj_res3_ssp3 <- st_as_sf(proj_res3_ssp3)
proj_res3_ssp3 <- st_transform(proj_res3_ssp3, 4326)

proj_res3_ssp1 <- merge(proj_res2_ssp1, counties_sf, by.x=c("COUNTY", "STATE"),  by.y=c("county_name", "state_name"))
proj_res3_ssp1 <- st_as_sf(proj_res3_ssp1)
proj_res3_ssp1 <- st_transform(proj_res3_ssp1, 4326)

#

r1_ssp5_2050_e <- exact_extract(r1_ssp5_2050, proj_res3_ssp5, "sum", max_cells_in_memory=1e09)

r1_ssp2_2050_e <- exact_extract(r1_ssp2_2050,proj_res3_ssp2, "sum", max_cells_in_memory=1e09)

# merge back to data
colnames(r1_ssp5_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
proj_res3_ssp5 <- bind_cols(proj_res3_ssp5, r1_ssp5_2050_e)

colnames(r1_ssp2_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
proj_res3_ssp2 <- bind_cols(proj_res3_ssp2, r1_ssp2_2050_e)

r1_ssp3_2050_e <- exact_extract(r1_ssp3_2050, proj_res3_ssp3, "sum", max_cells_in_memory=1e09)

r1_ssp1_2050_e <- exact_extract(r1_ssp1_2050,proj_res3_ssp1, "sum", max_cells_in_memory=1e09)

# merge back to data
colnames(r1_ssp3_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
proj_res3_ssp3 <- bind_cols(proj_res3_ssp3, r1_ssp3_2050_e)

colnames(r1_ssp1_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
proj_res3_ssp1 <- bind_cols(proj_res3_ssp1, r1_ssp1_2050_e)



# calculate MAPEs and MSEs

mape <- function(real,modeled){
  
  round(mean(abs((real - modeled)/real), na.rm=T), 2)
  
}


mse <- function(real,modeled){
  
  round(mean((real - modeled)^2, na.rm=T), 2)
  
}

#########

# a<- ggplot()+
#   geom_point(aes(y=proj_res3_ssp2$TPOP_age_gr_0_19, x=proj_res3_ssp2$age_gr_0_19), alpha=0.1) + 
#   scale_x_log10(limits= c(1000, 100000))+
#   scale_y_log10(limits= c(1000, 100000))+
#   geom_abline()+
#   xlab("This study (SSP2)")+
#   ylab("O' Neil et al. (SSP2)")+
#   annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp2$TPOP_age_gr_0_19, proj_res3_ssp2$age_gr_0_19), 2)))+
#   annotate("text", x=2.2e+03, y=0.8e+05, label= paste0("MAPE = ", mape(proj_res3_ssp2$TPOP_age_gr_0_19, proj_res3_ssp2$age_gr_0_19)*100, " %"))+
#   ggtitle("Age group 0-19")
# 
# b<- ggplot()+
#   geom_point(aes(y=proj_res3_ssp2$TPOP_age_gr_20-60, x=proj_res3_ssp2$age_gr_20-60), alpha=0.1) + 
#   scale_x_log10(limits= c(1000, 100000))+
#   scale_y_log10(limits= c(1000, 100000))+
#   geom_abline()+
#   xlab("This study (SSP2)")+
#   ylab("O' Neil et al. (SSP2)")+
#   annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp2$TPOP_age_gr_20-60, proj_res3_ssp2$age_gr_20-60), 2)))+
#   annotate("text", x=2.2e+03, y=0.8e+05, label= paste0("MAPE = ", mape(proj_res3_ssp2$TPOP_age_gr_20-60, proj_res3_ssp2$age_gr_20-60)*100, " %"))+
#   ggtitle("Age group 20-69")

proj_res3_ssp1 <- filter(proj_res3_ssp1, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
proj_res3_ssp1 <- filter(proj_res3_ssp1, `TPOP_age_gr_69+`>0 & !is.na(`TPOP_age_gr_69+`))
proj_res3_ssp2 <- filter(proj_res3_ssp2, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
proj_res3_ssp2 <- filter(proj_res3_ssp2, `TPOP_age_gr_69+`>0 & !is.na(`TPOP_age_gr_69+`))
proj_res3_ssp3 <- filter(proj_res3_ssp3, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
proj_res3_ssp3 <- filter(proj_res3_ssp3, `TPOP_age_gr_69+`>0 & !is.na(`TPOP_age_gr_69+`))
proj_res3_ssp5 <- filter(proj_res3_ssp5, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
proj_res3_ssp5 <- filter(proj_res3_ssp5, `TPOP_age_gr_69+`>0 & !is.na(`TPOP_age_gr_69+`))

a<- ggplot()+
  geom_point(aes(y=proj_res3_ssp2$`TPOP_age_gr_69+`, x=proj_res3_ssp2$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2020, SSP2), tot. = ", round(sum(proj_res3_ssp2$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("O' Neil et al. (2020, SSP2), tot. = ", round(sum(proj_res3_ssp2$`TPOP_age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp2$`TPOP_age_gr_69+`, proj_res3_ssp2$`age_gr_69+`), 2)))+
  annotate("text", x=1.8e+03, y=0.7e+05, label= paste0("MAPE = ", mape(proj_res3_ssp2$`TPOP_age_gr_69+`, proj_res3_ssp2$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

# d<- ggplot()+
#   geom_point(aes(y=proj_res3_ssp5$TPOP_age_gr_0_19, x=proj_res3_ssp5$age_gr_0_19), alpha=0.1) + 
#   scale_x_log10(limits= c(1000, 100000))+
#   scale_y_log10(limits= c(1000, 100000))+
#   geom_abline()+
#   xlab("This study (SSP5)")+
#   ylab("O' Neil et al. (SSP5)")+
#   annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp5$TPOP_age_gr_0_19, proj_res3_ssp5$age_gr_0_19), 2)))+
#   annotate("text", x=2.2e+03, y=0.8e+05, label= paste0("MAPE = ", mape(proj_res3_ssp5$TPOP_age_gr_0_19, proj_res3_ssp5$age_gr_0_19)*100, " %"))+
#   ggtitle("Age group 0-19")
# 
# e<- ggplot()+
#   geom_point(aes(y=proj_res3_ssp5$TPOP_age_gr_20-60, x=proj_res3_ssp5$age_gr_20-60), alpha=0.1) + 
#   scale_x_log10(limits= c(1000, 100000))+
#   scale_y_log10(limits= c(1000, 100000))+
#   geom_abline()+
#   xlab("This study (SSP5)")+
#   ylab("O' Neil et al. (SSP5)")+
#   annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp5$TPOP_age_gr_20-60, proj_res3_ssp5$age_gr_20-60), 2)))+
#   annotate("text", x=2.2e+03, y=0.8e+05, label= paste0("MAPE = ", mape(proj_res3_ssp5$TPOP_age_gr_20-60, proj_res3_ssp5$age_gr_20-60)*100, " %"))+
#   ggtitle("Age group 20-69")

b<- ggplot()+
  geom_point(aes(y=proj_res3_ssp5$`TPOP_age_gr_69+`, x=proj_res3_ssp5$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2020, SSP5), tot. = ", round(sum(proj_res3_ssp5$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("O' Neil et al. (2020, SSP5), tot. = ", round(sum(proj_res3_ssp5$`TPOP_age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp5$`TPOP_age_gr_69+`, proj_res3_ssp5$`age_gr_69+`), 2)))+
  annotate("text", x=1.8e+03, y=0.7e+05, label= paste0("MAPE = ", mape(proj_res3_ssp5$`TPOP_age_gr_69+`, proj_res3_ssp5$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

c<- ggplot()+
  geom_point(aes(y=proj_res3_ssp1$`TPOP_age_gr_69+`, x=proj_res3_ssp1$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2020, SSP1), tot. = ", round(sum(proj_res3_ssp1$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("O' Neil et al. (2020, SSP1), tot. = ", round(sum(proj_res3_ssp1$`TPOP_age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp1$`TPOP_age_gr_69+`, proj_res3_ssp1$`age_gr_69+`), 2)))+
  annotate("text", x=1.8e+03, y=0.7e+05, label= paste0("MAPE = ", mape(proj_res3_ssp1$`TPOP_age_gr_69+`, proj_res3_ssp1$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

d<- ggplot()+
  geom_point(aes(y=proj_res3_ssp3$`TPOP_age_gr_69+`, x=proj_res3_ssp3$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2020, SSP3), tot. = ", round(sum(proj_res3_ssp3$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("O' Neil et al. (2020, SSP3), tot. = ", round(sum(proj_res3_ssp3$`TPOP_age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.5e+03, y=1e+05, label= paste0("r = ", round(cor(proj_res3_ssp3$`TPOP_age_gr_69+`, proj_res3_ssp3$`age_gr_69+`), 2)))+
  annotate("text", x=1.8e+03, y=0.7e+05, label= paste0("MAPE = ", mape(proj_res3_ssp3$`TPOP_age_gr_69+`, proj_res3_ssp3$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

#

panel <- a + b + c + d + plot_layout(ncol=2) + plot_annotation(
  title = 'United States counties',
  theme = theme(plot.title = element_text(size = 15)))

ggsave("us_comparison_2050.png", panel, height = 8, width = 11)

###

c <- ggplot(proj_res3_ssp2)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp2$`TPOP_age_gr_69+`, proj_res3_ssp2$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP2")

proj_res3_ssp2_s <- proj_res3_ssp2 %>% group_by(STATE) %>% dplyr::summarise(`TPOP_age_gr_69+`=sum(`TPOP_age_gr_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

d <- ggplot(proj_res3_ssp2_s)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp2_s$`TPOP_age_gr_69+`, proj_res3_ssp2_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP2")

sum(proj_res3_ssp2$`TPOP_age_gr_69+`)/1e6
sum(proj_res3_ssp2$`age_gr_69+`)/1e6

g <- ggplot(proj_res3_ssp5)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp5$`TPOP_age_gr_69+`, proj_res3_ssp5$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP5")

proj_res3_ssp5_s <- proj_res3_ssp5 %>% group_by(STATE) %>% dplyr::summarise(`TPOP_age_gr_69+`=sum(`TPOP_age_gr_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

h <- ggplot(proj_res3_ssp5_s)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp5_s$`TPOP_age_gr_69+`, proj_res3_ssp5_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP5")

sum(proj_res3_ssp5$`TPOP_age_gr_69+`)/1e6
sum(proj_res3_ssp5$`age_gr_69+`)/1e6

a <- ggplot(proj_res3_ssp1)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp1$`TPOP_age_gr_69+`, proj_res3_ssp1$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP1")

proj_res3_ssp1_s <- proj_res3_ssp1 %>% group_by(STATE) %>% dplyr::summarise(`TPOP_age_gr_69+`=sum(`TPOP_age_gr_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

b <- ggplot(proj_res3_ssp1_s)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp1_s$`TPOP_age_gr_69+`, proj_res3_ssp1_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP1")


sum(proj_res3_ssp1$`TPOP_age_gr_69+`)/1e6
sum(proj_res3_ssp1$`age_gr_69+`)/1e6

e <- ggplot(proj_res3_ssp3)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp3$`TPOP_age_gr_69+`, proj_res3_ssp3$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP3")


proj_res3_ssp3_s <- proj_res3_ssp3 %>% group_by(STATE) %>% dplyr::summarise(`TPOP_age_gr_69+`=sum(`TPOP_age_gr_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

f <- ggplot(proj_res3_ssp3_s)+
  geom_sf(aes(fill=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`)/`TPOP_age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`TPOP_age_gr_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",    mape(proj_res3_ssp3_s$`TPOP_age_gr_69+`, proj_res3_ssp3_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP3")

sum(proj_res3_ssp3$`TPOP_age_gr_69+`)/1e6
sum(proj_res3_ssp3$`age_gr_69+`)/1e6

panel <- a + b + c + d + e + f + g + h +  plot_layout(ncol=2, guides = "collect") + plot_annotation(
  title = 'United States counties and states',
  theme = theme(plot.title = element_text(size = 15))) & xlab("") & ylab("")

ggsave("us_comparison_apes_2050.png", panel, height = 14, width = 11, scale=1.25)

###########
# also "validate"/compare with EU

eu_data <- read_tsv(gzfile("proj_19rp3.tsv.gz"))   

# merge to shapefile

eu_data <- dplyr::select(eu_data, 1, 52)
colnames(eu_data) <- c("var", "value")
eu_data <- filter(eu_data, grepl("BSL", eu_data$var))

# create age variable

val <- as.data.frame(do.call(rbind, strsplit(eu_data$var, ",")))
val <- dplyr::select(val, 2, 3, 5)
colnames(val) <- c("age", "gender", "NUTS_ID")

eu_data <- bind_cols(eu_data, val)

# exclude total (keep gender stratification)

eu_data <- filter(eu_data, eu_data$gender=="T")
eu_data$var <- NULL
eu_data$age <- gsub("Y", "", eu_data$age)
eu_data <- filter(eu_data, eu_data$age!="TOTAL" & eu_data$age!="_GE100" & eu_data$age!="_LT1")
eu_data$age <- as.numeric(eu_data$age)

eu_data$age <- ifelse(eu_data$age>0 & eu_data$age<20, "0-19", ifelse(eu_data$age>=20 & eu_data$age<70, "20-69", "69+"))

eu_data <- dplyr::group_by(eu_data, age, gender, NUTS_ID) %>% dplyr::summarise(value=sum(value, na.rm=T))

eu_data <- eu_data %>%
  pivot_wider(names_from = c(gender, age), values_from = c(value))

# sf

r <- read_sf("NUTS_RG_60M_2021_4326.geojson")

#

eu_data <-merge(eu_data, r, by="NUTS_ID") 
eu_data <- st_as_sf(eu_data)
st_crs(eu_data) <- 4326

#########################

# extract projections of pop

r1_ssp2_2050_e <- exact_extract(r1_ssp2_2050, eu_data, "sum")
colnames(r1_ssp5_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")

r1_ssp5_2050_e <- exact_extract(r1_ssp5_2050, eu_data, "sum")
colnames(r1_ssp5_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")

r1_ssp1_2050_e <- exact_extract(r1_ssp1_2050, eu_data, "sum")
colnames(r1_ssp3_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")

r1_ssp3_2050_e <- exact_extract(r1_ssp3_2050, eu_data, "sum")
colnames(r1_ssp3_2050_e) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")

#r1_ssp2_2050_e <- (r1_ssp2_2050_e + r1_ssp5_2050_e) / 2

eu_data_ssp1 <- bind_cols(eu_data, r1_ssp1_2050_e)
eu_data_ssp2 <- bind_cols(eu_data, r1_ssp2_2050_e)
eu_data_ssp3 <- bind_cols(eu_data, r1_ssp3_2050_e)
eu_data_ssp5 <- bind_cols(eu_data, r1_ssp5_2050_e)

colnames(eu_data_ssp1)[14:16] <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
colnames(eu_data_ssp2)[14:16] <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
colnames(eu_data_ssp3)[14:16] <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
colnames(eu_data_ssp5)[14:16] <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")

# calculate MAPEs and MSEs

mape <- function(real,modeled){
  
  round(mean(abs((real - modeled)/real), na.rm=T), 2)
  
}


mse <- function(real,modeled){
  
  round(mean((real - modeled)^2, na.rm=T), 2)
  
}

##

eu_data_ssp1 <- filter(eu_data_ssp1, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
eu_data_ssp1 <- filter(eu_data_ssp1, `T_69+`>0 & !is.na(`T_69+`))
eu_data_ssp2 <- filter(eu_data_ssp2, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
eu_data_ssp2 <- filter(eu_data_ssp2, `T_69+`>0 & !is.na(`T_69+`))
eu_data_ssp3 <- filter(eu_data_ssp3, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
eu_data_ssp3 <- filter(eu_data_ssp3, `T_69+`>0 & !is.na(`T_69+`))
eu_data_ssp5 <- filter(eu_data_ssp5, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
eu_data_ssp5 <- filter(eu_data_ssp5, `T_69+`>0 & !is.na(`T_69+`))


a <- ggplot()+
  geom_point(aes(y=eu_data_ssp1$`T_69+`, x=eu_data_ssp1$`age_gr_69+`), alpha=0.3) + 
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP1), tot. = ", round(sum(eu_data_ssp1$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Eurostat projections, tot. = ", round(sum(eu_data_ssp1$`T_69+`)/1e6, 2)))+
  annotate("text", x=1.2e+04, y=1e+05, label= paste0("r = ", round(cor(eu_data_ssp1$`T_69+`, eu_data_ssp1$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+04, y=0.8e+05, label= paste0("MAPE = ", mape(eu_data_ssp1$`T_69+`, eu_data_ssp1$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

b <- ggplot()+
  geom_point(aes(y=eu_data_ssp2$`T_69+`, x=eu_data_ssp2$`age_gr_69+`), alpha=0.3) + 
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP2), tot. = ", round(sum(eu_data_ssp2$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Eurostat projections, tot. = ", round(sum(eu_data_ssp2$`T_69+`)/1e6, 2)))+
  annotate("text", x=1.2e+04, y=1e+05, label= paste0("r = ", round(cor(eu_data_ssp2$`T_69+`, eu_data_ssp2$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+04, y=0.8e+05, label= paste0("MAPE = ", mape(eu_data_ssp2$`T_69+`, eu_data_ssp2$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

c <- ggplot()+
  geom_point(aes(y=eu_data_ssp3$`T_69+`, x=eu_data_ssp3$`age_gr_69+`), alpha=0.3) + 
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP3), tot. = ", round(sum(eu_data_ssp3$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Eurostat projections, tot. = ", round(sum(eu_data_ssp3$`T_69+`)/1e6, 2)))+
  annotate("text", x=1.2e+04, y=1e+05, label= paste0("r = ", round(cor(eu_data_ssp3$`T_69+`, eu_data_ssp3$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+04, y=0.8e+05, label= paste0("MAPE = ", mape(eu_data_ssp3$`T_69+`, eu_data_ssp3$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

d <- ggplot()+
  geom_point(aes(y=eu_data_ssp5$`T_69+`, x=eu_data_ssp5$`age_gr_69+`), alpha=0.3) + 
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP5), tot. = ", round(sum(eu_data_ssp5$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Eurostat projections, tot. = ", round(sum(eu_data_ssp5$`T_69+`)/1e6, 2)))+
  annotate("text", x=1.2e+04, y=1e+05, label= paste0("r = ", round(cor(eu_data_ssp5$`T_69+`, eu_data_ssp5$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+04, y=0.8e+05, label= paste0("MAPE = ", mape(eu_data_ssp5$`T_69+`, eu_data_ssp5$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")


a + b + c + d + plot_layout(ncol=2) + plot_annotation(
  title = 'Europe, NUTS-3 regions',
  theme = theme(plot.title = element_text(size = 15)))

ggsave("eu_comparison_2050.png", scale=1, height = 7, width = 11)

###

eu_data_m_ssp1 <- st_filter(eu_data_ssp1, bbox_to_poly(c(-20, 30, 40, 80), crs=4326))
eu_data_m_ssp2 <- st_filter(eu_data_ssp2, bbox_to_poly(c(-20, 30, 40, 80), crs=4326))
eu_data_m_ssp3 <- st_filter(eu_data_ssp3, bbox_to_poly(c(-20, 30, 40, 80), crs=4326))
eu_data_m_ssp5 <- st_filter(eu_data_ssp5, bbox_to_poly(c(-20, 30, 40, 80), crs=4326))

a <- ggplot(eu_data_m_ssp1)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-10, y=70, label= paste0("MAPE = ",  mape(eu_data_m_ssp1$`T_69+`, eu_data_m_ssp1$`age_gr_69+`)*100, " %"))+
  coord_sf(xlim=c(-20, 35), ylim=c(35, 70))

eu_data_m_s_ssp1 <- eu_data_m_ssp1 %>% group_by(CNTR_CODE) %>% dplyr::summarise(`T_69+`=sum(`T_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

b <- ggplot(eu_data_m_s_ssp1)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100, label=))+
  geom_sf_label_repel(aes(label=round(abs((`T_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=5, y=70, label= paste0("MAPE = ",  mape(eu_data_m_s_ssp1$`T_69+`, eu_data_m_s_ssp1$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")

sum(eu_data_m_s_ssp1$`T_69+`)/1e6
sum(eu_data_m_s_ssp1$`age_gr_69+`)/1e6

c <- ggplot(eu_data_m_ssp2)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-10, y=70, label= paste0("MAPE = ",  mape(eu_data_m_ssp2$`T_69+`, eu_data_m_ssp2$`age_gr_69+`)*100, " %"))+
  coord_sf(xlim=c(-20, 35), ylim=c(35, 70))

eu_data_m_s_ssp2 <- eu_data_m_ssp2 %>% group_by(CNTR_CODE) %>% dplyr::summarise(`T_69+`=sum(`T_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

d <- ggplot(eu_data_m_s_ssp2)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100, label=))+
  geom_sf_label_repel(aes(label=round(abs((`T_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=5, y=70, label= paste0("MAPE = ",  mape(eu_data_m_s_ssp2$`T_69+`, eu_data_m_s_ssp2$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")

sum(eu_data_m_s_ssp2$`T_69+`)/1e6
sum(eu_data_m_s_ssp2$`age_gr_69+`)/1e6

e <- ggplot(eu_data_m_ssp3)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-10, y=70, label= paste0("MAPE = ",  mape(eu_data_m_ssp3$`T_69+`, eu_data_m_ssp3$`age_gr_69+`)*100, " %"))+
  coord_sf(xlim=c(-20, 35), ylim=c(35, 70))

eu_data_m_s_ssp3 <- eu_data_m_ssp3 %>% group_by(CNTR_CODE) %>% dplyr::summarise(`T_69+`=sum(`T_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

f <- ggplot(eu_data_m_s_ssp3)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100, label=))+
  geom_sf_label_repel(aes(label=round(abs((`T_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=5, y=70, label= paste0("MAPE = ",  mape(eu_data_m_s_ssp3$`T_69+`, eu_data_m_s_ssp3$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")

sum(eu_data_m_s_ssp3$`T_69+`)/1e6
sum(eu_data_m_s_ssp3$`age_gr_69+`)/1e6

g <- ggplot(eu_data_m_ssp5)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-10, y=70, label= paste0("MAPE = ",  mape(eu_data_m_ssp5$`T_69+`, eu_data_m_ssp5$`age_gr_69+`)*100, " %"))+
  coord_sf(xlim=c(-20, 35), ylim=c(35, 70))

eu_data_m_s_ssp5 <- eu_data_m_ssp5 %>% group_by(CNTR_CODE) %>% dplyr::summarise(`T_69+`=sum(`T_69+`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

h <- ggplot(eu_data_m_s_ssp5)+
  geom_sf(aes(fill=round(abs((`T_69+` - `age_gr_69+`)/`T_69+`), 2)*100, label=))+
  geom_sf_label_repel(aes(label=round(abs((`T_69+` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=5, y=70, label= paste0("MAPE = ",  mape(eu_data_m_s_ssp5$`T_69+`, eu_data_m_s_ssp5$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")

sum(eu_data_m_s_ssp5$`T_69+`)/1e6
sum(eu_data_m_s_ssp5$`age_gr_69+`)/1e6

a + b + c + d + e + f + g + h + plot_layout(ncol=2, guides = "collect") + plot_annotation(
  title = 'EU NUTS and countries',
  theme = theme(plot.title = element_text(size = 15))) & xlab("") & ylab("")

ggsave("eu_comparison_apes_2050.png", height = 14, width = 8, scale=1.25)

##

# ###################
# 
r1_g <- stack(paste0(getwd(), "/r1_2020.tif"))
r1_g_245 <- stack(paste0(getwd(), "/r1_g_245_2050_final.tif"))
r1_g_585 <- stack(paste0(getwd(), "/r1_g_585_2050_final.tif"))
r1_g_126 <- stack(paste0(getwd(), "/r1_g_126_2050_final.tif"))
r1_g_370 <- stack(paste0(getwd(), "/r1_g_370_2050_final.tif"))

# validate for US

r1_ssp5_2050 <- r1_g_585
r1_ssp2_2050 <- r1_g_245 
r1_ssp1_2050 <- r1_g_126
r1_ssp3_2050 <- r1_g_370 

################

#1) # usa 2nd source

usa2 <- read_csv("new_comparison_data/SSP_asrc/SSP_asrc.csv")

# create age group counts

usa2_age_gr_69 <- dplyr::filter(usa2, AGE>14)
usa2_age_gr_69 <- dplyr::filter(usa2_age_gr_69, YEAR==2050)
usa2_age_gr_69 <- dplyr::group_by(usa2_age_gr_69, GEOID) %>% dplyr::summarise(SSP1=sum(SSP1, na.rm=T), SSP2=sum(SSP2, na.rm=T), SSP3=sum(SSP3, na.rm=T), SSP5=sum(SSP5, na.rm=T))

# merge 

fipslist <- read_csv(file="https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt", col_names = FALSE) %>%
  mutate(GEOID = paste0(X2, X3))

usa2_age_gr_69 <- merge(usa2_age_gr_69, fipslist, "GEOID")

counties_sf <- get_urbn_map("counties", sf = TRUE)

usa2_age_gr_69 <- merge(usa2_age_gr_69, counties_sf, by.x=c("X4", "X1"),  by.y=c("county_name", "state_abbv"))

usa2_age_gr_69 <- st_as_sf(usa2_age_gr_69)

usa2_age_gr_69 <- st_transform(usa2_age_gr_69, 4326)

usa2_age_gr_69_ssp5_2050 <- exact_extract(r1_ssp5_2050, usa2_age_gr_69, "sum", max_cells_in_memory=1e09)
usa2_age_gr_69_ssp2_2050 <- exact_extract(r1_ssp2_2050,usa2_age_gr_69 , "sum", max_cells_in_memory=1e09)
usa2_age_gr_69_ssp1_2050 <- exact_extract(r1_ssp1_2050, usa2_age_gr_69, "sum", max_cells_in_memory=1e09)
usa2_age_gr_69_ssp3_2050 <- exact_extract(r1_ssp3_2050,usa2_age_gr_69, "sum", max_cells_in_memory=1e09)

# merge back to data
colnames(usa2_age_gr_69_ssp5_2050) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
usa2_age_gr_69_ssp5_2050 <- bind_cols(usa2_age_gr_69, usa2_age_gr_69_ssp5_2050)

colnames(usa2_age_gr_69_ssp2_2050) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
usa2_age_gr_69_ssp2_2050 <- bind_cols(usa2_age_gr_69, usa2_age_gr_69_ssp2_2050)

colnames(usa2_age_gr_69_ssp1_2050) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
usa2_age_gr_69_ssp1_2050 <- bind_cols(usa2_age_gr_69, usa2_age_gr_69_ssp1_2050)

colnames(usa2_age_gr_69_ssp3_2050) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
usa2_age_gr_69_ssp3_2050 <- bind_cols(usa2_age_gr_69, usa2_age_gr_69_ssp3_2050)

###

# calculate MAPEs and MSEs

mape <- function(real,modeled){
  
  round(mean(abs((real - modeled)/real), na.rm=T), 2)
  
}


mse <- function(real,modeled){
  
  round(mean((real - modeled)^2, na.rm=T), 2)
  
}


usa2_age_gr_69_ssp1_2050 <- filter(usa2_age_gr_69_ssp1_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
usa2_age_gr_69_ssp1_2050 <- filter(usa2_age_gr_69_ssp1_2050, SSP1>0 & !is.na(SSP1))
usa2_age_gr_69_ssp2_2050 <- filter(usa2_age_gr_69_ssp2_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
usa2_age_gr_69_ssp2_2050 <- filter(usa2_age_gr_69_ssp2_2050, SSP2>0 & !is.na(SSP2))
usa2_age_gr_69_ssp3_2050 <- filter(usa2_age_gr_69_ssp3_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
usa2_age_gr_69_ssp3_2050 <- filter(usa2_age_gr_69_ssp3_2050, SSP3>0 & !is.na(SSP3))
usa2_age_gr_69_ssp5_2050 <- filter(usa2_age_gr_69_ssp5_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
usa2_age_gr_69_ssp5_2050 <- filter(usa2_age_gr_69_ssp5_2050, SSP5>0 & !is.na(SSP5))

a<- ggplot()+
  geom_point(aes(y=usa2_age_gr_69_ssp1_2050$SSP1, x=usa2_age_gr_69_ssp1_2050$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP1), tot. = ", round(sum(usa2_age_gr_69_ssp1_2050$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Hauer (2050, SSP1), tot. = ", round(sum(usa2_age_gr_69_ssp1_2050$SSP1)/1e6, 2)))+
  annotate("text", x=1.2e+03, y=1e+05, label= paste0("r = ", round(cor(usa2_age_gr_69_ssp1_2050$SSP1, usa2_age_gr_69_ssp1_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.5e+03, y=0.7e+05, label= paste0("MAPE = ", mape(usa2_age_gr_69_ssp1_2050$SSP1, usa2_age_gr_69_ssp1_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

b<- ggplot()+
  geom_point(aes(y=usa2_age_gr_69_ssp2_2050$SSP2, x=usa2_age_gr_69_ssp2_2050$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP2), tot. = ", round(sum(usa2_age_gr_69_ssp2_2050$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Hauer (2050, SSP2), tot. = ", round(sum(usa2_age_gr_69_ssp2_2050$SSP2)/1e6, 2)))+
  annotate("text", x=1.2e+03, y=1e+05, label= paste0("r = ", round(cor(usa2_age_gr_69_ssp2_2050$SSP2, usa2_age_gr_69_ssp2_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.5e+03, y=0.7e+05, label= paste0("MAPE = ", mape(usa2_age_gr_69_ssp2_2050$SSP2, usa2_age_gr_69_ssp2_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

c<- ggplot()+
  geom_point(aes(y=usa2_age_gr_69_ssp3_2050$SSP3, x=usa2_age_gr_69_ssp3_2050$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP3), tot. = ", round(sum(usa2_age_gr_69_ssp3_2050$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Hauer (2050, SSP3), tot. = ", round(sum(usa2_age_gr_69_ssp3_2050$SSP3)/1e6, 2)))+
  annotate("text", x=1.2e+03, y=1e+05, label= paste0("r = ", round(cor(usa2_age_gr_69_ssp3_2050$SSP3, usa2_age_gr_69_ssp3_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.5e+03, y=0.7e+05, label= paste0("MAPE = ", mape(usa2_age_gr_69_ssp3_2050$SSP3, usa2_age_gr_69_ssp3_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

d<- ggplot()+
  geom_point(aes(y=usa2_age_gr_69_ssp5_2050$SSP5, x=usa2_age_gr_69_ssp5_2050$`age_gr_69+`), alpha=0.1) + 
  scale_x_log10(limits= c(1000, 100000))+
  scale_y_log10(limits= c(1000, 100000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP5), tot. = ", round(sum(usa2_age_gr_69_ssp5_2050$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("Hauer (2050, SSP5), tot. = ", round(sum(usa2_age_gr_69_ssp5_2050$SSP5)/1e6, 2)))+
  annotate("text", x=1.2e+03, y=1e+05, label= paste0("r = ", round(cor(usa2_age_gr_69_ssp5_2050$SSP5, usa2_age_gr_69_ssp5_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.5e+03, y=0.7e+05, label= paste0("MAPE = ", mape(usa2_age_gr_69_ssp5_2050$SSP5, usa2_age_gr_69_ssp5_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

box <- a + b + c + d + plot_layout(ncol=2) + plot_annotation(
  title = 'United States counties',
  theme = theme(plot.title = element_text(size = 15)))

ggsave("us_comparison_hauer_2050.png", box, height = 8, width = 11)

###

c <- ggplot(usa2_age_gr_69_ssp2_2050)+
  geom_sf(aes(fill=round(abs((SSP2 - `age_gr_69+`)/SSP2), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(usa2_age_gr_69_ssp2_2050$SSP2, usa2_age_gr_69_ssp2_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP2")

usa2_age_gr_69_ssp2_2050_s <- usa2_age_gr_69_ssp2_2050 %>% group_by(state_fips) %>% dplyr::summarise(SSP2=sum(SSP2, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

d <- ggplot(usa2_age_gr_69_ssp2_2050_s)+
  geom_sf(aes(fill=round(abs((SSP2 - `age_gr_69+`)/SSP2), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((SSP2 - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",   mape(usa2_age_gr_69_ssp2_2050_s$SSP2, usa2_age_gr_69_ssp2_2050_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP2")

sum(usa2_age_gr_69_ssp2_2050$SSP2)/1e6
sum(usa2_age_gr_69_ssp2_2050$`age_gr_69+`)/1e6

g <- ggplot(usa2_age_gr_69_ssp5_2050)+
  geom_sf(aes(fill=round(abs((SSP5 - `age_gr_69+`)/SSP5), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(usa2_age_gr_69_ssp5_2050$SSP5, usa2_age_gr_69_ssp5_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP5")

usa2_age_gr_69_ssp5_2050_s <- usa2_age_gr_69_ssp5_2050 %>% group_by(state_fips) %>% dplyr::summarise(SSP5=sum(SSP5, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

h <- ggplot(usa2_age_gr_69_ssp5_2050_s)+
  geom_sf(aes(fill=round(abs((SSP5 - `age_gr_69+`)/SSP5), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((SSP5 - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",   mape(usa2_age_gr_69_ssp5_2050_s$SSP5, usa2_age_gr_69_ssp5_2050_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP5")

sum(usa2_age_gr_69_ssp5_2050$SSP5)/1e6
sum(usa2_age_gr_69_ssp5_2050$`age_gr_69+`)/1e6

a <- ggplot(usa2_age_gr_69_ssp1_2050)+
  geom_sf(aes(fill=round(abs((SSP1 - `age_gr_69+`)/SSP1), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(usa2_age_gr_69_ssp1_2050$SSP1, usa2_age_gr_69_ssp1_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP1")

usa2_age_gr_69_ssp1_2050_s <- usa2_age_gr_69_ssp1_2050 %>% group_by(state_fips) %>% dplyr::summarise(SSP1=sum(SSP1, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

b <- ggplot(usa2_age_gr_69_ssp1_2050_s)+
  geom_sf(aes(fill=round(abs((SSP1 - `age_gr_69+`)/SSP1), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((SSP1 - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",   mape(usa2_age_gr_69_ssp1_2050_s$SSP1, usa2_age_gr_69_ssp1_2050_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP1")

sum(usa2_age_gr_69_ssp1_2050$SSP1)/1e6
sum(usa2_age_gr_69_ssp1_2050$`age_gr_69+`)/1e6

e <- ggplot(usa2_age_gr_69_ssp3_2050)+
  geom_sf(aes(fill=round(abs((SSP3 - `age_gr_69+`)/SSP3), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-115, y=52, label= paste0("MAPE = ",    mape(usa2_age_gr_69_ssp3_2050$SSP3, usa2_age_gr_69_ssp3_2050$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP3")

usa2_age_gr_69_ssp3_2050_s <- usa2_age_gr_69_ssp3_2050 %>% group_by(state_fips) %>% dplyr::summarise(SSP3=sum(SSP3, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

f <- ggplot(usa2_age_gr_69_ssp3_2050_s)+
  geom_sf(aes(fill=round(abs((SSP3 - `age_gr_69+`)/SSP3), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((SSP3 - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-105, y=52, label= paste0("MAPE = ",   mape(usa2_age_gr_69_ssp3_2050_s$SSP3, usa2_age_gr_69_ssp3_2050_s$`age_gr_69+`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP3")

sum(usa2_age_gr_69_ssp3_2050$SSP3)/1e6
sum(usa2_age_gr_69_ssp3_2050$`age_gr_69+`)/1e6

box <- a + b + c + d + e + f + g + h +  plot_layout(ncol=2, guides = "collect") + plot_annotation(
  title = 'United States counties and states',
  theme = theme(plot.title = element_text(size = 15))) & xlab("") & ylab("")

ggsave("us_comparison_hauer_apes_2050.png", box, height = 14, width = 11, scale=1.25)


################

#2) 

uk <- read.csv("new_comparison_data/tot_pop_uk.csv", sep=";", dec=".")

# create age group counts

uk$AGE.GROUP <- ifelse(uk$AGE.GROUP %in% unique(uk$AGE.GROUP)[1:4], "age_gr_0_19", ifelse(uk$AGE.GROUP %in% unique(uk$AGE.GROUP)[5:14], "age_gr_20_69", ifelse(uk$AGE.GROUP %in% unique(uk$AGE.GROUP)[15:19], "age_gr_69+", NA)))

uk <- na.omit(uk)
uk <- dplyr::filter(uk, AREA!="England")

# merge to a shapefile
uk <- dplyr::select(uk , 1, 2, 3, 29)

uk$AREA <- substr(uk$CODE, 1, 4)

uk$X2043 <- gsub("\\.", "", uk$X2043)
uk$X2043 <- as.numeric(uk$X2043)

uk <- uk %>% group_by(AREA, CODE, AGE.GROUP) %>% dplyr::summarise(value=sum(X2043, na.rm=T))

uk <- pivot_wider(uk,  names_from = AGE.GROUP, values_from = value)

uk_sf <- read_sf("new_comparison_data/Local_Authority_Districts_May_2022_UK_BFE_V3_2022_-4981034458221456209.gpkg")

uk2 <- merge(uk, uk_sf %>% dplyr::select(LAD22CD), by.x="CODE", by.y="LAD22CD", all.y=T)

uk3 <- st_as_sf(uk2)

uk3 <- st_transform(uk3, 4326)

uk3_ssp5_2050 <- exact_extract(r1_ssp5_2050, uk3, "sum", max_cells_in_memory=1e09)
uk3_ssp2_2050 <- exact_extract(r1_ssp2_2050,uk3  , "sum", max_cells_in_memory=1e09)
uk3_ssp1_2050 <- exact_extract(r1_ssp1_2050, uk3, "sum", max_cells_in_memory=1e09)
uk3_ssp3_2050 <- exact_extract(r1_ssp3_2050,uk3 , "sum", max_cells_in_memory=1e09)

# merge back to data
colnames(uk3_ssp5_2050) <- c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
uk3_ssp5_2050 <- bind_cols(uk3 , uk3_ssp5_2050)

colnames(uk3_ssp2_2050) <-  c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
uk3_ssp2_2050 <- bind_cols(uk3  , uk3_ssp2_2050)

colnames(uk3_ssp1_2050) <- c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
uk3_ssp1_2050 <- bind_cols(uk3, uk3_ssp1_2050)

colnames(uk3_ssp3_2050) <- c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
uk3_ssp3_2050 <- bind_cols(uk3, uk3_ssp3_2050)

uk3_ssp1_2050 <- filter(uk3_ssp1_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
uk3_ssp1_2050 <- filter(uk3_ssp1_2050, `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))
uk3_ssp2_2050 <- filter(uk3_ssp2_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
uk3_ssp2_2050 <- filter(uk3_ssp2_2050, `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))
uk3_ssp3_2050 <- filter(uk3_ssp3_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
uk3_ssp3_2050 <- filter(uk3_ssp3_2050, `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))
uk3_ssp5_2050 <- filter(uk3_ssp5_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
uk3_ssp5_2050 <- filter(uk3_ssp5_2050, `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))

a<- ggplot()+
  geom_point(aes(y=uk3_ssp1_2050$`age_gr_69+_mod`, x=uk3_ssp1_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  xlab(paste0("This study (2050, SSP1), tot. = ", round(sum(uk3_ssp1_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Nash/ONS (2043), tot. = ", round(sum(uk3_ssp1_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.1e+04, y=1e+05, label= paste0("r = ", round(cor(uk3_ssp1_2050$`age_gr_69+_mod`, uk3_ssp1_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.25e+04, y=0.8e+05, label= paste0("MAPE = ", mape(uk3_ssp1_2050$`age_gr_69+`, uk3_ssp1_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")

b<- ggplot()+
  geom_point(aes(y=uk3_ssp2_2050$`age_gr_69+_mod`, x=uk3_ssp2_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  xlab(paste0("This study (2050, SSP2), tot. = ", round(sum(uk3_ssp2_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Nash/ONS (2043), tot. = ", round(sum(uk3_ssp2_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.1e+04, y=1e+05, label= paste0("r = ", round(cor(uk3_ssp2_2050$`age_gr_69+_mod`, uk3_ssp2_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.25e+04, y=0.8e+05, label= paste0("MAPE = ", mape(uk3_ssp2_2050$`age_gr_69+`, uk3_ssp2_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")

c<- ggplot()+
  geom_point(aes(y=uk3_ssp3_2050$`age_gr_69+_mod`, x=uk3_ssp3_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  xlab(paste0("This study (2050, SSP3), tot. = ", round(sum(uk3_ssp3_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Nash/ONS (2043), tot. = ", round(sum(uk3_ssp3_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.1e+04, y=1e+05, label= paste0("r = ", round(cor(uk3_ssp3_2050$`age_gr_69+_mod`, uk3_ssp3_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.25e+04, y=0.8e+05, label= paste0("MAPE = ", mape(uk3_ssp3_2050$`age_gr_69+`, uk3_ssp3_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")

d<- ggplot()+
  geom_point(aes(y=uk3_ssp5_2050$`age_gr_69+_mod`, x=uk3_ssp5_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(10000, 100000))+
  scale_y_log10(limits= c(10000, 100000))+
  xlab(paste0("This study (2050, SSP5), tot. = ", round(sum(uk3_ssp5_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Nash/ONS (2043), tot. = ", round(sum(uk3_ssp5_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.1e+04, y=1e+05, label= paste0("r = ", round(cor(uk3_ssp5_2050$`age_gr_69+_mod`, uk3_ssp5_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.25e+04, y=0.8e+05, label= paste0("MAPE = ", mape(uk3_ssp5_2050$`age_gr_69+`, uk3_ssp5_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")


box <- a + b + c + d + plot_layout(ncol=2) + plot_annotation(
  title = 'UK local authorities',
  theme = theme(plot.title = element_text(size = 15)))


ggsave("uk_comparison_2050.png", box, height = 8, width = 11)

####

c <- ggplot(uk3_ssp2_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-4, y=56, label= paste0("MAPE = ",   mape(uk3_ssp2_2050$`age_gr_69+`, uk3_ssp2_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP2")+
  xlab("")+
  ylab("")

uk3_ssp2_2050_s <- uk3_ssp2_2050 %>% group_by(AREA) %>% dplyr::summarise(`age_gr_69+_mod`=sum(`age_gr_69+_mod`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))

uk3_ssp2_2050_s <- st_as_sf(uk3_ssp2_2050_s)
uk3_ssp2_2050_s <- st_transform(uk3_ssp2_2050_s, 4326)
st_geometry(uk3_ssp2_2050_s) = "geometry"

d <- ggplot(uk3_ssp2_2050_s)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-2, y=56, label= paste0("MAPE = ",   mape(uk3_ssp2_2050_s$`age_gr_69+`, uk3_ssp2_2050_s$`age_gr_69+_mod`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP2")

g <- ggplot(uk3_ssp5_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-4, y=56, label= paste0("MAPE = ",   mape(uk3_ssp5_2050$`age_gr_69+`, uk3_ssp5_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP5")+
  xlab("")+
  ylab("")

uk3_ssp5_2050_s <- uk3_ssp5_2050 %>% group_by(AREA) %>% dplyr::summarise(`age_gr_69+_mod`=sum(`age_gr_69+_mod`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))
uk3_ssp5_2050_s <- st_transform(uk3_ssp5_2050_s, 4326)
st_geometry(uk3_ssp5_2050_s) = "geometry"

h <- ggplot(uk3_ssp5_2050_s)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-2, y=56, label= paste0("MAPE = ",   mape(uk3_ssp5_2050_s$`age_gr_69+`, uk3_ssp5_2050_s$`age_gr_69+_mod`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP5")

a <- ggplot(uk3_ssp1_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-4, y=56, label= paste0("MAPE = ",   mape(uk3_ssp1_2050$`age_gr_69+`, uk3_ssp1_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP1")+
  xlab("")+
  ylab("")

uk3_ssp1_2050_s <- uk3_ssp1_2050 %>% group_by(AREA) %>% dplyr::summarise(`age_gr_69+_mod`=sum(`age_gr_69+_mod`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))
uk3_ssp1_2050_s <- st_transform(uk3_ssp1_2050_s, 4326)
st_geometry(uk3_ssp1_2050_s) = "geometry"

b <- ggplot(uk3_ssp1_2050_s)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-2, y=56, label= paste0("MAPE = ",   mape(uk3_ssp1_2050_s$`age_gr_69+`, uk3_ssp1_2050_s$`age_gr_69+_mod`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP1")

e <- ggplot(uk3_ssp3_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=-4, y=56, label= paste0("MAPE = ",   mape(uk3_ssp3_2050$`age_gr_69+`, uk3_ssp3_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP3")+
  xlab("")+
  ylab("")

uk3_ssp3_2050_s <- uk3_ssp3_2050 %>% group_by(AREA) %>% dplyr::summarise(`age_gr_69+_mod`=sum(`age_gr_69+_mod`, na.rm=T), `age_gr_69+`=sum(`age_gr_69+`, na.rm=T))
uk3_ssp3_2050_s <- st_transform(uk3_ssp3_2050_s, 4326)
st_geometry(uk3_ssp3_2050_s) = "geometry"

f <- ggplot(uk3_ssp3_2050_s)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  annotate("text", x=-2, y=56, label= paste0("MAPE = ",   mape(uk3_ssp3_2050_s$`age_gr_69+`, uk3_ssp3_2050_s$`age_gr_69+_mod`)*100, " %"))+
  xlab("")+
  ylab("")+
  ggtitle("SSP3")

box <- a + b + c + d + e + f + g + h +  plot_layout(ncol=2, guides = "collect") + plot_annotation(
  title = 'UK local authorities and regions',
  theme = theme(plot.title = element_text(size = 15))) & xlab("") & ylab("")

ggsave("uk_comparison_apes_2050.png", box, height = 14, width = 7, scale=1.35)


################

#3) 

china <- list.files("new_comparison_data/china", recursive = T, full.names = T, pattern = "2050")
sheets = excel_sheets(china[1])
china <- lapply(china, function(X){lapply(sheets, function(Y){r <- read_xlsx(X, sheet=Y); r$province <- Y; r$ssp <- substr(basename(X), 5, 8); r})})
china2 <- bind_rows(china)

china2$value <- rowSums(china2[,c(2:15)], na.rm = T) 
china2[,c(2:15)] <- NULL
china2 <- filter(china2, ...1!="总人数")

#china2 %>% group_by(ssp) %>% dplyr::summarise(value=sum(value, na.rm=T))

# create age group counts

china2$...1 <- ifelse(china2$...1=="100+", "100", china2$...1)
china2$...1 <- as.numeric(as.character(china2$...1))

china2$agegr <- NA
china2$agegr[china2$...1>=0 & china2$...1<20] <- "age_gr_0_19"
china2$agegr[china2$...1>=20 & china2$...1<69] <- "age_gr_20_69"
china2$agegr[china2$...1>=69] <- "age_gr_69+"

china2$province <- ifelse(china2$province=="Ningxia", "Ningxia Hui", china2$province)
china2$province <- ifelse(china2$province=="NeiMonggol", "Nei Mongol", china2$province)
china2$province <- ifelse(china2$province=="Xinjiang", "Xinjiang Uygur", china2$province)

china3 <- china2 %>% group_by(province, ssp, agegr) %>% dplyr::summarise(value=sum(value, na.rm=T))

china3 <- pivot_wider(china3,  names_from = agegr, values_from = value)


# merge to a shapefile

china_sf <- st_as_sf(raster::getData('GADM', country='CHN', level=1))

china3 <- merge(china3, china_sf, by.x="province", by.y="NAME_1")

china3 <- st_as_sf(china3)

china3_ssp5_2050 <- exact_extract(r1_ssp5_2050, china3 %>% dplyr::filter(ssp=="SSP5"), "sum", max_cells_in_memory=1e09)
china3_ssp2_2050 <- exact_extract(r1_ssp2_2050,china3  %>% dplyr::filter(ssp=="SSP2") , "sum", max_cells_in_memory=1e09)
china3_ssp1_2050 <- exact_extract(r1_ssp1_2050, china3  %>% dplyr::filter(ssp=="SSP1"), "sum", max_cells_in_memory=1e09)
china3_ssp3_2050 <- exact_extract(r1_ssp3_2050,china3  %>% dplyr::filter(ssp=="SSP3"), "sum", max_cells_in_memory=1e09)

# merge back to data
colnames(china3_ssp5_2050) <- c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
china3_ssp5_2050 <- bind_cols(china3  %>% dplyr::filter(ssp=="SSP5"), china3_ssp5_2050)

colnames(china3_ssp2_2050) <-  c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
china3_ssp2_2050 <- bind_cols(china3  %>% dplyr::filter(ssp=="SSP2"), china3_ssp2_2050)

colnames(china3_ssp1_2050) <- c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
china3_ssp1_2050 <- bind_cols(china3 %>% dplyr::filter(ssp=="SSP1"), china3_ssp1_2050)

colnames(china3_ssp3_2050) <- c("age_gr_0_19_mod", "age_gr_20_69_mod", "age_gr_69+_mod")
china3_ssp3_2050 <- bind_cols(china3 %>% dplyr::filter(ssp=="SSP3"), china3_ssp3_2050)

###

# calculate MAPEs and MSEs

mape <- function(real,modeled){
  
  round(mean(abs((real - modeled)/real), na.rm=T), 2)
  
}


mse <- function(real,modeled){
  
  round(mean((real - modeled)^2, na.rm=T), 2)
  
}


china3_ssp1_2050 <- filter(china3_ssp1_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`) & `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))
china3_ssp2_2050 <- filter(china3_ssp2_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`) & `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))
china3_ssp3_2050 <- filter(china3_ssp3_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`) &  `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))
china3_ssp5_2050 <- filter(china3_ssp5_2050, `age_gr_69+`>0 & !is.na(`age_gr_69+`) &  `age_gr_69+_mod`>0 & !is.na(`age_gr_69+_mod`))

a<- ggplot()+
  geom_point(aes(y=china3_ssp1_2050$`age_gr_69+_mod`, x=china3_ssp1_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(1000000, 30000000))+
  scale_y_log10(limits= c(1000000, 30000000))+
  xlab(paste0("This study (2050, SSP1), tot. = ", round(sum(china3_ssp1_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Chen (2050, SSP1), tot. = ", round(sum(china3_ssp1_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.15e+06, y=2.75e+07, label= paste0("r = ", round(cor(china3_ssp1_2050$`age_gr_69+_mod`, china3_ssp1_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+06, y=2e+07, label= paste0("MAPE = ", mape(china3_ssp1_2050$`age_gr_69+`, china3_ssp1_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")


b <- ggplot()+
  geom_point(aes(y=china3_ssp2_2050$`age_gr_69+_mod`, x=china3_ssp2_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(1000000, 30000000))+
  scale_y_log10(limits= c(1000000, 30000000))+
  xlab(paste0("This study (2050, SSP2), tot. = ", round(sum(china3_ssp2_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Chen (2050, SSP2), tot. = ", round(sum(china3_ssp2_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.15e+06, y=2.75e+07, label= paste0("r = ", round(cor(china3_ssp2_2050$`age_gr_69+_mod`, china3_ssp2_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+06, y=2e+07, label= paste0("MAPE = ", mape(china3_ssp2_2050$`age_gr_69+`, china3_ssp2_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")


c<- ggplot()+
  geom_point(aes(y=china3_ssp3_2050$`age_gr_69+_mod`, x=china3_ssp3_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(1000000, 30000000))+
  scale_y_log10(limits= c(1000000, 30000000))+
  xlab(paste0("This study (2050, SSP3), tot. = ", round(sum(china3_ssp3_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Chen (2050, SSP3), tot. = ", round(sum(china3_ssp3_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.15e+06, y=2.75e+07, label= paste0("r = ", round(cor(china3_ssp3_2050$`age_gr_69+_mod`, china3_ssp3_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+06, y=2e+07, label= paste0("MAPE = ", mape(china3_ssp3_2050$`age_gr_69+`, china3_ssp3_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")

d<- ggplot()+
  geom_point(aes(y=china3_ssp5_2050$`age_gr_69+_mod`, x=china3_ssp5_2050$`age_gr_69+`), alpha=0.75) + 
  geom_abline()+
  scale_x_log10(limits= c(1000000, 30000000))+
  scale_y_log10(limits= c(1000000, 30000000))+
  xlab(paste0("This study (2050, SSP5), tot. = ", round(sum(china3_ssp5_2050$`age_gr_69+_mod`)/1e6, 2)))+
  ylab(paste0("Chen (2050, SSP5), tot. = ", round(sum(china3_ssp5_2050$`age_gr_69+`)/1e6, 2)))+
  annotate("text", x=1.15e+06, y=2.75e+07, label= paste0("r = ", round(cor(china3_ssp5_2050$`age_gr_69+_mod`, china3_ssp5_2050$`age_gr_69+`), 2)))+
  annotate("text", x=1.35e+06, y=2e+07, label= paste0("MAPE = ", mape(china3_ssp5_2050$`age_gr_69+`, china3_ssp5_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("Age group 69+")

box <- a + b + c + d + plot_layout(ncol=2) + plot_annotation(
  title = 'China provinces',
  theme = theme(plot.title = element_text(size = 15)))

ggsave("china_comparison_2050.png", box, height = 8, width = 11)

####

b <- ggplot(china3_ssp2_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=80, y=52, label= paste0("MAPE = ",    mape(china3_ssp2_2050$`age_gr_69+`, china3_ssp2_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP2")+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  xlab("")+
  ylab("")
  

a <- ggplot(china3_ssp1_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=80, y=52, label= paste0("MAPE = ",    mape(china3_ssp1_2050$`age_gr_69+`, china3_ssp1_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP1")+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  xlab("")+
  ylab("")
  

c <- ggplot(china3_ssp3_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=80, y=52, label= paste0("MAPE = ",    mape(china3_ssp3_2050$`age_gr_69+`, china3_ssp3_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP3")+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  xlab("")+
  ylab("")
  

d <- ggplot(china3_ssp5_2050)+
  geom_sf(aes(fill=round(abs((`age_gr_69+` - `age_gr_69+_mod`)/`age_gr_69+`), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((`age_gr_69+_mod` - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=80, y=52, label= paste0("MAPE = ",    mape(china3_ssp5_2050$`age_gr_69+`, china3_ssp5_2050$`age_gr_69+_mod`)*100, " %"))+
  ggtitle("SSP5")+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  xlab("")+
  ylab("")
  
box <- a + b + c + d +  plot_layout(ncol=2, guides = "collect") + plot_annotation(
  title = 'China provinces',
  theme = theme(plot.title = element_text(size = 15))) & xlab("") & ylab("")

ggsave("china_comparison_apes_2050.png", box, height = 7, width = 11, scale=1.25)


#4) 

india <- read.csv("new_comparison_data/Sub-national_Population_and_Human_Capital_Projections_India/files_to_share/popprojPOP_share_SSP2_2017-03-07_0741.csv")
india <- dplyr::filter(india, region!="india" & age>65 & year==2051)
india <- dplyr::group_by(india, region) %>% dplyr::summarise(value=sum(value, na.rm=T))
colnames(india)[1] <- "NAME_1"

# merge to a shapefile

india_gadm <- st_as_sf(raster::getData('GADM', country='IND', level=1))

india <- fuzzyjoin::stringdist_inner_join(india, india_gadm, by="NAME_1", max_dist=3.5)
india <- st_as_sf(india)

###

india_out <- exact_extract(r1_ssp2_2050,india , "sum", max_cells_in_memory=1e09)

colnames(india_out) <- c("age_gr_0_19", "age_gr_20_69", "age_gr_69+")
india <- bind_cols(india, india_out)

india <- filter(india, `age_gr_69+`>0 & !is.na(`age_gr_69+`))
india <- filter(india, value>0 & !is.na(value))

a <- ggplot()+
  geom_point(aes(y=india$value, x=india$`age_gr_69+`), alpha=1) + 
  scale_x_log10(limits= c(10000, 10000000))+
  scale_y_log10(limits= c(10000, 10000000))+
  geom_abline()+
  xlab(paste0("This study (2050, SSP2), tot. = ", round(sum(india$`age_gr_69+`)/1e6, 2)))+
  ylab(paste0("KC (2050, SSP2), tot. = ", round(sum(india$value)/1e6, 2)))+
  annotate("text", x=1.5e+04, y=1e+07, label= paste0("r = ", round(cor(india$value, india$`age_gr_69+`), 2)))+
  annotate("text", x=2.2e+04, y=7e+06, label= paste0("MAPE = ", mape(india$value, india$`age_gr_69+`)*100, " %"))+
  ggtitle("Age group 69+")

ggsave("india_comparison_kc_2050.png", height = 5, width = 5)

a <- ggplot(india)+
  geom_sf(aes(fill=round(abs((value - `age_gr_69+`)/value), 2)*100))+
  geom_sf_label_repel(aes(label=round(abs((value - `age_gr_69+`))/1e6, 2)))+
  scale_fill_stepsn(name="APE", colors = RColorBrewer::brewer.pal(5, "Reds"), breaks=c(0, 5, 10, 25, 50, 100), limits=c(0, 100))+
  annotate("text", x=73.5, y=37, label= paste0("MAPE = ",    mape(india$value, india$`age_gr_69+`)*100, " %"))+
  ggtitle("SSP2")+
  labs(caption = "Label: absolute difference (million people aged >69)")+
  xlab("")+
  ylab("")

ggsave("india_comparison_kc_apes_2050.png", height = 5, width = 5, scale=1.35)

###

mape <- function(real,modeled){
  
  mean(abs((real - modeled)/real), na.rm=T)
  
}

# country # source # scenarios # mape (%) # sum of error (million people aged 69+)

###

d1 <- mape(proj_res3_ssp1$`TPOP_age_gr_69+`, proj_res3_ssp1$`age_gr_69+`)*100
d12 <- round(abs((sum(proj_res3_ssp1$`TPOP_age_gr_69+`) - sum(proj_res3_ssp1$`age_gr_69+`)) / sum(proj_res3_ssp1$`TPOP_age_gr_69+`))*100, 2)
d2 <- round(abs((sum(proj_res3_ssp1$`TPOP_age_gr_69+`) - sum(proj_res3_ssp1$`age_gr_69+`)) / 1e6), 2)

d3 <- mape(proj_res3_ssp2$`TPOP_age_gr_69+`, proj_res3_ssp2$`age_gr_69+`)*100
d34 <- round(abs((sum(proj_res3_ssp2$`TPOP_age_gr_69+`) - sum(proj_res3_ssp2$`age_gr_69+`)) / sum(proj_res3_ssp2$`TPOP_age_gr_69+`))*100, 2)
d4 <- round(abs((sum(proj_res3_ssp2$`TPOP_age_gr_69+`) - sum(proj_res3_ssp2$`age_gr_69+`)) / 1e6), 2)

d5 <- mape(proj_res3_ssp3$`TPOP_age_gr_69+`, proj_res3_ssp3$`age_gr_69+`)*100
d56 <- round(abs((sum(proj_res3_ssp3$`TPOP_age_gr_69+`) - sum(proj_res3_ssp3$`age_gr_69+`)) / sum(proj_res3_ssp2$`TPOP_age_gr_69+`))*100, 2)
d6 <- round(abs((sum(proj_res3_ssp3$`TPOP_age_gr_69+`) - sum(proj_res3_ssp3$`age_gr_69+`)) / 1e6), 2)

d7 <- mape(proj_res3_ssp5$`TPOP_age_gr_69+`, proj_res3_ssp5$`age_gr_69+`)*100
d78 <- round(abs((sum(proj_res3_ssp5$`TPOP_age_gr_69+`) - sum(proj_res3_ssp5$`age_gr_69+`)) / sum(proj_res3_ssp2$`TPOP_age_gr_69+`))*100, 2)
d8 <- round(abs((sum(proj_res3_ssp5$`TPOP_age_gr_69+`) - sum(proj_res3_ssp5$`age_gr_69+`)) / 1e6), 2)

d0_m <- data.frame(country="United States", source="Striessnig et al.", scenario=c("SSP1", "SSP2", "SSP3", "SSP5"), ref_scen=c("SSP1", "SSP2", "SSP3", "SSP5"), spat_unit ="Counties", sum_error=c(d2, d4, d6, d8), pct_error = c(d12, d34, d56, d78), mape=c(d1, d3, d5, d7))

###

d1 <- mape(eu_data_ssp1$`T_69+`, eu_data_ssp1$`age_gr_69+`)*100
d12 <- round(abs((sum(eu_data_ssp1$`T_69+`) - sum(eu_data_ssp1$`age_gr_69+`)) / sum(eu_data_ssp1$`T_69+`))*100, 2)
d2 <- round(abs((sum(eu_data_ssp1$`T_69+`) - sum(eu_data_ssp1$`age_gr_69+`)) / 1e6), 2)

d3 <- mape(eu_data_ssp2$`T_69+`, eu_data_ssp2$`age_gr_69+`)*100
d34 <- round(abs((sum(eu_data_ssp2$`T_69+`) - sum(eu_data_ssp2$`age_gr_69+`)) / sum(eu_data_ssp2$`T_69+`))*100, 2)
d4 <- round(abs((sum(eu_data_ssp2$`T_69+`) - sum(eu_data_ssp2$`age_gr_69+`)) / 1e6), 2)

d5 <- mape(eu_data_ssp3$`T_69+`, eu_data_ssp3$`age_gr_69+`)*100
d56 <- round(abs((sum(eu_data_ssp3$`T_69+`) - sum(eu_data_ssp3$`age_gr_69+`)) / sum(eu_data_ssp3$`T_69+`))*100, 2)
d6 <- round(abs((sum(eu_data_ssp3$`T_69+`) - sum(eu_data_ssp3$`age_gr_69+`)) / 1e6), 2)

d7 <- mape(eu_data_ssp5$`T_69+`, eu_data_ssp5$`age_gr_69+`)*100
d78 <- round(abs((sum(eu_data_ssp5$`T_69+`) - sum(eu_data_ssp5$`age_gr_69+`)) / sum(eu_data_ssp5$`T_69+`))*100, 2)
d8 <- round(abs((sum(eu_data_ssp5$`T_69+`) - sum(eu_data_ssp5$`age_gr_69+`)) / 1e6), 2)

d0b_m <- data.frame(country="European Union", source="Eurostat", scenario=c("SSP1", "SSP2", "SSP3", "SSP5"), ref_scen=c("Own modelling"), spat_unit ="NUTS3", sum_error=c(d2, d4, d6, d8), pct_error = c(d12, d34, d56, d78), mape=c(d1, d3, d5, d7))


###

d1 <- mape(usa2_age_gr_69_ssp1_2050$SSP1, usa2_age_gr_69_ssp1_2050$`age_gr_69+`)*100
d12 <- round(abs((sum(usa2_age_gr_69_ssp1_2050$SSP1) - sum(usa2_age_gr_69_ssp1_2050$`age_gr_69+`)) / sum(usa2_age_gr_69_ssp1_2050$SSP1))*100, 2)
d2 <- round(abs((sum(usa2_age_gr_69_ssp1_2050$SSP1) - sum(usa2_age_gr_69_ssp1_2050$`age_gr_69+`)) / 1e6), 2)

d3 <- mape(usa2_age_gr_69_ssp2_2050$SSP2, usa2_age_gr_69_ssp2_2050$`age_gr_69+`)*100
d34 <- round(abs((sum(usa2_age_gr_69_ssp2_2050$SSP2) - sum(usa2_age_gr_69_ssp2_2050$`age_gr_69+`)) / sum(usa2_age_gr_69_ssp2_2050$SSP2))*100, 2)
d4 <- round(abs((sum(usa2_age_gr_69_ssp2_2050$SSP2) - sum(usa2_age_gr_69_ssp2_2050$`age_gr_69+`)) / 1e6), 2)

d5 <- mape(usa2_age_gr_69_ssp3_2050$SSP3, usa2_age_gr_69_ssp3_2050$`age_gr_69+`)*100
d56 <- round(abs((sum(usa2_age_gr_69_ssp3_2050$SSP3) - sum(usa2_age_gr_69_ssp3_2050$`age_gr_69+`)) / sum(usa2_age_gr_69_ssp3_2050$SSP3))*100, 2)
d6 <- round(abs((sum(usa2_age_gr_69_ssp3_2050$SSP3) - sum(usa2_age_gr_69_ssp3_2050$`age_gr_69+`)) / 1e6), 2)

d7 <- mape(usa2_age_gr_69_ssp5_2050$SSP5, usa2_age_gr_69_ssp5_2050$`age_gr_69+`)*100
d78 <- round(abs((sum(usa2_age_gr_69_ssp5_2050$SSP5) - sum(usa2_age_gr_69_ssp5_2050$`age_gr_69+`)) / sum(usa2_age_gr_69_ssp5_2050$SSP5))*100, 2)
d8 <- round(abs((sum(usa2_age_gr_69_ssp5_2050$SSP5) - sum(usa2_age_gr_69_ssp5_2050$`age_gr_69+`)) / 1e6), 2)

d1_m <- data.frame(country="United States", source="Hauer", scenario=c("SSP1", "SSP2", "SSP3", "SSP5"), ref_scen=c("SSP1", "SSP2", "SSP3", "SSP5"), spat_unit ="Counties", sum_error=c(d2, d4, d6, d8), pct_error = c(d12, d34, d56, d78), mape=c(d1, d3, d5, d7))

###

d1 <- mape(uk3_ssp1_2050$`age_gr_69+`, uk3_ssp1_2050$`age_gr_69+_mod`)*100
d12 <- round(abs((sum(uk3_ssp1_2050$`age_gr_69+`) - sum(uk3_ssp1_2050$`age_gr_69+_mod`)) / sum(uk3_ssp1_2050$`age_gr_69+`))*100, 2)
d2 <- round(abs((sum(uk3_ssp1_2050$`age_gr_69+`) - sum(uk3_ssp1_2050$`age_gr_69+_mod`)) / 1e6), 2)

d3 <- mape(uk3_ssp2_2050$`age_gr_69+`, uk3_ssp2_2050$`age_gr_69+_mod`)*100
d34 <- round(abs((sum(uk3_ssp2_2050$`age_gr_69+`) - sum(uk3_ssp2_2050$`age_gr_69+_mod`)) / sum(uk3_ssp2_2050$`age_gr_69+`))*100, 2)
d4 <- round(abs((sum(uk3_ssp2_2050$`age_gr_69+`) - sum(uk3_ssp2_2050$`age_gr_69+_mod`)) / 1e6), 2)

d5 <- mape(uk3_ssp3_2050$`age_gr_69+`, uk3_ssp3_2050$`age_gr_69+_mod`)*100
d56 <- round(abs((sum(uk3_ssp3_2050$`age_gr_69+`) - sum(uk3_ssp3_2050$`age_gr_69+_mod`)) / sum(uk3_ssp3_2050$`age_gr_69+`))*100, 2)
d6 <- round(abs((sum(uk3_ssp3_2050$`age_gr_69+`) - sum(uk3_ssp3_2050$`age_gr_69+_mod`)) / 1e6), 2)

d7 <- mape(uk3_ssp5_2050$`age_gr_69+`, uk3_ssp5_2050$`age_gr_69+_mod`)*100
d78 <- round(abs((sum(uk3_ssp5_2050$`age_gr_69+`) - sum(uk3_ssp5_2050$`age_gr_69+_mod`)) / sum(uk3_ssp5_2050$`age_gr_69+`))*100, 2)
d8 <- round(abs((sum(uk3_ssp5_2050$`age_gr_69+`) - sum(uk3_ssp5_2050$`age_gr_69+_mod`)) / 1e6), 2)

d2_m <- data.frame(country="United Kingdom", source="Nash", scenario=c("SSP1", "SSP2", "SSP3", "SSP5"), ref_scen=c("Own modelling"), spat_unit ="Local authorities", sum_error=c(d2, d4, d6, d8), pct_error = c(d12, d34, d56, d78), mape=c(d1, d3, d5, d7))

###

d1 <- mape(china3_ssp1_2050$`age_gr_69+`, china3_ssp1_2050$`age_gr_69+_mod`)*100
d12 <- round(abs((sum(china3_ssp1_2050$`age_gr_69+`) - sum(china3_ssp1_2050$`age_gr_69+_mod`)) / sum(china3_ssp1_2050$`age_gr_69+`))*100, 2)
d2 <- round(abs((sum(china3_ssp1_2050$`age_gr_69+`) - sum(china3_ssp1_2050$`age_gr_69+_mod`)) / 1e6), 2)

d3 <- mape(china3_ssp2_2050$`age_gr_69+`, china3_ssp2_2050$`age_gr_69+_mod`)*100
d34 <- round(abs((sum(china3_ssp2_2050$`age_gr_69+`) - sum(china3_ssp2_2050$`age_gr_69+_mod`)) / sum(china3_ssp2_2050$`age_gr_69+`))*100, 2)
d4 <- round(abs((sum(china3_ssp2_2050$`age_gr_69+`) - sum(china3_ssp2_2050$`age_gr_69+_mod`)) / 1e6), 2)

d5 <- mape(china3_ssp3_2050$`age_gr_69+`, china3_ssp3_2050$`age_gr_69+_mod`)*100
d56 <- round(abs((sum(china3_ssp3_2050$`age_gr_69+`) - sum(china3_ssp3_2050$`age_gr_69+_mod`)) / sum(china3_ssp3_2050$`age_gr_69+`))*100, 2)
d6 <- round(abs((sum(china3_ssp3_2050$`age_gr_69+`) - sum(china3_ssp3_2050$`age_gr_69+_mod`)) / 1e6), 2)


d7 <- mape(china3_ssp5_2050$`age_gr_69+`, china3_ssp5_2050$`age_gr_69+_mod`)*100
d78 <- round(abs((sum(china3_ssp5_2050$`age_gr_69+`) - sum(china3_ssp5_2050$`age_gr_69+_mod`)) / sum(china3_ssp5_2050$`age_gr_69+`))*100, 2)
d8 <- round(abs((sum(china3_ssp5_2050$`age_gr_69+`) - sum(china3_ssp5_2050$`age_gr_69+_mod`)) / 1e6), 2)

d3_m <- data.frame(country="China", source="Chen et al.", scenario=c("SSP1", "SSP2", "SSP3", "SSP5"), ref_scen=c("SSP1", "SSP2", "SSP3", "SSP5"), spat_unit ="Provinces", sum_error=c(d2, d4, d6, d8), pct_error = c(d12, d34, d56, d78), mape=c(d1, d3, d5, d7))

###

d1 <- mape(india$value, india$`age_gr_69+`)*100
d12 <- round(abs((sum(india$value) - sum(india$`age_gr_69+`)) / sum(india$value))*100, 2)
d2 <- round(abs((sum(india$value) - sum(india$`age_gr_69+`)) / 1e6), 2)

d4_m <- data.frame(country="India", source="KC et al.", scenario=c("SSP2"), spat_unit ="Regions", ref_scen="SSP2", sum_error=c(d2), pct_error = c(d12), mape=c(d1))

d_all <- bind_rows(d0_m, d1_m, d0b_m, d2_m, d3_m, d4_m)

###

colnames(d_all) <- c("Country", "Source", "Scenario", "Reference scenario", "Spat. unit", "Tot. err. (mil. people 69+)", "Perc. err. (%)", "Avg. local pctg. err. (%)")

stargazer(d_all, type="latex", out = "compar_projections_other_data_sources.tex", summary = F, title="A", digits=1, rownames = F)
