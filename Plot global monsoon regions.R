

library(ncdf4)
library(tidyverse)
library(iterators)
library(foreach)
library(sf)
library(doParallel)
library(raster)
library(ggplot2)
library(parallel)
theme_set(theme_bw())
library(zoo)
library(tools)
library(tidyr)
library('ggedit')
library(ggplot2)
library(paletteer)
# library(XploreR)
library(BAMMtools)
library(paletteer)
# library("rjson")
#modified DET
sf_use_s2(FALSE)
library(jsonlite)



#get monsoon regions
WGI <- st_read("NEW_MonsoonRegions.shp")

WGI$Cor <- NA

#计算矢量面的中心坐标
centroids <- st_centroid(WGI)
coordinates <- st_coordinates(centroids)
WGI$lon <- coordinates[, "X"]
WGI$lat <- coordinates[, "Y"]


# #get global country bundary
# countryGlobal <- st_read("countryGlobal.shp")


color <- c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd')


P <- 
  ggplot()+
  # geom_sf(data = WGI,aes(fill=Cor),color = NA)+
  geom_sf(data = countryGlobal,fill="#e0e0e0",color= NA,size=0.01,alpha=0.8)+
  geom_sf(data = WGI,fill=color,color='#4575b4',size=3)+
  geom_text(data = WGI,
            aes(x = lon, y = lat, label = Monsoon),color = "black",
            size = 3) +
  coord_sf(ylim = c(-60,60),expand = F)+
  # guides(guide = "none")+
  theme(
    text=element_text(family="Helvetica"),
    # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold"),
    legend.title = element_text(family = "Helvetica",size = 12,face = "bold")
    # legend.title = element_blank()
  )



ggsave("Monsoon-reference-regions.png",plot=P, width = 8, height = 3, dpi = 300,units="in")









































































