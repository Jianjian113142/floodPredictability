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
#modified DET
sf_use_s2(FALSE)




getHDetSHP <- function(){
  input <- "MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv"
  
  DETall <- read.csv(file = input,header = F)
  DETall$V1[is.nan(DETall$V1)] <- NA
  
  DETall$bina <- '0'
  DETall$bina[!is.na(DETall$V1)] <- '1'
  
  shpmodelSec$Value <- NA
  shpmodelSec$bina <- NA
  shpmodelSec$Value[!nodata] <- as.numeric(DETall$V1)
  shpmodelSec$bina[!nodata] <- DETall$bina
  
  shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
  shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$bina),]
  return(shpmodelSecNEW)
  
}



shpmodelSecNEW <- getHDetSHP()

# Draw the spatial distribution of flood events according to the selected time--------------------------------------------------------


getFldLoc <- function(Model,type){
  
  MonsSec <- Mons[Mons$Monsoon==type,]
  plot(MonsSec)
  
  inter <- st_intersection(Model,MonsSec)
  
  inter1 <- inter[inter$bina=="1",]

  # 计算每个 polygon 的中心点
  centroids <- st_centroid(inter1)
  
  # plot(centroids)
  
  return(centroids)
  
}



Mons <- st_read("NEW_MonsoonRegions.shp")

#------AusMCM------------------------------------------------

type <- "AusMCM"
xmin <- 90
ymin <- -30
xmax <- 170
ymax <- 20
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系

rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)


rasterAusMCM <- brick("AusMCM降水均值相对季风时段差值.tif")

Cuttif <- crop(rasterAusMCM,rectangle_sp)

plot(Cuttif)



Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 
Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]


AusMCM_locat <- getFldLoc(shpmodelSecNEW,"AusMCM")

scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}

# twscols <- paletteer_d("ggsci::green_material")
twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')


Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)

coords <- st_coordinates(AusMCM_locat)
AusMCM_locat$lon <- coords[, "X"]  
AusMCM_locat$lat <- coords[, "Y"]  

pp <-
  ggplot()+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  scale_fill_distiller_custom(twscols,name = "Pre",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  geom_sf(data = countryGlobal,fill=NA,color='Black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  # geom_sf(data = AusMCM_locat,fill=NA,color= "black",alpha=0.9,size = 0.0001)+
  geom_point(data = AusMCM_locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
  coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
  xlab("")+ylab("")+
  theme(
    text=element_text(family="Helvetica"),
    # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
    # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
    legend.title = element_blank(),
    # legend.position = "bottom"
    legend.key.width = unit(0.3, "cm"),
    legend.text.align = 1
  )

# pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
output <- paste(type,".png",sep = "")
ggsave(output,
       plot=pp, width = 4, height = 2.3,dpi=300,units = "in")



# EqAmerM -----------------------------------------------------------------

xmin <- -110
ymin <- -10
xmax <- -40
ymax <- 30
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系

rectangle_sp <- as(rectangle_sf, "Spatial")


# plot(rectangle_sf)

drawPRE_all_EqAmerM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')

  
    rasterAusMCM <- brick("EqAmerM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    

    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    coords <- st_coordinates(AusMCM_locat)
    AusMCM_locat$lon <- coords[, "X"]  
    AusMCM_locat$lat <- coords[, "Y"]  
    
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "precipitation Anomaly (mm)",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='Black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_point(data = AusMCM_locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
    output <- paste(type,".png",sep = "")
    
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
  
}


EqAmerM_locat <- getFldLoc(shpmodelSecNEW,"EqAmerM")

plot(EqAmerM_locat)

drawPRE_all_EqAmerM(EqAmerM_locat,"EqAmerM")




# SAsiaM ------------------------------------------------------------------

xmin <- 48
ymin <- -6
xmax <- 153
ymax <- 45
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  

rectangle_sp <- as(rectangle_sf, "Spatial")


drawPRE_all_SAsiaM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  # twscols <- paletteer_d("ggsci::green_material")
  # twscols <- rev(c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4'))
  
  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')
  
  
    rasterAusMCM <- brick("SAsiaM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    
    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    coords <- st_coordinates(AusMCM_locat)
    AusMCM_locat$lon <- coords[, "X"] 
    AusMCM_locat$lat <- coords[, "Y"] 
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "Pre",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='Black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_point(data = AusMCM_locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    output <- paste(type,".png",sep = "")
    
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
    
    

  
}


SAsiaM_locat <- getFldLoc(shpmodelSecNEW,"SAsiaM")

plot(SAsiaM_locat)

drawPRE_all_SAsiaM(SAsiaM_locat,"SAsiaM")




# SAmerM  ------------------------------------------------------------
xmin <- -100
ymin <- -40
xmax <- -20
ymax <- 9
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  

rectangle_sp <- as(rectangle_sf, "Spatial")


drawPRE_all_SAmerM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')
  
  Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
  

    
    rasterAusMCM <- brick("SAmerM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    
    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    coords <- st_coordinates(AusMCM_locat)
    AusMCM_locat$lon <- coords[, "X"]  
    AusMCM_locat$lat <- coords[, "Y"]  
    
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "precipitation Anomaly (mm)",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_point(data = AusMCM_locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
    output <- paste(type,".png",sep = "")
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
    
    
  
}


SAmerM_locat <- getFldLoc(shpmodelSecNEW,"SAmerM")

plot(SAmerM_locat)

drawPRE_all_SAmerM(SAmerM_locat,"SAmerM")






# SAfriM ------------------------------------------------------------------
xmin <- 0
ymin <- -38
xmax <- 55
ymax <- 3
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))

rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  

rectangle_sp <- as(rectangle_sf, "Spatial")


drawPRE_all_SAfriM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')
  
    
    rasterAusMCM <- brick("SAfriM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    
    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    coords <- st_coordinates(AusMCM_locat)
    AusMCM_locat$lon <- coords[, "X"] 
    AusMCM_locat$lat <- coords[, "Y"]  
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "precipitation Anomaly (mm)",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_point(data = AusMCM_locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),

        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
    output <- paste(type,".png",sep = "")
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
    
  
}


SAfriM_locat <- getFldLoc(shpmodelSecNEW,"SAfriM")

plot(SAfriM_locat)

drawPRE_all_SAfriM(SAfriM_locat,"SAfriM")




# WAfriM ------------------------------------------------------------------

xmin <- -30
ymin <- -16
xmax <- 65
ymax <- 28
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  

rectangle_sp <- as(rectangle_sf, "Spatial")


drawPRE_all_WAfriM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  

  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')




    rasterAusMCM <- brick("WAfriM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    

    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE) 
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    
    
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "precipitation Anomaly (mm)",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_sf(data = AusMCM_locat,fill=NA,color= "black",alpha=0.9,size = 0.0001)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
    output <- paste(type,".png",sep = "")
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
    
    
  
}


WAfriM_locat <- getFldLoc(shpmodelSecNEW,"WAfriM")

plot(WAfriM_locat)

drawPRE_all_WAfriM(WAfriM_locat,"WAfriM")



# EAsiaM ------------------------------------------------------------------

xmin <- 90
ymin <- 14
xmax <- 145
ymax <- 50
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  

rectangle_sp <- as(rectangle_sf, "Spatial")


drawPRE_all_EAsiaM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  # twscols <- paletteer_d("ggsci::green_material")
  # twscols <- rev(c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4'))
  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')

    rasterAusMCM <- brick("~/Desktop/洪水日数95th-morethan1mm/GPM_geotif_result/EAsiaM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    
    # 将 raster 数据转换为 data frame
    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    
    coords <- st_coordinates(AusMCM_locat)
    AusMCM_locat$lon <- coords[, "X"]  
    AusMCM_locat$lat <- coords[, "Y"]  
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "precipitation Anomaly (mm)",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_point(data = AusMCM_locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
    output <- paste(type,".png",sep = "")
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
    
}


EAsiaM_locat <- getFldLoc(shpmodelSecNEW,"EAsiaM")

plot(EAsiaM_locat)

drawPRE_all_EAsiaM(EAsiaM_locat,"EAsiaM")


# NAmerM ------------------------------------------------------------------

xmin <- -134
ymin <- 13
xmax <- -77
ymax <- 44
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  

rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sp)


drawPRE_all_NAmerM <- function(AusMCM_locat,type){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  twscols <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494')

    rasterAusMCM <- brick("NAmerM降水均值相对季风时段差值.tif")
    
    Cuttif <- crop(rasterAusMCM,rectangle_sp)
    

    Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  
    Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]
    Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
    
    
    pp <-
      ggplot()+
      # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
      geom_tile(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
      scale_fill_distiller_custom(twscols,name = "precipitation Anomaly (mm)",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
      geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
      geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
      geom_sf(data = AusMCM_locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
      coord_sf(ylim = c(ymin,ymax),xlim = c(xmin,xmax),expand = F)+
      xlab("")+ylab("")+
      theme(
        text=element_text(family="Helvetica"),
        # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
        axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
        axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
        axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
        # legend.title = element_text(family = "Helvetica",size = 10,face = "bold"),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.text.align = 1
      )
    
    # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
    output <- paste(type,".png",sep = "")
    ggsave(output,
           plot=pp, width = 4, height = 2.3,dpi=300,units = "in")
    
    
}


NAmerM_locat <- getFldLoc(shpmodelSecNEW,"NAmerM")

plot(NAmerM_locat)

drawPRE_all_NAmerM(NAmerM_locat,"NAmerM")



















