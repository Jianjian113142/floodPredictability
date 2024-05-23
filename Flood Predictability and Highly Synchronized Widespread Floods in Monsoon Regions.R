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


library(seasonal)
library(stlplus)
library(lubridate)
library(tidyverse)



# Significance test of the difference between the flood events 
# in the selected time period and the mean value in the monsoon period-----------------


time <- seq(as.Date("2002-04-01"),as.Date("2016-08-31"),by="day")

#Choose a leap year to cycle through every day. Here, choose 2004.
Year2004 <- seq(as.Date("2004-01-01"),as.Date("2004-12-31"),by="day")

#date
dateTws <- seq(as.Date('2002-04-01'), as.Date('2016-8-31'), by='day')


#Find the monsoon period, which is mainly from June to September in the Northern Hemisphere, and from December to March in the Southern Hemisphere.
loc2 <- month(dateTws) %in% c(6,7,8,9)
#south hemisphere
loc3 <-  month(dateTws) %in% c(12,1,2,3)


##average in Northern hemisphere monsoon period 
Mean_fld2 <- apply(floodevents_morethan1mm_95th_2002_2016[loc2,], 2, mean)
std_fld2 <- apply(floodevents_morethan1mm_95th_2002_2016[loc2,], 2, sd)


#average in Southern hemisphere monsoon period 
Mean_fld3 <-  apply(floodevents_morethan1mm_95th_2002_2016[loc3,], 2, mean)
std_fld3 <- apply(floodevents_morethan1mm_95th_2002_2016[loc3,], 2, sd)



# for NAmerM-------------------------------------------

#Northern hemisphere monsoon period

NorthComp <- function(MonsT,Mean_fld2,std_fld2,type,loc2){
  
  loc <- match(MonsT,dateTws)
  

  Mean_fld1 <-  apply(floodevents_morethan1mm_95th_2002_2016[loc,], 2, mean)
  std_fld1 <- apply(floodevents_morethan1mm_95th_2002_2016[loc,], 2, sd)
  
  n1 <- length(loc)
  n2 <- sum(loc2)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_fld1,std_fld1,Mean_fld2,std_fld2,n1,n2)
  
  Mons_Test[is.nan(Mons_Test)] <- NA
  
  deltafld <- Mean_fld1 - Mean_fld2
  
  ModelTIf <- rep(NA, length(deltafld))
  
  dx <- Mons_Test>CriticalV
  dx[is.na(dx)] <- FALSE #Set location in NA as False
  
  ModelTIf[dx] <- deltafld[dx]
  
  return(ModelTIf)
}



#Southern hemisphere monsoon period 
SouthComp <- function(MonsT,Mean_fld3,std_fld3,type,loc3){

  loc <- match(MonsT,dateTws)

  Mean_fld1 <-  apply(floodevents_morethan1mm_95th_2002_2016[loc,], 2, mean)
  std_fld1 <- apply(floodevents_morethan1mm_95th_2002_2016[loc,], 2, sd)
  
  n1 <- length(loc)
  n2 <- sum(loc3)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_fld1,std_fld1,Mean_fld3,std_fld3,n1,n2)
  Mons_Test[is.nan(Mons_Test)] <- NA
  
  deltafld <- Mean_fld1 - Mean_fld3 
  
  ModelTIf <- rep(NA, length(deltafld))
  
  dx <- Mons_Test>CriticalV
  dx[is.na(dx)] <- FALSE 
  
  ModelTIf[dx] <- deltafld[dx]
  
  return(ModelTIf)
  
}


# Function: Calculate whether each grid point is significant based on degrees of freedom and threshold ------------------------------------------------------

TTest <- function(Mean1,Std1,Mean2,Std2,n1,n2){
  
  Ttest <- (Mean1-Mean2)/(sqrt((n1*Std1^2+n2*Std2^2)/(n1+n2-2))*sqrt(1/n1+1/n2))
  
  return(Ttest)
}


# Function: calculation of critical value ---------------------------------------------------------------
critical_value <- function(n1,n2){

  df <- n1 + n2 - 2  # 
  alpha <- 0.05  # 
  
  critical_value <- qt(1 - alpha, df)
  return(critical_value)
}





# get model data ------------------------------------------------------------------

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


#Find the position of the selected grid point in shpmodelSecNEW, and then retain the corresponding grid point according to the selected time
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

#read monsoon regions boundary
Mons <-  st_read("NEW_MonsoonRegions.shp")

# NAmerM ------------------------------------------------------------------

xmin <- -134
ymin <- 13
xmax <- -77
ymax <- 44
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系


inter <- st_intersection(shpmodelSecNEW,rectangle_sf)
plot(inter)



indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)


fld_NAmerM <- NorthComp(NAmerM,Mean_fld2,std_fld2,"NAmerM",loc2)
fldAnomDiff <- fld_NAmerM[indx]


# plot ---

drawFlds_NAmerM <- function(locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  # Jbrs[2] <- 1
  
  # plot(locat)
  coords <- st_coordinates(locat)
  locat$lon <- coords[, "X"]  # 提取经度
  locat$lat <- coords[, "Y"]  # 提取纬度
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
    geom_sf(data = countryGlobal,fill=NA,color='Black',linewidth = 0.1,alpha=0.9)+
    geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
    geom_point(data = locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
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
      # legend.key.height = unit(2, "cm"),  # 调整图例键高度
      legend.key.width = unit(0.3, "cm"),
      legend.text.align = 1
    )
  
  # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
  output <- paste(type,".png",sep = "")
  ggsave(output,
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}



NAmerM_locat <- getFldLoc(shpmodelSecNEW,"NAmerM")
drawFlds_NAmerM(NAmerM_locat,"NAmerM",fldAnomDiff,indxLoc,shpmodelSecNEW)



# EqAmerM -----------------------------------------------------------------

xmin <- -110
ymin <- -10
xmax <- -40
ymax <- 30
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系


inter <- st_intersection(shpmodelSecNEW,rectangle_sf)


indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)


fld_EqAmerM <- NorthComp(EqAmerM,Mean_fld2,std_fld2,"EqAmerM",loc2)
fldAnomDiff <- fld_EqAmerM[indx]



drawFlds_EqAmerM <- function(locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(locat)
  coords <- st_coordinates(locat)
  locat$lon <- coords[, "X"]  # 提取经度
  locat$lat <- coords[, "Y"]  # 提取纬度
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
    geom_sf(data = countryGlobal,fill=NA,color='Black',linewidth = 0.1,alpha=0.9)+
    geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
    geom_point(data = locat,aes(x = lon, y = lat),colour="black",alpha=0.5,size = 0.0001,stroke = 0.5)+
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
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}


EqAmerM_locat <- getFldLoc(shpmodelSecNEW,"EqAmerM")

drawFlds_EqAmerM(EqAmerM_locat,"EqAmerM",fldAnomDiff,indxLoc,shpmodelSecNEW)



# WAfriM ------------------------------------------------------------------
xmin <- -30
ymin <- -16
xmax <- 65
ymax <- 28
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系



inter <- st_intersection(shpmodelSecNEW,rectangle_sf)
plot(inter)


indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)


fld_WAfriM <- NorthComp(WAfriM,Mean_fld2,std_fld2,"WAfriM",loc2)
fldAnomDiff <- fld_WAfriM[indx]



drawFlds_WAfriM <- function(AusMCM_locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(AusMCM_locat)
  coords <- st_coordinates(AusMCM_locat)
  AusMCM_locat$lon <- coords[, "X"]  # 
  AusMCM_locat$lat <- coords[, "Y"]  # 
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
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
      # legend.position = "bottom"
      legend.key.width = unit(0.3, "cm"),
      legend.text.align = 1
    )
  
  # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
  output <- paste(type,".png",sep = "")
  ggsave(output,
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}



WAfriM_locat <- getFldLoc(shpmodelSecNEW,"WAfriM")
drawFlds_WAfriM(WAfriM_locat,"WAfriM",fldAnomDiff,indxLoc,shpmodelSecNEW)




# SAsiaM ------------------------------------------------------------------

xmin <- 60
ymin <- 0
xmax <- 135
ymax <- 40
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系

inter <- st_intersection(shpmodelSecNEW,rectangle_sf)
plot(inter)

indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)

fld_SAsiaM <- NorthComp(SAsiaM,Mean_fld2,std_fld2,"SAsiaM",loc2)
fldAnomDiff <- fld_SAsiaM[indx]


drawFlds_SAsiaM <- function(AusMCM_locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(AusMCM_locat)
  coords <- st_coordinates(AusMCM_locat)
  AusMCM_locat$lon <- coords[, "X"]  # 提取经度
  AusMCM_locat$lat <- coords[, "Y"]  # 提取纬度
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
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
      # legend.position = "bottom"
      legend.key.width = unit(0.3, "cm"),
      legend.text.align = 1
    )
  
  # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
  output <- paste(type,".png",sep = "")
  ggsave(output,
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}


SAsiaM_locat <- getFldLoc(shpmodelSecNEW,"SAsiaM")
drawFlds_SAsiaM(SAsiaM_locat,"SAsiaM",fldAnomDiff,indxLoc,shpmodelSecNEW)




# EAsiaM ------------------------------------------------------------------
xmin <- 90
ymin <- 14
xmax <- 145
ymax <- 50
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系



inter <- st_intersection(shpmodelSecNEW,rectangle_sf)
plot(inter)



indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)



fld_EAsiaM <- NorthComp(EAsiaM,Mean_fld2,std_fld2,"EAsiaM",loc2)
fldAnomDiff <- fld_EAsiaM[indx]




drawFlds_EAsiaM <- function(AusMCM_locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(AusMCM_locat)
  coords <- st_coordinates(AusMCM_locat)
  AusMCM_locat$lon <- coords[, "X"]  # 
  AusMCM_locat$lat <- coords[, "Y"]  # 
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
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
      # legend.position = "bottom"
      legend.key.width = unit(0.3, "cm"),
      legend.text.align = 1
    )
  
  # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
  output <- paste(type,".png",sep = "")
  ggsave(output,
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}



EAsiaM_locat <- getFldLoc(shpmodelSecNEW,"EAsiaM")
drawFlds_EAsiaM(EAsiaM_locat,"EAsiaM",fldAnomDiff,indxLoc,shpmodelSecNEW)





# SAmerM  ------------------------------------------------------------
xmin <- -100
ymin <- -40
xmax <- -20
ymax <- 9
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系


inter <- st_intersection(shpmodelSecNEW,rectangle_sf)
plot(inter)


indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)


fld_SAmerM <- SouthComp(SAmerM,Mean_fld3,std_fld3,"SAmerM",loc3)
fldAnomDiff <- fld_SAmerM[indx]



drawFlds_SAmerM <- function(AusMCM_locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(AusMCM_locat)
  coords <- st_coordinates(AusMCM_locat)
  AusMCM_locat$lon <- coords[, "X"]  # 
  AusMCM_locat$lat <- coords[, "Y"]  #
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
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
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}



SAmerM_locat <- getFldLoc(shpmodelSecNEW,"SAmerM")
drawFlds_SAmerM(SAmerM_locat,"SAmerM",fldAnomDiff,indxLoc,shpmodelSecNEW)




# SAfriM ------------------------------------------------------------------

xmin <- 0
ymin <- -38
xmax <- 55
ymax <- 3
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系



inter <- st_intersection(shpmodelSecNEW,rectangle_sf)
plot(inter)


indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)



fld_SAfriM <- SouthComp(SAfriM,Mean_fld3,std_fld3,"SAfriM",loc3)
fldAnomDiff <- fld_SAfriM[indx]



drawFlds_SAfriM <- function(AusMCM_locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(AusMCM_locat)
  coords <- st_coordinates(AusMCM_locat)
  AusMCM_locat$lon <- coords[, "X"]  
  AusMCM_locat$lat <- coords[, "Y"] 
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
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
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}



SAfriM_locat <- getFldLoc(shpmodelSecNEW,"SAfriM")
drawFlds_SAfriM(SAfriM_locat,"SAfriM",fldAnomDiff,indxLoc,shpmodelSecNEW)



# Display the corresponding monsoon area-precipitation anomaly-flood distribution ------AusMCM------------------------------------------------
# Define the coordinates of the four corner points of the rectangle
xmin <- 90
ymin <- -30
xmax <- 170
ymax <- 20
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 



inter <- st_intersection(shpmodelSecNEW,rectangle_sf)

plot(inter)



indx <- match(inter$id,as.numeric(id13882))

indxLoc <- match(inter$id,shpmodelSecNEW$id)


fld_AusMCM <- SouthComp(AusMCM,Mean_fld3,std_fld3,"AusMCM",loc3)
fldAnomDiff <- fld_AusMCM[indx]


drawFlds_AusMCM <- function(AusMCM_locat,type,fldAnomDiff,indxLoc,shpmodelSecNEW){
  
  scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
    binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
  }
  
  
  shpmodelSecNEW$fldEnt <- NA
  shpmodelSecNEW$fldEnt[indxLoc] <- fldAnomDiff
  shpmodelSecNEWnew <- shpmodelSecNEW[!is.na(shpmodelSecNEW$fldEnt),]
  
  
  twscols <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  
  Jbrs <- round(getJenksBreaks(shpmodelSecNEWnew$fldEnt, 7, subset = NULL),digits = 4)
  
  # plot(AusMCM_locat)
  coords <- st_coordinates(AusMCM_locat)
  AusMCM_locat$lon <- coords[, "X"]  # 提取经度
  AusMCM_locat$lat <- coords[, "Y"]  # 提取纬度
  
  pp <-
    ggplot()+
    geom_sf(data = shpmodelSecNEWnew,aes(fill=fldEnt),color = NA)+
    # geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
    scale_fill_distiller_custom(twscols,name = "Number of floods",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
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
      # legend.position = "bottom"
      legend.key.width = unit(0.3, "cm"),
      legend.text.align = 1
    )
  
  # pp <- pp + guides(fill = guide_colorbar(title.position = "right"))
  output <- paste(type,".png",sep = "")
  ggsave(output,
         plot=pp, width = 4, height = 2.3, dpi = 300,units="in")
  
}

AusMCM_locat <- getFldLoc(shpmodelSecNEW,"AusMCM")
drawFlds_AusMCM(AusMCM_locat,"AusMCM",fldAnomDiff,indxLoc,shpmodelSecNEW)
































