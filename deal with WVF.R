# code GB2312
library(tidyr)
library(dplyr)
library(readxl)
library(stats)
library(paletteer)
library(ggplot2)
library(sf)
theme_set(theme_bw())
library(raster)
library(zoo)
library(extrafont)
library(ncdf4)
library(seasonal)
library(lubridate)
library(pracma)
library(tidyr)
library(dplyr)
library(readxl)
library(stats)
library(paletteer)
library(ggplot2)
library(sf)
theme_set(theme_bw())
library(raster)
library(zoo)
library(extrafont)
library(lubridate)
library(foreach)
library(doParallel)
library(sf)
library(parallel)
library(fasstr)
library(BAMMtools)
library(ggarchery)
library(tools)
# loadfonts(dev="win")

sf_use_s2(FALSE)

Time <- seq(from = as.Date("2002-01-01"), to = as.Date("2017-12-31"), by = 'day')

# read v-wind --------------------------------------------------------------

vfile <- list_files_with_exts("vwind","nc")


# initial a Var.
stacked_V <- NULL

for (i in 1:length(vfile)){
  
    v <- nc_open(vfile[i])
    vwnd <- ncvar_get(v,"vwnd")
    # vlevel <- ncvar_get(v,"level")
    # t <- ncvar_get(v,"time")
    vwnd850 <- vwnd[,,3,]
    
    if (is.null(stacked_V)) {
      stacked_V <- vwnd850
    } else {
    
      stacked_V <- array(c(stacked_V, vwnd850), dim = c(dim(stacked_V)[1], dim(stacked_V)[2],dim(vwnd850)[3] + dim(stacked_V)[3]))
      
    }
    
    
    
}




#read u-wind -------------------------------------------------------------

ufile <-  list_files_with_exts("uwind","nc")


# 初始化一个空变量
stacked_U <- NULL

for (i in 1:length(ufile)){
  
  u <- nc_open(ufile[i])
  uwnd <- ncvar_get(u,"uwnd")
  uwnd850 <- uwnd[,,3,]
  
  if (is.null(stacked_U)) {
    stacked_U <- uwnd850
  } else {
    
    stacked_U <- array(c(stacked_U, uwnd850), dim = c(dim(stacked_U)[1], dim(stacked_U)[2],dim(uwnd850)[3] + dim(stacked_U)[3]))
    
  }
  
}


# Read specific humidity data ----------------------------------------------------------------

#read data function

readnc <- function(folder){
  
  
  ufile <-  list_files_with_exts(folder,"nc")
  
  # initial Var.
  stacked_U <- NULL
  
  for (i in 1:length(ufile)){
    
    u <- nc_open(ufile[i])
    uwnd <- ncvar_get(u,"shum")  
    uwnd850 <- uwnd[,,3,]
    
    if (is.null(stacked_U)) {
      stacked_U <- uwnd850
    } else {
      
      stacked_U <- array(c(stacked_U, uwnd850), dim = c(dim(stacked_U)[1], dim(stacked_U)[2],dim(uwnd850)[3] + dim(stacked_U)[3]))
      
    }
    
  }
  
  return(stacked_U)
  
}


stacked_shum <- readnc("shum")


# Then calculate the magnitude of the water vapor flux--------------------------------------------------------------

Q <- sqrt((stacked_U*stacked_shum)^2+(stacked_V*stacked_shum)^2)

Qsec <- Q[,,Time>=as.Date("2002-04-01")&Time<=as.Date("2016-08-31")]
dim(Qsec)

# The specified time average of the water vapor flux relative to the average during the monsoon period. --------------------------------------------------


#time length
dateTws <- seq(as.Date('2002-04-01'), as.Date('2016-8-31'), by='day')
# 



# begin ------------------------------------------------------------
NorthComp(NAmerM,"NAmerM",Qsec)
NorthComp(WAfriM,"WAfriM",Qsec)
NorthComp(SAsiaM,"SAsiaM",Qsec)
NorthComp(EAsiaM,"EAsiaM",Qsec)
NorthComp(EqAmerM,"EqAmerM",Qsec)

test <- SouthComp(SAmerM,"SAmerM",Qsec)
SouthComp(SAfriM,"SAfriM",Qsec)
SouthComp(AusMCM,"AusMCM",Qsec)

# hist(test)

#for NAmerM-------------------------------------------
#For the Northern Hemisphere, 
NorthComp <- function(MonsT,type,Qsec){
  #monsoon periods
  loc2 <- month(dateTws) %in% c(6,7,8,9)
  
  ##
  MSD2 <- Mean_STD(Qsec,loc2)
  Mean_pre2 <- MSD2[[1]]
  std_pre2 <- MSD2[[2]]

  loc <- match(MonsT,dateTws)
  
  #mean and std
  MSD <- Mean_STD(Qsec,loc)
  Mean_pre1 <- MSD[[1]]
  std_pre1 <- MSD[[2]]
  
  n1 <- length(loc)
  n2 <- sum(loc2)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_pre1,std_pre1,Mean_pre2,std_pre2,n1,n2)
  
  deltaPre <- Mean_pre1 - Mean_pre2
  
  ModelTIf <- array(NA, dim = dim(deltaPre))
  
  ModelTIf[Mons_Test>CriticalV] <- deltaPre[Mons_Test>CriticalV]
  
  print(dim(ModelTIf))
  
  TIF <- Mat2Tif(ModelTIf)
  plot(TIF)
  writeRaster(x = TIF,filename = paste(type,"WVF_composite_analysis_0.05.tif",sep=""))
  
}



#for the Southern Hemisphere

SouthComp <- function(MonsT,type,Qsec){
  
  #
  loc3 <-  month(dateTws) %in% c(12,1,2,3)
  #
  MSD3 <- Mean_STD(Qsec,loc3)
  Mean_pre3 <- MSD3[[1]]
  std_pre3 <- MSD3[[2]]
#index
  loc <- match(MonsT,dateTws)
  
  #mean and std
  MSD <- Mean_STD(Qsec,loc)
  Mean_pre1 <- MSD[[1]]
  std_pre1 <- MSD[[2]]
  
  n1 <- length(loc)
  n2 <- sum(loc3)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_pre1,std_pre1,Mean_pre3,std_pre3,n1,n2)
  
  deltaPre <- Mean_pre1 - Mean_pre3 #均值差异，是我们需要保留的地方
  
  ModelTIf <- array(NA, dim = dim(deltaPre))
  
  ModelTIf[Mons_Test>CriticalV] <- deltaPre[Mons_Test>CriticalV]
  
  print(dim(ModelTIf))
  
  TIF <- Mat2Tif(ModelTIf)
  plot(TIF)
  writeRaster(x = TIF,filename = paste(type,"WVF_composite_analysis_0.05.tif",sep=""))
  return(deltaPre)
}



#function _significant analysis using T-test-------------------

TTest <- function(Mean1,Std1,Mean2,Std2,n1,n2){
  
  Ttest <- (Mean1-Mean2)/(sqrt((n1*Std1^2+n2*Std2^2)/(n1+n2-2))*sqrt(1/n1+1/n2))
  
  return(Ttest)
}




# function critical value ---------------------------------------------------------------
critical_value <- function(n1,n2){
  # 计算右尾临界值
  df <- n1 + n2 - 2  # 自由度
  alpha <- 0.05  # 显著性水平
  
  critical_value <- qt(1 - alpha, df)
  return(critical_value)
}



# mean and std -----------------------------------------------------------------

Mean_STD <- function(Matrx,loc){
  
  dat <- Matrx[,,loc]
  
  mean_values <- apply(dat,c(1, 2),mean)
  
  std_dev <- apply(dat, c(1, 2), sd)
  
  return(list(mean_values,std_dev)) #返回均值
  
}



# convert matrix to geotiff---------------------------------------------------


Mat2Tif <- function(Arr){
  
  QfinalArr <- t(Arr)
  
  #将0-360变成-180~180
  QfinalArrMove <- array(data = NA,dim = c(73,144))
  QfinalArrMove[,1:72] <- QfinalArr[,73:144]
  QfinalArrMove[,73:144] <- QfinalArr[,1:72]
  
  
  r <- raster(QfinalArrMove, xmn=-180, xmx=180, ymn=-90, ymx=90, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # r <- flip(r, direction='y')
  return(r)
  
}



# significant direction for WVF --------------------------------------------------------------

#index for uq and vq
stacked_USec <- stacked_U[,,Time>=as.Date("2002-04-01")&Time<=as.Date("2016-08-31")]
stacked_VSec <- stacked_V[,,Time>=as.Date("2002-04-01")&Time<=as.Date("2016-08-31")]


# Significance-tested grid index for northern and southern hemispheres-------------------------------------------------------------------------

startFNorth <- function(MonsT,type,stacked){
  
  #
  dateTws <- seq(as.Date('2002-04-01'), as.Date('2016-8-31'), by='day')
  
  #Find the monsoon period, which is mainly from June to September in the Northern Hemisphere, and from December to March in the Southern Hemisphere.
  loc2 <- month(dateTws) %in% c(6,7,8,9)

  #
  MSD2 <- Mean_STD(stacked,loc2)
  Mean_pre2 <- MSD2[[1]]
  std_pre2 <- MSD2[[2]]

  idx <- NorthCompIndx(MonsT,Mean_pre2,std_pre2,type,stacked,loc2)
  
  return(idx)
  
  
}



startFSouth <- function(MonsT,type,stacked){
  
  #
  dateTws <- seq(as.Date('2002-04-01'), as.Date('2016-8-31'), by='day')
  
  #Find the monsoon period, which is mainly from June to September in the Northern Hemisphere, and from December to March in the Southern Hemisphere.
  loc3 <-  month(dateTws) %in% c(12,1,2,3)
  
  #
  MSD3 <- Mean_STD(stacked,loc3)
  Mean_pre3 <- MSD3[[1]]
  std_pre3 <- MSD3[[2]]
  
  
  idx <- SouthCompIndx(MonsT,Mean_pre3,std_pre3,type,stacked,loc3)

  return(idx)

}





#Index for Northern Hemisphere difference valid values-----------

NorthCompIndx <- function(MonsT,Mean_pre2,std_pre2,type,dat,loc2){
  
  loc <- match(MonsT,dateTws)
  
  #均值和标准差
  MSD <- Mean_STD(dat,loc)
  Mean_pre1 <- MSD[[1]]
  std_pre1 <- MSD[[2]]
  
  n1 <- length(loc)
  n2 <- sum(loc2)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_pre1,std_pre1,Mean_pre2,std_pre2,n1,n2)
  
  deltaPre <- Mean_pre1 - Mean_pre2
  
  # ModelTIf <- array(NA, dim = dim(deltaPre))
  
  # ModelTIf[Mons_Test>CriticalV] <- deltaPre[Mons_Test>CriticalV]
  IDX <- Mons_Test>CriticalV
  
  return(list(IDX,deltaPre)) #Returns the specific index value and mean difference
}



#Index for Southern Hemisphere difference valid values

SouthCompIndx <- function(MonsT,Mean_pre3,std_pre3,type,dat,loc3){

  
  loc <- match(MonsT,dateTws)
  
 
  MSD <- Mean_STD(dat,loc)
  Mean_pre1 <- MSD[[1]]
  std_pre1 <- MSD[[2]]
  
  n1 <- length(loc)
  n2 <- sum(loc3)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_pre1,std_pre1,Mean_pre3,std_pre3,n1,n2)
  
  deltaPre <- Mean_pre1 - Mean_pre3 #
  
  # ModelTIf <- array(NA, dim = dim(deltaPre))
  
  # ModelTIf[Mons_Test>CriticalV] <- deltaPre[Mons_Test>CriticalV]
  IDX <- Mons_Test>CriticalV
  
  return(list(IDX,deltaPre)) #Returns the specific index value and mean difference
  
}



# function compute the center of high-DET------------------------------------------------------------

getFldLoc <- function(Model,type){
  
  MonsSec <- Mons[Mons$Monsoon==type,]
  plot(MonsSec)
  
  inter <- st_intersection(Model,MonsSec)
  
  inter1 <- inter[inter$bina=="1",]
  
  # get every center of polygon 
  centroids <- st_centroid(inter1)
  
  # plot(centroids)
  
  return(centroids)
  
}


# get model and DET 
getHDetSHP <- function(){
  input <- "MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv"
  
  DETall <- read.csv(file = input,header = F)
  DETall$V1[is.nan(DETall$V1)] <- NA
  
  
  #
  DETall$bina <- '0'
  DETall$bina[!is.na(DETall$V1)] <- '1'
  
  shpmodelSec$Value <- NA
  shpmodelSec$bina <- NA
  shpmodelSec$Value[!nodata] <- as.numeric(DETall$V1)
  shpmodelSec$bina[!nodata] <- DETall$bina
  
  #
  shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
  shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$bina),]
  return(shpmodelSecNEW)
  
}


# function for both rasters of two direction -------------------------------------------------

uvRasterNorth <- function(MonsType,type){
    
    uidx <- startFNorth(MonsType,type,stacked_USec)
    vidx <- startFNorth(MonsType,type,stacked_VSec)
    
    #combine two index
    fnlIdx <- uidx[[1]] | vidx[[1]]
    
    #Then output tif separately to facilitate subsequent drawing of arrows.
    fnlU <- uidx[[2]]
    fnlU[!fnlIdx] <- NA
    
    fnlV <- vidx[[2]]
    fnlV[!fnlIdx] <- NA
    
    #convert to tif
    tifU <- Mat2Tif(fnlU)
    
    plot(tifU)
    
    tifV <- Mat2Tif(fnlV)
    plot(tifV)
    return(list(tifU,tifV))
}



uvRasterSouth <- function(MonsType,type){
  
  uidx <- startFSouth(MonsType,type,stacked_USec)
  vidx <- startFSouth(MonsType,type,stacked_VSec)

  #combine two index
  fnlIdx <- uidx[[1]] | vidx[[1]]
  
  #Then output tif separately to facilitate subsequent drawing of arrows.
  fnlU <- uidx[[2]]
  fnlU[!fnlIdx] <- NA
  
  fnlV <- vidx[[2]]
  fnlV[!fnlIdx] <- NA
  
  #convert to tif
  tifU <- Mat2Tif(fnlU)
  
  plot(tifU)
  
  tifV <- Mat2Tif(fnlV)
  plot(tifV)
  return(list(tifU,tifV))
}



# begin ----------------------------------------------------------------


shpmodelSecNEW <- getHDetSHP()
plot(shpmodelSecNEW)


time <- seq(as.Date("2002-04-01"),as.Date("2016-08-31"),by="day")


#monsoon boundary
Mons <- st_read("NEW_MonsoonRegions.shp")


# Specialized drawing function construction --- Take NAmerM as an example first -------------------------------------------------------------
#Mapping elements include high-DET followed by monsoon zone boundary followed by q and v direction

#
MonsType <- NAmerM
type <- "NAmerM"


# "NAmerM"
# "WAfriM"
# "SAsiaM"
#"EAsiaM"
# "EqAmerM"
# 
# "SAmerM"
# "SAfriM"
# "AusMCM"

#--fix display boundary--
xmin <- -134
ymin <- 13
xmax <- -77
ymax <- 44
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))

# 
rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 


# 
rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)


#mask precipitation

raster <- brick("NAmerM_WVF_composite_analysis_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)


#get monsoon regions location
locat <- getFldLoc(shpmodelSecNEW,type)

#Then integrate Q, qu and qv into a data frame for subsequent drawing.

Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  #convert it to data.frame
# Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]


#Get two component vectors
vc <- uvRasterNorth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]

#The coefficients here are calculated separately to ensure that all arrows have the same length, which is 1.8 deg.
tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)

# Then merge the data frames
Cuttif_df$xend <- Cuttif_df$x+a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]

#plot
scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}

#global country borders

# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')

# Jbrs <- round(c(-9,-6,-4,0,7,13,19,25,31,37.3),digits = 1)
# Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
# hist(Cuttif_df[[3]])

pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.05, "inches"),ends = "last", type = "open")) +
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




# -WAfriM------------------------------------------------------------------------

MonsType <- WAfriM
type <- "WAfriM"

xmin <- -30
ymin <- -16
xmax <- 65
ymax <- 28
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 



rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)



raster <- brick("WAfriM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)


#
locat <- getFldLoc(shpmodelSecNEW,type)

#

Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 转换为 data frame，并包含坐标信息
# Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]

#
vc <- uvRasterNorth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]

#

tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)

#
Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]

#
scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')

# Jbrs <- round(c(-9,-6,-4,0,7,13,19,25,31,37.3),digits = 1)
# Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
# hist(Cuttif_df[[3]])

pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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
       plot=pp, width = 4, height = 1.8, dpi = 300,units="in")



# SAsiaM ------------------------------------------------------------------

MonsType <- SAsiaM
type <- "SAsiaM"

xmin <- 60
ymin <- 0
xmax <- 135
ymax <- 40
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))

#
rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系


#
rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)


raster <- brick("SAsiaM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)


#
locat <- getFldLoc(shpmodelSecNEW,type)

#

Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 
# Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]

#
vc <- uvRasterNorth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]

#
tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 2/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)

# 
Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]

#
scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')


pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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


# EAsiaM ------------------------------------------------------------------

MonsType <- EAsiaM
type <- "EAsiaM"

xmin <- 90
ymin <- 14
xmax <- 145
ymax <- 50
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 

# 
rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)



raster <- brick("EAsiaM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)



locat <- getFldLoc(shpmodelSecNEW,type)


Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 转换为 data frame，并包含坐标信息


vc <- uvRasterNorth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]


tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)


Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]


scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')

# Jbrs <- round(c(-9,-6,-4,0,7,13,19,25,31,37.3),digits = 1)
# Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
# hist(Cuttif_df[[3]])

pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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



# EqAmerM -----------------------------------------------------------------

MonsType <- EqAmerM
type <- "EqAmerM"

xmin <- -110
ymin <- -10
xmax <- -40
ymax <- 30
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系


rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)



raster <- brick("EqAmerM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)



locat <- getFldLoc(shpmodelSecNEW,type)


Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 转换为 data frame，并包含坐标信息
# Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]

vc <- uvRasterNorth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]


tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)


Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]


scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')

# Jbrs <- round(c(-9,-6,-4,0,7,13,19,25,31,37.3),digits = 1)
# Jbrs <- round(getJenksBreaks(Cuttif_df[[3]], 7, subset = NULL),digits = 1)
# hist(Cuttif_df[[3]])

pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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



# SAmerM-------------------------------------------------------------------


MonsType <- SAmerM
type <- "SAmerM"

xmin <- -100
ymin <- -40
xmax <- -20
ymax <- 9
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 



rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)


raster <- brick("SAmerM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)


locat <- getFldLoc(shpmodelSecNEW,type)


Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 
# Cuttif_df <- Cuttif_df[!is.na(Cuttif_df[[3]]),]

#
vc <- uvRasterSouth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]



tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)


Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]


scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}

# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')

pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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




# SAfriM ------------------------------------------------------------------


MonsType <- SAfriM
type <- "SAfriM"

xmin <- 0
ymin <- -38
xmax <- 55
ymax <- 3
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系



rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)




raster <- brick("SAfriM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)


locat <- getFldLoc(shpmodelSecNEW,type)


Cuttif_df <- as.data.frame(Cuttif, xy = TRUE) 

#
vc <- uvRasterSouth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]


tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)

# 
Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]

#
scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')



pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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


# AusMCM ------------------------------------------------------------------


MonsType <- AusMCM
type <- "AusMCM"

xmin <- 90
ymin <- -30
xmax <- 170
ymax <- 20
rectangle <- st_polygon(list(rbind(c(xmin, ymin), c(xmin, ymax), c(xmax, ymax), c(xmax, ymin), c(xmin, ymin))))


rectangle_sf <- st_sfc(rectangle, crs = st_crs(4326))  # 4326 是一个常用的坐标系，您可以根据需要选择不同的坐标系



rectangle_sp <- as(rectangle_sf, "Spatial")

plot(rectangle_sf)




raster <- brick("AusMCM水汽通量均值相对季风时段差值_0.05.tif")

Cuttif <- crop(raster,rectangle_sp)

plot(Cuttif)



locat <- getFldLoc(shpmodelSecNEW,type)


Cuttif_df <- as.data.frame(Cuttif, xy = TRUE)  # 


vc <- uvRasterSouth(MonsType,type)
tifU <- vc[[1]]
tifV <- vc[[2]]


tifU_cut <-  crop(tifU,rectangle_sp)
tifV_cut <-  crop(tifV,rectangle_sp)
tifU_cutDF <- as.data.frame(tifU_cut, xy = TRUE)
tifV_cutDF <- as.data.frame(tifV_cut, xy = TRUE)

a <- 1.8/sqrt(tifU_cutDF[[3]]^2+tifV_cutDF[[3]]^2)


Cuttif_df$xend <- Cuttif_df$x+ a*tifU_cutDF[[3]]
Cuttif_df$yend <- Cuttif_df$y+ a*tifV_cutDF[[3]]


scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}



# twscols <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))
twscols <- c('#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c')


pp <-
  ggplot(data = Cuttif_df)+
  # geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA)+
  geom_raster(data = Cuttif_df, aes(x = x, y = y, fill = Cuttif_df[[3]])) +
  # scale_fill_distiller_custom(twscols,name = "WVF",breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  scale_fill_gradientn(name = "WVF",colours = twscols,na.value = "transparent") +
  geom_sf(data = countryGlobal,fill=NA,color='black',linewidth = 0.1,alpha=0.9)+
  geom_sf(data = Mons[Mons$Monsoon==type,],fill=NA,color= "#F033FF",alpha=0.9,linewidth = 0.5)+
  geom_sf(data = locat,fill=NA,color= "black",alpha=0.9,size = 0.000001,linewidth = 0.00001)+
  geom_arrowsegment(aes(x = x, xend = xend, y = y, yend = yend),arrows = arrow(length = unit(0.03, "inches"),ends = "last", type = "open")) +
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


























