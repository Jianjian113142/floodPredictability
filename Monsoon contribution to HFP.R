
library(raster)
library(tools)
library(zoo)
library(seasonal)
library(stlplus)
library(lubridate)
library(tidyverse)
library(ncdf4)
library(sf)
library(parallel)
theme_set(theme_bw())
library(zoo)
library(tools)
# library(tidyr)
library('ggedit')
library(ggplot2)
library(paletteer)
library(exploreR)
library(BAMMtools)
library(readxl)
library(paletteer)
sf_use_s2(FALSE)
library(lubridate)

library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)


#read monsoon boundary
Mons <- st_read("NEW_MonsoonRegions.shp")


input <- "MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv"

DETall <- read.csv(file = input,header = F)
DETall$V1[is.nan(DETall$V1)] <- NA

#add a binary var.
DETall$bina <- '0'
DETall$bina[!is.na(DETall$V1)] <- '1'

shpmodelSec$bina <- NA
shpmodelSec$bina[!nodata] <- DETall$bina

#remove non-land grid points
shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]



#consider spring and summer separately for North and South hemisphere

fractionForWinter <- function(Model,MonsSec){
  inter <- st_intersection(Model,MonsSec)
  
  inter1 <- inter[inter$bina=="1",]
  # inter0 <- inter[inter$bina=="0",]
  
  #get index of grid points
  idx1 <- match(inter1$id,as.numeric(id13882))
  fldSEC1 <- floodevents_morethan1mm_95th_2002_2016[,idx1]
  
  OneVal <- c()
  
  for (i in 1:dim(fldSEC1)[2]){
    
    time <- seq(as.Date("2002-04-01"),as.Date("2016-08-31"),by="day")
    
    TimeSec <- time[fldSEC1[,i]>=1]
    
    OneVal <- c(OneVal,sum((month(TimeSec) %in% c(12,1,2,3)))/length(TimeSec))
    
  }
  return(OneVal)
} 



fractionForEachSummer <- function(Model,MonsSec){
  inter <- st_intersection(Model,MonsSec)
  
  inter1 <- inter[inter$bina=="1",]
  # inter0 <- inter[inter$bina=="0",]
  
  #get drid point index
  idx1 <- match(inter1$id,as.numeric(id13882))
  fldSEC1 <- floodevents_morethan1mm_95th_2002_2016[,idx1]
  
  OneVal <- c()
  
  for (i in 1:dim(fldSEC1)[2]){
    
    time <- seq(as.Date("2002-04-01"),as.Date("2016-08-31"),by="day")
    
    TimeSec <- time[fldSEC1[,i]>=1]
    
    OneVal <- c(OneVal,sum((month(TimeSec) %in% c(6,7,8,9)))/length(TimeSec))
    
  }
  
  return(OneVal)
  
} 




DRAW <- function(type){
  
  # type <- "NAmerM"
  MonsSec <- Mons[Mons$Monsoon==type,]
  # plot(MonsSec)
  
  if (type == "SAmerM" || type == "SAfriM" || type == "AusMCM"){
  OneVal <- fractionForWinter(shpmodelSecNEW,MonsSec)}else{
    
  OneVal <- fractionForEachSummer(shpmodelSecNEW,MonsSec)
  }
  
  # # ZeroVal <- fractionForEach0(shpmodelSecNEW)
  # hist(OneVal)
  # sum(OneVal>=0.5)/length(OneVal)
  
  OneValDF <- data.frame(Fraction =OneVal)
  
  OneValDF$type <- type
  
  return(OneValDF)
  
}


DF_AusMCM <- DRAW("AusMCM")
DF_AusMCM$contribution <- sum(DF_AusMCM$Fraction>=0.5)/dim(DF_AusMCM)[1] *100


DF_SAmerM <- DRAW("SAmerM")
DF_SAmerM$contribution <- sum(DF_SAmerM$Fraction>=0.5)/dim(DF_SAmerM)[1]*100


DF_SAfriM <- DRAW("SAfriM")
DF_SAfriM$contribution <- sum(DF_SAfriM$Fraction>=0.5)/dim(DF_SAfriM)[1]*100

DF_EAsiaM <- DRAW("EAsiaM")
DF_EAsiaM$contribution <- sum(DF_EAsiaM$Fraction>=0.5)/dim(DF_EAsiaM)[1]*100

DF_SAsiaM <- DRAW("SAsiaM")
DF_SAsiaM$contribution <- sum(DF_SAsiaM$Fraction>=0.5)/dim(DF_SAsiaM)[1]*100

DF_NAmerM <- DRAW("NAmerM")
DF_NAmerM$contribution <- sum(DF_NAmerM$Fraction>=0.5)/dim(DF_NAmerM)[1]*100

DF_EqAmerM <- DRAW("EqAmerM")
DF_EqAmerM$contribution <- sum(DF_EqAmerM$Fraction>=0.5)/dim(DF_EqAmerM)[1]*100

DF_WAfriM <- DRAW("WAfriM")
DF_WAfriM$contribution <- sum(DF_WAfriM$Fraction>=0.5)/dim(DF_WAfriM)[1]*100


frac_comb <- rbind(DF_NAmerM,DF_EAsiaM,DF_SAfriM,DF_SAsiaM,DF_WAfriM,DF_SAmerM,DF_EqAmerM,DF_AusMCM)

frac_comb$type <- factor(frac_comb$type, levels = rev(c("NAmerM","EAsiaM","SAfriM","SAsiaM","WAfriM","SAmerM","EqAmerM","AusMCM")))

ggplot(frac_comb, aes(x = Fraction, y = type, fill = contribution)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Contribution(%)", option = "C") +
  # labs(title = 'Temperatures in Lincoln NE in 2016') +
  # theme_ipsum() +
  xlim(0,1)+
  ylab("Monsoon regions")+
  theme(
    # legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    axis.title.x = element_text(family = "Helvetica",size = 16,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 16,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold"),
    legend.title = element_text(family = "Helvetica",size = 12,face = "bold")
  )


















