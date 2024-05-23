
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
library(exploreR)
library(BAMMtools)

library(paletteer)
#modified DET
sf_use_s2(FALSE)

#
DETall <- read.csv(file = "MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv",header = F)
DETall$V1[is.nan(DETall$V1)] <- NA



max(DETall$V1,na.rm = T)
min(DETall$V1,na.rm = T)
hist(DETall$V1)

#add a binary column
DETall$bina <- '0'
DETall$bina[!is.na(DETall$V1)] <- '1'

shpmodelSec$Value <- NA
shpmodelSec$bina <- NA
shpmodelSec$Value[!nodata] <- as.numeric(DETall$V1)
shpmodelSec$bina[!nodata] <- DETall$bina


#remove non-land grids
shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$bina),]


scale_fill_fermenter_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "fermenter", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}

scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


twscols <- c('#ffeda0','#f03b20')
# twscols <- c('#3B9AB2','#F21A00')

pp <-
  ggplot()+
  geom_sf(data = shpmodelSecNEW,aes(fill=bina),color = NA)+
  geom_sf(data = countryGlobal,fill=NA,color='Black',size=0.01,alpha=0.8)+
  # scale_fill_gradientn(name = "Lags",colours = twscols,limits=c(0,12)) +
  # scale_fill_distiller(twscols,name = "Lags",breaks=c(0,1)) +
  scale_fill_manual(values = twscols,na.value = "transparent",labels = c('0'='Low','1'='High'))+
  coord_sf(ylim = c(-60,60),expand = F)+
  # guides(guide = "none")+
  theme(
    text=element_text(family="Helvetica"),
    # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_text(family = "Helvetica",size = 16,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 16,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold"),
    # legend.title = element_text(family = "Helvetica",size = 12,face = "bold")
    legend.title = element_blank()
  )

# pp

ggsave("MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.png",
       plot=pp, width = 8, height = 3, dpi = 300,units="in")



#draw significant DET-------------------------------------------------------------


DETall <- read.csv(file = "MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv",header = F)
DETall$V1[is.nan(DETall$V1)] <- NA


max(DETall$V1,na.rm = T)
min(DETall$V1,na.rm = T)
hist(DETall$V1)


shpmodelSec$Value <- NA
shpmodelSec$Value[!nodata] <- as.numeric(DETall$V1)

#使用讲非陆地边界扣除的数据
shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$Value),]

# twscols <- paletteer_d("dichromat::BluetoOrange_10")[c(1,3,5,7,9,10)]
twscols <- c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026")

#寻找对应的自然间断法的breaks
Jbrs <- round(getJenksBreaks(shpmodelSecNEW$Value, 7, subset = NULL), 2)
# Jbrs <- c(0.76,0.80,0.84,0.88,0.92,0.96,1)

scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}


pp <- 
  ggplot()+
  geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA,size=0.01)+
  geom_sf(data = countryGlobal,fill=NA,color='Black',size=0.01,alpha=0.8)+
  # scale_fill_gradientn(name = "MED_RP_DET",colours = twscols,limits=c(0,max(DETall$V1)),breaks=Jbrs) +
  # scale_fill_viridis_c(option = "A")+
  # scale_fill_gradientn(name = "MED_RP_DET_Ratio_LAM",colours = twscols,limits=c(min(DETall$V1),max(DETall$V1))) +
  scale_fill_distiller_custom(twscols,name = "DET",limits = c(min(Jbrs),max(Jbrs)),breaks=Jbrs) +
  coord_sf(ylim = c(-60,60),expand = F)+
  theme(
    text=element_text(family="Helvetica"),
    # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_text(family = "Helvetica",size = 16,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 16,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold"),
    legend.title = element_text(family = "Helvetica",size = 12,face = "bold")
  )


ggsave("MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst-DET.png",
       plot=pp, width = 8, height = 3, dpi = 300,units="in")





# plot distribution of best Tau -------------------------------------------

DETall <- read.csv(file = "MT1mm_fixP0_005_bestTau.csv",header = F)
DETall$V1[is.nan(DETall$V1)] <- NA

max(DETall$V1,na.rm = T)
min(DETall$V1,na.rm = T)
hist(DETall$V1)

shpmodelSec$Value <- NA
shpmodelSec$Value[!nodata] <- as.numeric(DETall$V1)


#remove non-land grids
shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$Value),]
hist(shpmodelSecNEW$Value)
# twscols <- paletteer_d("tvthemes::Regular")[1:5]
twscols <- rev(paletteer_c("grDevices::YlOrBr", 30))



pp <- 
  ggplot()+
  geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA,size=0.01)+
  geom_sf(data = countryGlobal,fill=NA,color='Black',size=0.01,alpha=0.8)+
  # scale_fill_gradientn(name = "MED_RP_DET",colours = twscols,limits=c(0,max(DETall$V1)),breaks=Jbrs) +
  # scale_fill_viridis_c(option = "A")+
  scale_fill_gradientn(name = "Tau",colours = twscols,limits=c(min(shpmodelSecNEW$Value),max(shpmodelSecNEW$Value))) +
  # scale_fill_distiller_custom(twscols,name = "Tau",breaks= Jbrs) +
  coord_sf(ylim = c(-60,60),expand = F)+
  theme(
    text=element_text(family="Helvetica"),
    # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_text(family = "Helvetica",size = 16,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 16,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold"),
    legend.title = element_text(family = "Helvetica",size = 12,face = "bold")
  )

ggsave("MT1mm_fixP0_005_bestTau.png",
       plot=pp, width = 8, height = 3, dpi = 300,units="in")



# global distribution of average Tau --------------------------------------

DETall <- read.csv(file = "DET_significance_95th/globalfld_Tau_meanTimeInterval.csv",header = F)
DETall$V1[is.nan(DETall$V1)] <- NA

max(DETall$V1,na.rm = T)
min(DETall$V1,na.rm = T)
hist(DETall$V1)

shpmodelSec$Value <- NA
shpmodelSec$Value[!nodata] <- as.numeric(DETall$V1)


#remove non-land grids
shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$Value),]
hist(shpmodelSecNEW$Value)
# twscols <- paletteer_d("tvthemes::Regular")[1:5]
twscols <- rev(paletteer_c("grDevices::YlOrBr", 30))


pp <- 
  ggplot()+
  geom_sf(data = shpmodelSecNEW,aes(fill=Value),color = NA,size=0.01)+
  geom_sf(data = countryGlobal,fill=NA,color='Black',size=0.01,alpha=0.8)+
  # scale_fill_gradientn(name = "MED_RP_DET",colours = twscols,limits=c(0,max(DETall$V1)),breaks=Jbrs) +
  # scale_fill_viridis_c(option = "A")+
  scale_fill_gradientn(name = "Tau",colours = twscols,limits=c(min(shpmodelSecNEW$Value),max(shpmodelSecNEW$Value))) +
  # scale_fill_distiller_custom(twscols,name = "Tau",breaks= Jbrs) +
  coord_sf(ylim = c(-60,60),expand = F)+
  theme(
    text=element_text(family="Helvetica"),
    # plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_text(family = "Helvetica",size = 16,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 16,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold"),
    legend.title = element_text(family = "Helvetica",size = 12,face = "bold")
  )

ggsave("globalfld_Tau_meanTimeInterval.tiff",
       plot=pp, width = 8, height = 3, dpi = 300,units="in")



# compute the histogram of global best Tau, then choose 20th, 50th and 90th quantile values to calculate the distribution--------
#of global DET separately.

bestTau <- read.csv(file = "MT1mm_fixP0_005_bestTau.csv",header = F)
ValidBestTau <- bestTau$V1[!is.nan(bestTau$V1)]

ValidBestTau <- data.frame(Tau=ValidBestTau)
max(ValidBestTau)


pp <- 
  ValidBestTau %>%
  ggplot( aes(x=Tau)) +
  geom_histogram( binwidth=10, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Global optimal Tau distribution") +
  ylab("Count")+
  theme(
    # plot.title = element_text(size=15),
    text=element_text(family="Helvetica"),
    plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 16,face = "bold"),
    axis.title.x = element_text(family = "Helvetica",size = 16,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 16,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 14,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 14,face = "bold")
  )


ggsave("OptimalTauDistribution.png",
       plot=pp, width = 6.21, height = 5.67, dpi = 300,units="in")


#29th quantile 

quan20 <- quantile(ValidBestTau$Tau,0.2)




#50th quantile

quan50 <- quantile(ValidBestTau$Tau,0.5)




#90th quantile

quan90 <- quantile(ValidBestTau$Tau,0.9)































