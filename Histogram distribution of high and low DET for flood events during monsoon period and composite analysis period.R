
library(raster)
library(tools)
library(zoo)
library(seasonal)
library(stlplus)
library(lubridate)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(readxl)
library(sf)
library(extrafont)


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
library(tidyr)
library('ggedit')
library(ggplot2)
library(paletteer)
library(exploreR)
library(BAMMtools)
library(readxl)
library(paletteer)
sf_use_s2(FALSE)



# get data model ------------------------------------------------------------------

getHDetSHP <- function(){
  input <- "~/Desktop/洪水日数95th-morethan1mm/DET_significance_95th/MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv"
  
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

# get high DET grids in monsoon regions -----------------------------------

getFldLoc <- function(Model,type){
  
  MonsSec <- Mons[Mons$Monsoon==type,]
  plot(MonsSec)
  
  inter <- st_intersection(Model,MonsSec)
  
  return(inter)
  
}





#Obtain the global spatial distribution of high and low predictability

shpmodelSecNEW <- getHDetSHP()


#define time
TT <- seq(as.Date("2002-04-01"),as.Date("2016-08-31"),by="day")

#get boundary of monsoon regions
Mons <- st_read("NEW_MonsoonRegions.shp")



# The number of floods in the composite period VS. in the  monsoon period. --------

North_CompositePeriods <- function(type,MonsT){
  # type <- "SAsiaM"
  
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #date index
  NorthMonsIdx <- match(MonsT,TT) #
  #spatial location index
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)
  spatial_idx0 <- match(MonsDET$id[MonsDET$bina=='0'],id13882)
  # length(spatial_idx1)
  
  lgh1 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% 6:9,spatial_idx1], 2, sum)+exp(-20)
  lgh0 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% 6:9,spatial_idx0], 2, sum)+exp(-20)
  
  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)/lgh1
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep("High",length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  low_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx0], 2, sum)/lgh0
  low_mean <- mean(low_rt)
  
  low_chr <-  rep("Low",length(low_rt))
  low_comb <- cbind(low_rt,low_chr)
  
  final <- rbind(high_comb,low_comb)
  
  
  
  dtt <- data.frame(Frequency=as.numeric(final[,1]),Type = final[,2])
  
  p <- dtt %>%
    ggplot( aes(x=Frequency, fill=Type)) +
    geom_density( aes(y = ..scaled..),color="#e9ecef", alpha=0.6) +
    scale_fill_manual(values=c("#404080","#69b3a2"))+
    geom_vline(xintercept = high_mean, linetype = "dashed", color = "#4A3D8D", size = 1) +  # Add mean line
    geom_vline(xintercept = low_mean, linetype = "dashed", color = "#22763F", size = 1) +  # Add mean line
    
    labs(fill="")+
    ylab("Density")+
    xlab("Fraction")+
    ggtitle(type)+
    theme(
      text=element_text(family="Helvetica"),
      axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
      axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
      axis.text.x = element_text(family = "Helvetica",size = 10,face = "bold"),
      axis.text.y = element_text(family = "Helvetica",size = 10,face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
    )
  p
  
  output <- paste(type,".png",sep="")
  
  ggsave(output,
         plot=p, width = 4.79, height = 3.67, dpi = 300,units="in")
  
  
}



South_CompositePeriods <- function(type,MonsT){
  # type <- "SAsiaM"
  
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #get date index
  
  NorthMonsIdx <- match(MonsT,TT) #

  #spatial index
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)
  spatial_idx0 <- match(MonsDET$id[MonsDET$bina=='0'],id13882)
  # length(spatial_idx1)
  
  lgh1 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% c(12,1,2,3),spatial_idx1], 2, sum)+exp(-20)
  lgh0 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% c(12,1,2,3),spatial_idx0], 2, sum)+exp(-20)
  
  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)/lgh1
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep("High",length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  low_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx0], 2, sum)/lgh0
  
  low_mean <- mean(low_rt)
  low_chr <-  rep("Low",length(low_rt))
  low_comb <- cbind(low_rt,low_chr)
  
  final <- rbind(high_comb,low_comb)
  
  dtt <- data.frame(Frequency=as.numeric(final[,1]),Type = final[,2])
  
  
  
  p <- dtt %>%
    ggplot( aes(x=Frequency, fill=Type)) +
    geom_density( aes(y = ..scaled..), color="#e9ecef", alpha=0.6) +
    scale_fill_manual(values=c("#404080","#69b3a2"))+
    geom_vline(xintercept = high_mean, linetype = "dashed", color = "#4A3D8D", size = 1) +  # Add mean line
    geom_vline(xintercept = low_mean, linetype = "dashed", color = "#22763F", size = 1) +  # Add mean line
    labs(fill="")+
    ylab("Density")+
    xlab("Fraction")+
    ggtitle(type)+
    theme(
      text=element_text(family="Helvetica"),
      axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
      axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
      axis.text.x = element_text(family = "Helvetica",size = 10,face = "bold"),
      axis.text.y = element_text(family = "Helvetica",size = 10,face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
    )
  # p
  
  output <- paste(type,".png",sep="")
  
  ggsave(output,
         plot=p, width = 4.79, height = 3.67, dpi = 300,units="in")
  
  
}


# font_import(pattern = "Helvetica")

North_CompositePeriods("NAmerM",NAmerM)
North_CompositePeriods("EqAmerM",EqAmerM)
North_CompositePeriods("WAfriM",WAfriM)
North_CompositePeriods("SAsiaM",SAsiaM)
North_CompositePeriods("EAsiaM",EAsiaM)

South_CompositePeriods("SAmerM",SAmerM)
South_CompositePeriods("SAfriM",SAfriM)
South_CompositePeriods("AusMCM",AusMCM)



length(NAmerM)
length(EqAmerM)
length(WAfriM)
length(SAsiaM)
length(EAsiaM)

length(SAmerM)
length(SAfriM)
length(AusMCM)




# draw boxplot --------------------------------------------------------------

North_boxplot <- function(type,MonsT){
  # type <- "SAsiaM"
  
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #get date index
  
  NorthMonsIdx <- match(MonsT,TT) 
  
  #spatial index
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)
  spatial_idx0 <- match(MonsDET$id[MonsDET$bina=='0'],id13882)
  # length(spatial_idx1)
  
  lgh1 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% 6:9,spatial_idx1], 2, sum)+exp(-20)
  lgh0 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% 6:9,spatial_idx0], 2, sum)+exp(-20)
  
  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)/lgh1
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep("High",length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  low_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx0], 2, sum)/lgh0
  low_mean <- mean(low_rt)
  
  low_chr <-  rep("Low",length(low_rt))
  low_comb <- cbind(low_rt,low_chr)
  
  final <- rbind(high_comb,low_comb)
  
  dtt <- data.frame(Frequency=as.numeric(final[,1]),Type = final[,2])
  
  dtt$Type <- factor(dtt$Type, levels = c("Low", "High"))
  
  TTST <- t.test(dtt$Frequency[dtt$Type=="Low"],dtt$Frequency[dtt$Type=="High"])
  print(TTST)
  
  p <- dtt %>%
    ggplot( aes(x=Frequency, fill=Type)) +
    geom_boxplot(aes(x=Type, y=Frequency , fill=Type),color = "#2F2F2F",width = 0.2,alpha=0.6,size = 0.2) +
    scale_fill_manual(values=c("#69b3a2","#404080"))+
    # scale_color_manual(values=c("#A53D25",""))+
    # geom_vline(xintercept = high_mean, linetype = "dashed", color = "#4A3D8D", size = 1) +  # Add mean line
    # geom_vline(xintercept = low_mean, linetype = "dashed", color = "#22763F", size = 1) +  # Add mean line
    # 
    labs(fill="")+
    ylab("Frequency")+
    xlab("Predictability")+
    # ggtitle(type)+
    theme(
      text=element_text(family="Helvetica"),
      axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
      axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
      axis.text.x = element_text(family = "Helvetica",size = 10,face = "bold"),
      axis.text.y = element_text(family = "Helvetica",size = 10,face = "bold"),
      # legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
    )

  
  output <- paste(type,".png",sep="")
  
  ggsave(output,
         plot=p, width = 2.5, height = 2, dpi = 300,units="in")
  
  
}



South_boxplot <- function(type,MonsT){
  # type <- "SAsiaM"
  
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #get date index
  
  NorthMonsIdx <- match(MonsT,TT) 
  #spatial index
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)
  spatial_idx0 <- match(MonsDET$id[MonsDET$bina=='0'],id13882)
  # length(spatial_idx1)
  
  lgh1 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% c(12,1,2,3),spatial_idx1], 2, sum)+exp(-20)
  lgh0 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% c(12,1,2,3),spatial_idx0], 2, sum)+exp(-20)
  
  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)/lgh1
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep("High",length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  low_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx0], 2, sum)/lgh0
  
  low_mean <- mean(low_rt)
  low_chr <-  rep("Low",length(low_rt))
  low_comb <- cbind(low_rt,low_chr)
  
  final <- rbind(high_comb,low_comb)
  
  dtt <- data.frame(Frequency=as.numeric(final[,1]),Type = final[,2])
  dtt$Type <- factor(dtt$Type, levels = c("Low", "High"))
  
  TTST <- t.test(dtt$Frequency[dtt$Type=="Low"],dtt$Frequency[dtt$Type=="High"])
  print(TTST)
  
  
  p <- dtt %>%
    ggplot( aes(x=Frequency, fill=Type)) +
    # geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +
    geom_boxplot(aes(x=Type, y=Frequency , fill=Type),color = "#2F2F2F",width = 0.2,alpha=0.6,size = 0.2) +
    scale_fill_manual(values=c("#69b3a2","#404080"))+
    # scale_fill_manual(values=c("#404080","#69b3a2"))+
    # geom_vline(xintercept = high_mean, linetype = "dashed", color = "#4A3D8D", size = 1) +  # Add mean line
    # geom_vline(xintercept = low_mean, linetype = "dashed", color = "#22763F", size = 1) +  # Add mean line
    labs(fill="")+
    ylab("Frequency")+
    xlab("Predictability")+
    # ggtitle()+
    theme(
      text=element_text(family="Helvetica"),
      axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
      axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
      axis.text.x = element_text(family = "Helvetica",size = 10,face = "bold"),
      axis.text.y = element_text(family = "Helvetica",size = 10,face = "bold"),
      # legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
    )
  # p
  
  output <- paste(type,".png",sep="")
  
  ggsave(output,
         plot=p, width = 2.5, height = 2, dpi = 300,units="in")
  
  
}


North_boxplot("NAmerM",NAmerM)
North_boxplot("EqAmerM",EqAmerM)
North_boxplot("WAfriM",WAfriM)
North_boxplot("SAsiaM",SAsiaM)
North_boxplot("EAsiaM",EAsiaM)

South_boxplot("SAmerM",SAmerM)
South_boxplot("SAfriM",SAfriM)
South_boxplot("AusMCM",AusMCM)






# The proportion of flood events in all monsoon areas in the composite period relative to the monsoon period is plotted together as a PDF --------------------------------------

#plot function

getNorth_CompositePeriods <- function(type,MonsT){
  
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #
  
  NorthMonsIdx <- match(MonsT,TT) 
  

  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)
  spatial_idx0 <- match(MonsDET$id[MonsDET$bina=='0'],id13882)
  # length(spatial_idx1)
  
  lgh1 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% 6:9,spatial_idx1], 2, sum)+exp(-20)
  lgh0 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% 6:9,spatial_idx0], 2, sum)+exp(-20)
  
  
  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)/lgh1
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep("High",length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  low_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx0], 2, sum)/lgh0
  low_mean <- mean(low_rt)
  
  low_chr <-  rep("Low",length(low_rt))
  low_comb <- cbind(low_rt,low_chr)
  
  final <- rbind(high_comb,low_comb)

  dtt <- data.frame(Frequency=as.numeric(final[,1]),Type = final[,2])
  
  return(dtt)
}


getSouth_CompositePeriods <- function(type,MonsT){
  # type <- "SAsiaM"
  
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #
  
  NorthMonsIdx <- match(MonsT,TT) 
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)
  spatial_idx0 <- match(MonsDET$id[MonsDET$bina=='0'],id13882)
  # length(spatial_idx1)
  
  lgh1 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% c(12,1,2,3),spatial_idx1], 2, sum)+exp(-20)
  lgh0 <- apply(floodevents_morethan1mm_95th_2002_2016[month(TT) %in% c(12,1,2,3),spatial_idx0], 2, sum)+exp(-20)
  
  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)/lgh1
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep("High",length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  low_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx0], 2, sum)/lgh0
  
  low_mean <- mean(low_rt)
  low_chr <-  rep("Low",length(low_rt))
  low_comb <- cbind(low_rt,low_chr)
  
  final <- rbind(high_comb,low_comb)
  
  dtt <- data.frame(Frequency=as.numeric(final[,1]),Type = final[,2])
  return(dtt)
}


Na <- getNorth_CompositePeriods("NAmerM",NAmerM)
Eq <- getNorth_CompositePeriods("EqAmerM",EqAmerM)
Wa <- getNorth_CompositePeriods("WAfriM",WAfriM)
Sa <- getNorth_CompositePeriods("SAsiaM",SAsiaM)
Ea <- getNorth_CompositePeriods("EAsiaM",EAsiaM)

Sam <- getSouth_CompositePeriods("SAmerM",SAmerM)
Saf <- getSouth_CompositePeriods("SAfriM",SAfriM)
Au <- getSouth_CompositePeriods("AusMCM",AusMCM)

dt_all <- rbind(Na,Eq,Wa,Sa,Ea,Sam,Saf,Au)


low_mean <- mean(dt_all$Frequency[dt_all$Type=="Low"])
high_mean <- mean(dt_all$Frequency[dt_all$Type=="High"])



p <- dt_all %>%
  ggplot( aes(x=Frequency, fill=Type)) +
  geom_density( aes(y = ..scaled..),color="#e9ecef", alpha=0.6) +
  scale_fill_manual(values=c("#404080","#69b3a2"))+
  geom_vline(xintercept = high_mean, linetype = "dashed", color = "#4A3D8D", size = 1) +  # Add mean line
  geom_vline(xintercept = low_mean, linetype = "dashed", color = "#22763F", size = 1) +  # Add mean line
  
  labs(fill="")+
  ylab("Density")+
  xlab("Fraction")+
  ggtitle("All Monsoon Regions")+
  theme(
    text=element_text(family="Helvetica"),
    axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 10,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 10,face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
  )
p

output <- ("PDF_all_monsoon_regions.png")

ggsave(output,
       plot=p, width = 4.79, height = 3.67, dpi = 300,units="in")


#接着计算boxplot

dt_all$Type <- factor(dt_all$Type, levels = c("Low", "High"))

TTST <- t.test(dt_all$Frequency[dt_all$Type=="Low"],dt_all$Frequency[dt_all$Type=="High"])
print(TTST)

p <- dt_all %>%
  ggplot( aes(x=Frequency, fill=Type)) +
  # geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +
  geom_boxplot(aes(x=Type, y=Frequency , fill=Type),color = "#2F2F2F",width = 0.2,alpha=0.6,size = 0.2) +
  scale_fill_manual(values=c("#69b3a2","#404080"))+
  labs(fill="")+
  ylab("Frequency")+
  xlab("Predictability")+
  # ggtitle()+
  theme(
    text=element_text(family="Helvetica"),
    axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 10,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 10,face = "bold"),
    # legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
  )
p

output <- ("boxplot_all_monsoon_regions.png")

ggsave(output,
       plot=p, width = 2.5, height = 2, dpi = 300,units="in")



