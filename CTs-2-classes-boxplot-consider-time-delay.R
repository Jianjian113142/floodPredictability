
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
# library(XploreR)
library(BAMMtools)
library(readxl)
library(paletteer)
sf_use_s2(FALSE)



##Note that:
#below is a examople for NINO34 climate index (CTs). other CTs will perform the same process

getNINO34 <- function(delay){
  NINO34 <- read_xlsx("Nino34/nino34.xlsx")
  NINO34.sec <- NINO34[NINO34$YEAR>=2000&NINO34$YEAR<=2016,-1]
  NINO34.sec.vec <- as.numeric(t(NINO34.sec))
  
  DATE <- seq(as.Date("2000-01-01"),as.Date("2016-12-31"),by="month")
  
  indx1 <- which(as.Date("2002-04-01")==DATE)
    
  indx2 <- which(as.Date("2016-08-01")==DATE)
  
  return(NINO34.sec.vec[(indx1-delay):(indx2-delay)])
}

#remove seasonal influence
remSea <- function(Flds_dat_Mon){
  
  TT <- seq(as.Date("2002-04-01"),as.Date("2016-08-31"),by="month")
  
  result1 <- rep(NA,length(TT))
  
  for (i in 1:12){
    
    meanTmp <- mean(Flds_dat_Mon[month(TT)==i])
    
    result1[month(TT)==i] <- Flds_dat_Mon[month(TT)==i] - meanTmp
    
  }
  
  # result2 <- detrend(result1)
  
  return(result1)
}




#Find the corresponding final time index according to the two categories.
find_matching_dates <- function(A, B) {
  # Extract year and month from A
  A_year_month <- format(A, "%Y-%m")
  
  # Extract year and month from B
  B_year_month <- format(B, "%Y-%m")
  
  # Find the date in A that has the same year and month as B
  matching_dates <- A[A_year_month %in% B_year_month]
  
  return(matching_dates)
}


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
  
  #
  shpmodelSecNEW <- shpmodelSec[match(outputModel$id,shpmodelSec$id),]
  shpmodelSecNEW <- shpmodelSecNEW[!is.na(shpmodelSecNEW$bina),]
  return(shpmodelSecNEW)
  
}

##get global HFP distribution
getFldLoc <- function(Model,type){
  
  MonsSec <- Mons[Mons$Monsoon==type,]
  plot(MonsSec)
  
  inter <- st_intersection(Model,MonsSec)
  
  return(inter)
  
}


# only consider high-DET ----------------------------------------------

North_highDET <- function(type,TComp,Osilation){
  # type <- "SAsiaM"
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  #date index
  
  NorthMonsIdx <- match(TComp,dateTws) #
  
  #spatial index
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)

  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)#
  
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep(Osilation,length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  return(high_comb)
  
}




South_highDET <- function(type,TComp,Osilation){ #Osilation is the character replacement that indicates whether it is positive or negative.
  # type <- "SAsiaM"
  MonsDET <- getFldLoc(shpmodelSecNEW,type)
  
  plot(MonsDET)
  
  NorthMonsIdx <- match(TComp,dateTws) #
  
  #需要的空间索引
  spatial_idx1 <- match(MonsDET$id[MonsDET$bina=='1'],id13882)

  # dim(floodevents_morethan1mm_95th_2002_2016)
  high_rt <- apply(floodevents_morethan1mm_95th_2002_2016[NorthMonsIdx,spatial_idx1], 2, sum)#
  
  high_mean <- mean(high_rt)
  
  high_chr <- rep(Osilation,length(high_rt))
  
  high_comb <- cbind(high_rt,high_chr)
  dim(high_comb)
  
  return(high_comb)
  
}



plotFunction <- function(EL,LA,type){
  
  final <- rbind(EL,LA)
  
  dtt <- data.frame(Fraction=as.numeric(final[,1]),Type = final[,2])
  
  dtt$Type <- factor(dtt$Type, levels = c("Positive","Negative"))
  
  p <- dtt %>%
    ggplot() +
    geom_boxplot(aes(x=Type, y=Fraction , fill=Type),color = "#2F2F2F",width = 0.2,alpha=0.6,size = 0.2) +
    scale_fill_manual(values=c("#fc8d59","#91bfdb"))+
    # scale_color_manual(values=c("#A53D25",""))+

    labs(fill="")+
    ylab("Flood events")+
    xlab("Category")+
    # ylim(0,1)+
    ggtitle(type)+
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
  
  
  output <- paste(type,"_Boxplot-",".png",sep="")
  
  ggsave(output,
         plot=p, width = 3.75, height = 3, dpi = 300,units="in")
  
  
}



Fun_delay <- function(type,DHSFEMT){
  #Create a new matrix to calculate the results of each monsoon zone-----------------------------
  results <- data.frame(NAME = rep(NA,13),MEAN = rep(NA,13),p = rep(NA,13),delay = rep(NA,13)) #Consider the p-value corresponding to a delay of 0-12 months
  for (i in 0:12){
    
    #Calculate the p-value under different delay sequences
    NINO34 <- remSea(getNINO34(i))
    #months
    TT <-  seq(as.Date("2002-04-01"),as.Date("2016-8-31"),by="month")
    
    #find El Niño months
    TT_El <- TT[NINO34>=0]
    #find La Niña months
    TT_La <- TT[NINO34< 0]
    
    El_Comp <- find_matching_dates(DHSFEMT,TT_El)
    
    La_Comp <- find_matching_dates(DHSFEMT,TT_La)
    
    EL <- North_highDET(type,El_Comp,"Positive")
    
    LA <- North_highDET(type,La_Comp,"Negative")
    
    TTST <- t.test(as.numeric(EL[,1]),as.numeric(LA[,1]))
    
    results$NAME[i+1] <- type
    results$MEAN[i+1] <- mean(as.numeric(EL[,1])) - mean(as.numeric(LA[,1]))
    results$p[i+1] <- TTST$p.value
    results$delay[i+1] <- i
    # plotFunction(EL,LA,type)
  }
  
  #Store results for a certain monsoon zone
  write.csv(x = results,file = paste(type,".csv",sep = ""),row.names = F)
  
}
#for south hemisphere
Fun_delay_south <- function(type,DHSFEMT){
  #
  results <- data.frame(NAME = rep(NA,13),MEAN = rep(NA,13),p = rep(NA,13),delay = rep(NA,13)) #
  for (i in 0:12){
    
    #
    NINO34 <- remSea(getNINO34(i))
    #
    TT <-  seq(as.Date("2002-04-01"),as.Date("2016-8-31"),by="month")
    
    #
    TT_El <- TT[NINO34>=0]
    #
    TT_La <- TT[NINO34< 0]
    
    El_Comp <- find_matching_dates(DHSFEMT,TT_El)
    
    La_Comp <- find_matching_dates(DHSFEMT,TT_La)
    
    EL <- South_highDET(type,El_Comp,"Positive")
    
    LA <- South_highDET(type,La_Comp,"Negative")
    
    TTST <- t.test(as.numeric(EL[,1]),as.numeric(LA[,1]))
    
    results$NAME[i+1] <- type
    results$MEAN[i+1] <- mean(as.numeric(EL[,1])) - mean(as.numeric(LA[,1]))
    results$p[i+1] <- TTST$p.value
    results$delay[i+1] <- i
    # plotFunction(EL,LA,type)
  }
  
  write.csv(x = results,file = paste(type,".csv",sep = ""),row.names = F)
  
}



##Plot based on the found optimal p-value and corresponding intermediate data

FinalPlt <- function(i,type,DHSFEMT){
  #Calculate the p-value under different delay sequences
  NINO34 <- remSea(getNINO34(i))
  #months
  TT <-  seq(as.Date("2002-04-01"),as.Date("2016-8-31"),by="month")
  
  # find El Niño months
  TT_El <- TT[NINO34>=0]
  #find La Niña months
  TT_La <- TT[NINO34< 0]
  
  
  El_Comp <- find_matching_dates(DHSFEMT,TT_El)
  
  La_Comp <- find_matching_dates(DHSFEMT,TT_La)
  
  EL <- North_highDET(type,El_Comp,"Positive")
  
  LA <- North_highDET(type,La_Comp,"Negative")
  
  TTST <- t.test(as.numeric(EL[,1]),as.numeric(LA[,1]))
  
  plotFunction(EL,LA,type)
  
}




# begin -------------------------------------------------------------------

#
Mons <- st_read("NEW_MonsoonRegions.shp")

shpmodelSecNEW <- getHDetSHP()


#calculate WAfriM
Fun_delay("WAfriM",WAfriM)


#calculate SAsiaM
Fun_delay("SAsiaM",SAsiaM)

#calculate EAsiaM
Fun_delay("EAsiaM",EAsiaM)

#calculate NAmerM
Fun_delay("NAmerM",NAmerM)

#calculate EqAmerM
Fun_delay("EqAmerM",EqAmerM)


#calculate SAmerM
Fun_delay_south("SAmerM",SAmerM)

#calculate SAfriM
Fun_delay_south("SAfriM",SAfriM)

#calculate AusMCM
Fun_delay_south("AusMCM",AusMCM)




# save results ---------------------------------------------------------------
# results$type <- rep("Nino3.4",8)
# write.csv(x = results,file = "Nino34_mean_p.csv",row.names = F)
# 
# 

# draw codes ---------boxplot----------------------------------------------------------

FinalPlt(i = 0,type = "AusMCM",AusMCM)

FinalPlt(i = 7,type = "EAsiaM",EAsiaM)

FinalPlt(i = 3,type = "EqAmerM",EqAmerM)

FinalPlt(i = 2,type = "NAmerM",NAmerM)

FinalPlt(i = 5,type = "SAfriM",SAfriM)

FinalPlt(i = 10,type = "SAmerM",SAmerM)

FinalPlt(i = 3,type = "SAsiaM",SAsiaM)

FinalPlt(i = 4,type = "WAfriM",WAfriM)










