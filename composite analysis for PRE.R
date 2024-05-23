
# Here we focus on analyzing the significance level of the precipitation mean corresponding to the extracted time period relative to the mean value of the entire monsoon period.------------------------------------------------------------------


dateTws <- seq(as.Date('2002-04-01'), as.Date('2016-8-31'), by='day')

filesNew <- list_files_with_exts("H:/GPM_day/geotif/GPM2002_2017_day","tif")



#Find the monsoon period, which is mainly from June to September in the Northern Hemisphere, and from December to March in the Southern Hemisphere.
loc2 <- month(dateTws) %in% c(6,7,8,9)
loc3 <-  month(dateTws) %in% c(12,1,2,3)



##mean and std in Northern hemisphere monsoon period
MSD2 <- Mean_STD(filesNew,loc2)
Mean_pre2 <- MSD2[[1]]
std_pre2 <- MSD2[[2]]

#mean and std in Southern Hemisphere Monsoon period
MSD3 <- Mean_STD(filesNew,loc3)
Mean_pre3 <- MSD3[[1]]
std_pre3 <- MSD3[[2]]



# begin ------------------------------------------------------------
NorthComp(NAmerM,Mean_pre2,std_pre2,"NAmerM")
NorthComp(WAfriM,Mean_pre2,std_pre2,"WAfriM")
NorthComp(SAsiaM,Mean_pre2,std_pre2,"SAsiaM")
NorthComp(EAsiaM,Mean_pre2,std_pre2,"EAsiaM")
NorthComp(EqAmerM,Mean_pre2,std_pre2,"EqAmerM")

SouthComp(SAmerM,Mean_pre3,std_pre3,"SAmerM")
SouthComp(SAfriM,Mean_pre3,std_pre3,"SAfriM")
SouthComp(AusMCM,Mean_pre3,std_pre3,"AusMCM")


#Results calculated for NAmerM-------------------------------------------


#Function calculated for the Northern Hemisphere

NorthComp <- function(MonsT,Mean_pre2,std_pre2,type){

loc <- match(MonsT,dateTws)

MSD <- Mean_STD(filesNew,loc)
Mean_pre1 <- MSD[[1]]
std_pre1 <- MSD[[2]]

n1 <- length(loc)
n2 <- sum(loc2)
CriticalV <- critical_value(n1,n2)

Mons_Test <- TTest(Mean_pre1,std_pre1,Mean_pre2,std_pre2,n1,n2)

deltaPre <- Mean_pre1 - Mean_pre2

ModelTIf <- raster(deltaPre)

ModelTIf[Mons_Test>CriticalV] <- deltaPre[Mons_Test>CriticalV]

plot(ModelTIf)

writeRaster(x = ModelTIf,filename = paste(type,"降水均值相对季风时段差值.tif",sep=""))
}



#Functions calculated for the Southern Hemisphere

SouthComp <- function(MonsT,Mean_pre3,std_pre3,type){
  loc <- match(MonsT,dateTws)
  
  MSD <- Mean_STD(filesNew,loc)
  Mean_pre1 <- MSD[[1]]
  std_pre1 <- MSD[[2]]
  
  n1 <- length(loc)
  n2 <- sum(loc3)
  CriticalV <- critical_value(n1,n2)
  
  Mons_Test <- TTest(Mean_pre1,std_pre1,Mean_pre3,std_pre3,n1,n2)
  
  deltaPre <- Mean_pre1 - Mean_pre3 
  
  ModelTIf <- raster(deltaPre)
  
  ModelTIf[Mons_Test>CriticalV] <- deltaPre[Mons_Test>CriticalV]
  
  plot(ModelTIf)
  
  writeRaster(x = ModelTIf,filename = paste(type,"降水均值相对季风时段差值.tif",sep=""))
}




# Function: Calculate whether each grid point is significant based on degrees of freedom and threshold ------------------------------------------------------


TTest <- function(Mean1,Std1,Mean2,Std2,n1,n2){
  
  Ttest <- (Mean1-Mean2)/(sqrt((n1*Std1^2+n2*Std2^2)/(n1+n2-2))*sqrt(1/n1+1/n2))
  
  return(Ttest)
}







# Function: calculation of critical value---------------------------------------------------------------
critical_value <- function(n1,n2){

  df <- n1 + n2 - 2  
  alpha <- 0.05  
  
  critical_value <- qt(1 - alpha, df)
  return(critical_value)
}



# mean and std function -----------------------------------------------------------------

Mean_STD <- function(filesNew,loc){
  
  filesNewSec <- filesNew[loc]
  GPM_MOdel <- brick(filesNewSec[1])
  values(GPM_MOdel) <- 0
  
  sum_of_squares <- GPM_MOdel 
  sum_values <- GPM_MOdel
  
  
  for (i in 1:length(filesNewSec)){
    
    
    GPM <- brick(filesNewSec[i])

    sum_of_squares <- sum_of_squares + GPM^2
    sum_values <- sum_values + GPM
  
  }

  mean_values <- sum_values / length(filesNewSec)
  std_dev <- sqrt((sum_of_squares / length(filesNewSec)) - (mean_values^2))
  
  return(list(mean_values,std_dev))
  
}













