
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



osciFIles <- list_files_with_exts("heatmap_data","csv")

finalRt <- data.frame()

for (i in 1:length(osciFIles)){
  
  tmp <- read.csv(osciFIles[i])
  
  finalRt <- rbind(finalRt,tmp)
  
  
}


# Determine whether it is significant based on PVALUE and distinguish between different significance levels
finalRt$Significance <- ifelse(finalRt$p < 0.001, "***", 
                            ifelse(finalRt$p < 0.01, "**", 
                                   ifelse(finalRt$p < 0.05, "*", "")))

Color <- rev(c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac'))
Jbrs <- round(getJenksBreaks(finalRt$MEAN, 7, subset = NULL),digits = 2)



scale_fill_distiller_custom <- function(pal, na.value = "transparent", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "distiller", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)
}
# plot heatmap
pp <- ggplot(finalRt, aes(x = NAME, y = type, fill = MEAN)) +
  geom_tile() +
  geom_text(aes(label = Significance)) +
  # scale_fill_gradientn(colors=Color) +
  scale_fill_distiller_custom(Color,breaks=Jbrs,limits = c(min(Jbrs),max(Jbrs))) +
  labs(x = "Monsoon Regions", y = "Climate Indices") +
  theme(
    text=element_text(family="Helvetica"),
    axis.title.x = element_text(family = "Helvetica",size = 12,face = "bold"),
    axis.title.y = element_text(family = "Helvetica", size = 12,face = "bold"),
    axis.text.x = element_text(family = "Helvetica",size = 8,face = "bold"),
    axis.text.y = element_text(family = "Helvetica",size = 8,face = "bold"),
    legend.title = element_blank(),
    # legend.position = "none",
    legend.text.align = 1,      
    plot.title = element_text(hjust = 0.5,family = "Helvetica",size = 12,face = "bold")
  )


ggsave("Climate Indices and HSWFE.png",
       plot=pp, width = 7.22, height = 6.22, dpi = 300,units="in")





