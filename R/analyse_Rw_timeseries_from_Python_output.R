

cat("/014")
rm(list = ls())

library(ggplot2)
library(readr)
library(tidyr)
library(stringr)



datapath <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/CALCYOM/lakes_cci_database/Pyramid"
setwd(datapath)
myfile = "lake_reflectance_timeseries.csv"

data <- read.csv(myfile)
data$date <- as.Date(data$date)

data_Rw560 <- data[which(data$band=="Rw560"),]


ggplot(data_Rw560)+
  # geom_segment(aes(x = date, y=mean-std/2, yend = mean+std/2), size=1)+
  geom_segment(aes(x = date, y=p5, yend = p95), color = "grey", alpha=0.5, size=1)+
  # geom_point(aes(date, mean), size=1)+
  geom_point(aes(date, median), color = "black", alpha=0.5, size=1)+
  ylab("Rw560")+
  ggtitle("Pyramid Lake - 2017-2024")+
  theme_bw()

