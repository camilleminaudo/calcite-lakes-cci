





cat("/014")
rm(list = ls())

library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggpubr)


datapath <- "C:/Projects/myGit/calcite-lakes-cci/results/"
setwd(datapath)
tab <- read.csv("data_summary_L2_OLCI.csv")

head(tab)


ggplot(tab, aes(median_BGR_median ))+geom_boxplot()+theme_bw()

ggplot(tab, aes(median_BGR_median, n_dates_median_BGR_extreme ))+geom_point()

ggplot(tab, aes(median_BGR_median, BGR_threshold_mean_plus_3std ))+geom_point()


ggplot(tab)+
  geom_segment(aes(lake_name, y=median_BGR_p5, yend = median_BGR_p95 ))+
  geom_point(aes(lake_name, median_BGR_median))+
  xlab("")+
  ylab("BGR")+
  coord_flip()+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw()

