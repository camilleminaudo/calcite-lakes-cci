





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


ggplot(tab, aes(n_dates_p95_BGR_extreme))+geom_density()+theme_bw()

ggplot(tab, aes(p95_BGR_mean, p95_BGR_std ))+geom_point()

ggplot(tab, aes(p95_BGR_median, p95_BGR_p75 - p95_BGR_p25 ))+geom_point()

ggplot(tab, aes(median_BGR_median, p95_BGR_p95))+geom_point()

list_whitings <- c("GLWD00000411","HYLA00014167","GLWD00000327","GLWD0000067","HYLA00000772","HYLA00000773","GLWD0000236","GLWD0000135")
tab.whitings <- tab[tab$lake_id %in% list_whitings,]
tab.whitings$lake_name

ggplot(mapping = aes(median_BGR_median, p95_BGR_p95))+
  geom_point(data = tab)+
  geom_point(data = tab.whitings, aes(colour = "known whitings"), size=4)+theme_bw()


ggplot(mapping = aes(median_BGR_median, n_dates_p95_BGR_extreme))+
  geom_point(data = tab)+
  geom_point(data = tab.whitings, aes(colour = "known whitings"), size=4)+theme_bw()



tab$file_name[which(tab$median_BGR_median<2 & tab$p95_BGR_p95>5 & tab$p95_BGR_p95>tab$BGR_mean_plus_3std & tab$n_dates_p95_BGR_extreme/tab$n_unique_dates>0.1)]


# ggplot(tab, aes(median_BGR_median, BGR_threshold_mean_plus_3std ))+geom_point()

tab <- tab[order(tab$p95_BGR_median),]
tab$order <- as.factor(seq(1,nrow(tab)))

ggplot(tab)+
  geom_segment(aes(order, y=p95_BGR_p25, yend = p95_BGR_p75 ))+
  geom_point(aes(order, p95_BGR_median))+
  xlab("")+
  ylab("BGR")+
  coord_flip()+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw()

