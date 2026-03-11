# Extract data from google earth engine S2 data 
# https://code.earthengine.google.com/?accept_repo=users/camilleminaudo/whiting

# Camille Minaudo
# January 2023

# ------------------- Load packages -------------------------


cat("\014")
rm(list = ls())

library(ggplot2)
library(plot3D)
library(stringi)
library(stringr)
library(egg)
library(lubridate)


# ------------------- Load data -------------------------

setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing")

df_sau <- read.table("gee_BGRarea_Sau.csv", header = T, sep = ";", dec = ".")
df_susq <- read.table("gee_BGRarea_Susqueda.csv", header = T, sep = ";", dec = ".")
df_gerg <- read.table("gee_BGRarea_Gergal.csv", header = T, sep = ";", dec = ".")


df_sau$date  <- seq.Date(as.Date("2017-04-01"),by = "month", length.out = length(df_sau$Band.name))
df_susq$date  <- seq.Date(as.Date("2017-04-01"),by = "month", length.out = length(df_susq$Band.name))
df_gerg$date  <- seq.Date(as.Date("2017-04-01"),by = "month", length.out = length(df_gerg$Band.name))

df_sau$site  <- "Sau"
df_susq$site  <- "Susqueda"
df_gerg$site  <- "Gergal"



df_all <- rbind(df_sau, df_susq, df_gerg)

ggplot(df_all[-which(df_all$site=="Gergal"),], aes(date, Band.value/1000, colour = site))+
  geom_abline(slope = 0, intercept = 13)+
  geom_path(size=1.2)+theme_article()+ylab("median BGR area/1000")+xlab("")+scale_x_date(date_breaks = "1 year")





df_sau_sencast <- read.table("sencast_sau_c2rcc.csv", header = T, sep = ",")
df_sau_sencast <- df_sau_sencast[df_sau_sencast$n_valid_pixels_sau>300,]
# df_sau_sencast <- df_sau_sencast[!is.na(df_sau_sencast$mean_BGR_area_sau),]

ggplot()+
  geom_abline(slope = 0, intercept = 13)+
  geom_path(data = df_sau, aes(date, Band.value/1000, colour = "GEE"))+
  geom_point(data = df_sau, aes(date, Band.value/1000, colour = "GEE"))+
  geom_path(data = df_sau_sencast, aes(as.Date(date), med_BGR_area_sau/1000, colour = "Sencast"))+
  geom_point(data = df_sau_sencast, aes(as.Date(date), med_BGR_area_sau/1000, colour = "Sencast"))+
  theme_article()+ylab("BGR area/1000")+xlab("")+scale_x_date(date_breaks = "1 year")+
  scale_colour_viridis_d(option = "A", begin = 0.1, end = 0.6)+
  theme(legend.title = element_blank(), legend.position = c(0.05,0.8))



ggplot()+
  geom_abline(slope = 0, intercept = 13)+
  geom_path(data = df_sau_sencast, aes(as.Date(date), med_BGR_area_sau/1000, colour = "Sau"))+
  geom_point(data = df_sau_sencast, aes(as.Date(date), med_BGR_area_sau/1000, colour = "Sau"))+
  geom_path(data = df_sau_sencast, aes(as.Date(date), med_BGR_area_susqueda/1000, colour = "Susqueda"))+
  geom_point(data = df_sau_sencast, aes(as.Date(date), med_BGR_area_susqueda/1000, colour = "Susqueda"))+
  theme_article()+ylab("median BGR area/1000")+xlab("")+scale_x_date(date_breaks = "1 week", limits = as.Date(c("2017-06-01", "2017-09-30")))+
  scale_colour_viridis_d(option = "A", begin = 0.1, end = 0.6)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  ylim(c(0,30))+
  theme(axis.text.x = element_text(angle = 90))




