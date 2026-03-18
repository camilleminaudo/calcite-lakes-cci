

cat("/014")
rm(list = ls())

library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggpubr)


datapath <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/CALCYOM/lakes_cci_database/extracted_Rw_time_series/"


# ------------- Plot Pyramid Lake reprocessed ------------------

myfile <- "Rw_OLCI_L2_Pyramid_GLWD00000411_20160425_20241231.csv"
filename_split <- str_split(myfile, pattern = "_")[[1]]

sensor = filename_split[2]
datalevel = filename_split[3]
lakename = filename_split[4]
lakeID = filename_split[5]

data <- read.csv(myfile)
data$date <- as.Date(data$date)

year_start = min(year(data$date))
year_end = max(year(data$date))

myvar = "Rw560"
longterm_med <- median(data$median[which(data$band==myvar)], na.rm = T)
ggplot(data[which(data$band==myvar),])+
  geom_point(aes(date, p90/longterm_med), color = "black", alpha=0.5, size=3)+
  xlab("")+
  ggtitle("OLCI L2 - Pyramid lake, Nevada, USA")+
  theme_bw()






p_Rw560 <- ggplot(data[which(data$band=="Rw560"),])+
  # geom_segment(aes(x = date, y=p5, yend = p95), color = "grey", alpha=0.5, size=1)+
  geom_point(aes(date, median), color = "black", alpha=0.5, size=1)+
  xlab("")+
  ylab("Rw at 560 nm [dl]")+
  ggtitle("OLCI L2 - Pyramid lake, Nevada, USA")+
  theme_bw()

ggplot(data[which(data$band=="Rw490"),])+
  geom_segment(aes(x = date, y=p5, yend = p95), color = "grey", alpha=0.5, size=1)+
  geom_point(aes(date, median), color = "black", alpha=0.5, size=1)+
  xlab("")+
  ylab("Rw at 560 nm [dl]")+
  ggtitle("OLCI L2 - Pyramid lake, Nevada, USA")+
  theme_bw()


p_BGR <- ggplot(data[which(data$band=="BGR"),])+
  geom_segment(aes(x = date, y=p5, yend = p95), color = "grey", alpha=0.5, size=1)+
  geom_point(aes(date, median), color = "black", alpha=0.5, size=1)+
  xlab("")+
  ylab("BGR area [dl]")+
  ggtitle("OLCI L2 - Pyramid lake, Nevada, USA")+
  theme_bw()


data_sel <- data[!duplicated(data$date),]
data_sel$bright_fraction <- data_sel$n_bright_pixels/data_sel$count*100

p_bright <- ggplot(data_sel[data_sel$count>30,])+
  geom_area(aes(date, bright_fraction), fill = "black",color = "black", alpha=0.5, size=1)+
  xlab("")+
  ylab("Bright fraction [%]")+
  ylim(c(0,100))+
  theme_bw()

plt <- ggarrange(p_Rw560, p_BGR, p_bright, ncol = 1, align = "v")
plt

plt.name = paste0("OLCI_L2_Pyramid_Lake",".jpg")

ggsave(filename = plt.name, plot = plt, path = "C:/Projects/myGit/calcite-lakes-cci/results",
       width = 8, height = 7, units = "in", dpi = 300)








# ------------- Plot them all ------------------

setwd(datapath)

myfiles <- list.files(path = datapath, pattern = ".csv")

for (myfile in myfiles){
  print(myfile)

  filename_split <- str_split(myfile, pattern = "_")[[1]]

  datalevel = filename_split[2]
  lakename = filename_split[3]
  lakeID = filename_split[4]



  data <- read.csv(myfile)
  data$date <- as.Date(data$date)

  year_start = min(year(data$date))
  year_end = max(year(data$date))



  myvars <- c("BGR","Rw560")

  for (myvar in myvars){

    data_sel <- data[which(data$band==myvar),]

    # require a minimum of 30 observations to trust statistics
    data_sel <- data_sel[data_sel$count>30,]

    mytitle <- paste(lakename,datalevel,year_start,"to",year_end, myvar,sep = "_")

    p <- ggplot(data_sel)+
      # geom_segment(aes(x = date, y=mean-std/2, yend = mean+std/2), size=1)+
      geom_segment(aes(x = date, y=p5, yend = p95), color = "grey", alpha=0.5, size=1)+
      # geom_point(aes(date, mean), size=1)+
      geom_point(aes(date, median), color = "black", alpha=0.5, size=1)+
      xlab("")+
      ylab(myvar)+
      ggtitle(mytitle)+
      theme_bw()

    plt.name = paste0(mytitle,".jpg")

    ggsave(filename = plt.name, plot = p, path = "C:/Projects/myGit/calcite-lakes-cci/results",
           width = 8, height = 3, units = "in", dpi = 300)
  }
}


# ------------- L2 vs L3 analysis ------------------

lakes <- c("PYRAMID","BOURGET","GENEVA","KIVU")

for (lake in lakes){
  print(lake)
  myfiles <- list.files(path = datapath, pattern = paste0("_",lake))

  data_all <- NULL
  for (myfile in myfiles){
    print(myfile)

    filename_split <- str_split(myfile, pattern = "_")[[1]]

    datalevel = filename_split[2]
    lakename = filename_split[3]
    lakeID = filename_split[4]

    data <- read.csv(myfile)
    data$date <- as.Date(data$date)
    data$level = as.factor(datalevel)

    data_all <- rbind(data_all, data)
  }

  date_start <- as.Date("2020-01-01")
  date_end <- as.Date("2020-12-31")


  myvars <- c("BGR","Rw560")

  for (myvar in myvars){

    data_sel <- data_all[which(data_all$band==myvar),]

    mytitle <- paste(lakename,date_start,"to",date_end, myvar,sep = "_")

    p1 <- ggplot(data_sel[which(data_sel$date>=date_start & data_sel$date<=date_end),])+
      # geom_point(aes(date, p95, colour = level), alpha=0.5, size=2)+
      geom_line(aes(date, p95, colour = level), alpha=1, size=0.5)+
      xlab("")+
      ylab(paste0("p95 of ",myvar))+
      ggtitle(mytitle)+
      theme_bw()+
      theme(legend.position = c(0.1, 0.6))

    p2 <- ggplot(data_sel[which(data_sel$date>=date_start & data_sel$date<=date_end & data_sel$level=="L3"),], aes(date, count))+
      # geom_point(aes(date, count, colour = level), alpha=0.5, size=2)+
      geom_line(alpha=1, size=0.5)+
      xlab("")+
      ylab("Count")+
      ggtitle("L3 valid pixel count")+
      theme_bw()+
      theme(legend.position = "none")

    p12 <- ggarrange(p1, p2, ncol = 1, heights = c(0.7,0.3), align = "v")

    plt.name = paste0("L2_vs_L3_",mytitle,".jpg")

    ggsave(filename = plt.name, plot = p12, path = "C:/Projects/myGit/calcite-lakes-cci/results",
           width = 8, height = 6, units = "in", dpi = 300)

  }
}
