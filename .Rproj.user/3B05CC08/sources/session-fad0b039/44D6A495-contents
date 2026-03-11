


cat("\014")
rm(list = ls())

library(ggplot2)
library(readr)
library(tidyr)
library(stringr)


get_spectral_shapes <- function(myfile, plotit = F){

  myl <- read_lines(file = myfile)
  data <- read.csv(myfile,
                   sep = "\t", header = T, na.strings = "NaN", nrows = length(myl)-2)

  isTOA <- str_detect(basename(myfile), "Oa")

  lake_name <- str_split(basename(myfile), pattern = "_")[[1]][2]
  mydate <- gsub(pattern = ".csv",x = str_split(basename(myfile), pattern = "_")[[1]][3], replacement = "")

  data_gath <- gather(data, pin, Rw, -Wavelength)

  if(plotit){
    p <- ggplot(data_gath[!is.na(data_gath$Rw),], aes(Wavelength, Rw))+geom_path(aes(group = pin))+theme_bw()+
      geom_hline(yintercept = 0)+
      xlab("Wavelength [nm]")

    if(isTOA){
      p <- p +ylab("Rw [mW m-2 sr-1 nm-1]")+ggtitle(paste0(lake_name," ",mydate," - Top of atmosphere"))
    } else {
      p <- p + ylab("Rw [dl]")+ggtitle(paste0(lake_name," ",mydate," - Water leaving"))
    }
    print(p)
  }


  data_gath$lake <- lake_name
  data_gath$date <- mydate
  data_gath$isTOA <- isTOA

  return(data_gath)
}



datapath <- "C:/Projects/myGit/calcite-lakes-cci/data/"

fs <- list.files(path = datapath, pattern = ".csv", full.names = T)

data_all <- NULL
for(f in fs){
  data_l <- get_spectral_shapes(myfile = f, plotit = F)
  data_all <- rbind(data_all, data_l)
}

data_all$uniqID <- paste(data_all$lake,data_all$date,data_all$pin, data_all$isTOA, sep = "_")

ggplot(data_all, aes(Wavelength, Rw))+geom_path(aes(group = uniqID, colour = lake))+theme_bw()+
  geom_hline(yintercept = 0)+
  xlab("Wavelength [nm]")+facet_grid(isTOA~lake, scales = "free_y")



data_Rw <- data_all[which(data_all$isTOA==F),]

data_sprd <- spread(data_Rw, key = Wavelength, value = Rw)
data_sprd$date <- as.Date(data_sprd$date, format = "%Y%m%d")


data_sprd$BGRarea <- 0.5 *
  abs(490*data_sprd[["560"]] +
        560*data_sprd[["665"]] +
        665*data_sprd[["490"]] -
        560*data_sprd[["490"]] -
        665*data_sprd[["560"]] -
        490*data_sprd[["665"]]
  )

ggplot(data_sprd, aes(BGRarea, fill = lake))+geom_density(alpha=0.5)+theme_bw()+
  facet_wrap(.~lake, scales = "free")


ggplot(data_sprd, aes(BGRarea, `560`))+geom_point()+theme_bw()+
  facet_wrap(.~lake, scales = "free")


data_Pyramid <- data_sprd[data_sprd$lake=="Pyramid",]
data_Pyramid_d <- NULL
for(d in unique(data_Pyramid$date)){
  data_Pyramid_d.temp <- data.frame(date = as.Date(d),
                                    avg_BGR = mean(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], na.rm = T),
                                    sd_BGR = sd(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], na.rm = T),
                                    median_BGR = median(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], na.rm = T),
                                    p5_BGR = quantile(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], 0.05, na.rm = T),
                                    p10_BGR = quantile(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], 0.1, na.rm = T),
                                    p90_BGR = quantile(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], 0.9, na.rm = T),
                                    p95_BGR = quantile(data_Pyramid$BGRarea[which(data_Pyramid$date==d)], 0.95, na.rm = T)
  )
  data_Pyramid_d <- rbind(data_Pyramid_d, data_Pyramid_d.temp)
}



ggplot(data_Pyramid_d)+
  geom_segment(aes(x = date, y=avg_BGR-sd_BGR/2, yend = avg_BGR+sd_BGR/2), size=2)+
  geom_segment(aes(x = date, y=p5_BGR, yend = p95_BGR), color = "red", alpha=0.5, size=2)+
  geom_point(aes(date, avg_BGR), size=4)+
  geom_point(aes(date, median_BGR), color = "red", alpha=0.5, size=4)+
  ylab("BGR triangle area")+
  ggtitle("Pyramid Lake - Whiting 2020")+
  theme_bw()

