


cat("\014")
rm(list = ls())

library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(ggpubr)



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


p_TOA <- ggplot(data_all[which(data_all$isTOA==T),], aes(Wavelength, Rw))+geom_path(aes(group = uniqID, colour = lake))+theme_bw()+
  geom_hline(yintercept = 0)+
  xlab("Wavelength [nm]")+
  ylab("Oa [mW m-2 sr-1 nm-1]")+facet_wrap(.~lake, scales = "free_y")+scale_colour_viridis_d(option = "B", end = 0.8)+
  theme(legend.position = "none")


p_Rw <- ggplot(data_all[which(data_all$isTOA==F),], aes(Wavelength, Rw))+geom_path(aes(group = uniqID, colour = lake))+theme_bw()+
  geom_hline(yintercept = 0)+
  xlab("Wavelength [nm]")+
  # ylim(c(-.05,NA))+
  ylab("Rw [dl]")+facet_wrap(.~lake, scales = "free_y")+scale_colour_viridis_d(option = "B", end = 0.8)+
  theme(legend.position = "none")


p <- ggarrange(p_TOA, p_Rw, ncol = 1, align = "v")

plt.name <- "spectral_shapes_whitings.jpg"

ggsave(filename = plt.name, plot = p, path = "C:/Projects/myGit/calcite-lakes-cci/results",
       width = 8, height = 6, units = "in", dpi = 300)



# ------------------ spectral shape Pyramid 2020 ----------------

# dissociate bright pixels from the others
Rw_pyramid <- data_all[which(data_all$isTOA==F & data_all$lake=="Pyramid"),]
Rw_pyramid$uniqPins <- paste(Rw_pyramid$date,Rw_pyramid$pin, sep = "_")


list_suspicious <- unique(Rw_pyramid$uniqPins[which(Rw_pyramid$Wavelength> 900 & Rw_pyramid$Rw>0.17)])
Rw_pyramid <- Rw_pyramid[!Rw_pyramid$uniqPins %in% list_suspicious,]



list_bright <- unique(Rw_pyramid$uniqPins[which(Rw_pyramid$Rw>0.4)])
Rw_pyramid$bright <- F
Rw_pyramid$bright[Rw_pyramid$uniqPins %in% list_bright] <- T


Rw_pyramid$bright <- F
Rw_pyramid$bright[Rw_pyramid$uniqPins %in% list_bright] <- T

p <- ggplot(Rw_pyramid, aes(Wavelength, Rw))+
  geom_path(aes(group = uniqID, colour = bright), size=0.5, alpha=0.7)+theme_bw()+
  geom_hline(yintercept = 0)+
  xlab("Wavelength [nm]")+
  # ylim(c(-.05,NA))+
  ylab("Rw [dl]")+scale_colour_viridis_d(option = "B", end = 0.8)

ggsave(filename = "specxtral_shape_Lake_Pyramid.jpg", plot = p, path = "C:/Projects/myGit/calcite-lakes-cci/results",
       width = 6, height = 4, units = "in", dpi = 300)


data_sprd <- spread(Rw_pyramid, key = Wavelength, value = Rw)
data_sprd$date <- as.Date(data_sprd$date, format = "%Y%m%d")

data_sprd$G2R_ratio <- data_sprd[["560"]]/median(data_sprd[["560"]], na.rm = T)

ggplot(data_sprd, aes(date, G2R_ratio, colour = bright))+geom_point()+theme_bw()


# ------------------ temporal dynamics  ----------------

for (isTOA in c(TRUE,FALSE)){

  data_subset <- data_all[which(data_all$isTOA==isTOA),]

  data_sprd <- spread(data_subset, key = Wavelength, value = Rw)
  data_sprd$date <- as.Date(data_sprd$date, format = "%Y%m%d")


  selected_columns <- seq(which(names(data_sprd)=="400"),which(names(data_sprd)=="1020"))
  data_sprd$total_radiance <- apply(data_sprd[, selected_columns, drop = FALSE], 1, sum)

  myvars = c("560","total_radiance")

  for (myvar in myvars){
    lakes <- unique(data_sprd$lake)

    for (lake in lakes){
      data_l <- data_sprd[data_sprd$lake==lake,]
      data_l_d <- NULL
      for(d in unique(data_l$date)){
        data_l_d.temp <- data.frame(date = as.Date(d),
                                    avg = mean(data_l[[myvar]][which(data_l$date==d)], na.rm = T),
                                    sd = sd(data_l[[myvar]][which(data_l$date==d)], na.rm = T),
                                    median = median(data_l[[myvar]][which(data_l$date==d)], na.rm = T),
                                    p5 = quantile(data_l[[myvar]][which(data_l$date==d)], 0.05, na.rm = T),
                                    p10 = quantile(data_l[[myvar]][which(data_l$date==d)], 0.1, na.rm = T),
                                    p90 = quantile(data_l[[myvar]][which(data_l$date==d)], 0.9, na.rm = T),
                                    p95 = quantile(data_l[[myvar]][which(data_l$date==d)], 0.95, na.rm = T)
        )
        data_l_d <- rbind(data_l_d, data_l_d.temp)
      }

      if (isTOA){
        vartype = "TOA"
        units = "[mW m-2 sr-1 nm-1]"
      } else {
        vartype = "Rw"
        units = "[dl]"
      }


      p_percs <- ggplot(data_l_d)+
        # geom_segment(aes(x = date, y=avg-sd/2, yend = avg+sd/2), size=2)+
        geom_segment(aes(x = date, y=p5, yend = p95), color = "grey10", alpha=0.5, size=2)+
        # geom_point(aes(date, avg), size=4)+
        geom_point(aes(date, median), color = "red", alpha=0.5, size=4)+
        ylab(paste(vartype, myvars, units, sep = " "))+
        ggtitle(paste(lake,vartype, myvars, "over time",sep = " "))+
        theme_bw()

      plt.name <- paste0("spectral_shapes_whitings_temporal_",lake,"_",vartype,"_",myvar,".jpg")

      ggsave(filename = plt.name, plot = p_percs, path = "C:/Projects/myGit/calcite-lakes-cci/results",
             width = 8, height = 4, units = "in", dpi = 300)

    }
  }

}




# --------------- WATER LEAVING REFLECTANCE ---------------------

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
  ylab("Rw560 [dl]")+
  xlab("BGR area [dl]")+
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

