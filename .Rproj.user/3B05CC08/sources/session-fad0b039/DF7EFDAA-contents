



cat("\014")
rm(list = ls())

library(ggplot2)
library(gridExtra)
library(grid)
library(lubridate)
library(dplyr)
library(tidyr)
library(zoo)
library(egg)
library(openair)
library(stringi)
library(stringr)


tabBGR_sau <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/BGR_Sau_Landsat4578_Sentinel2_forR.csv",
                     header = T, sep = ",", na.strings = "NA")
tabBGR_sau$Date <- as.Date(tabBGR_sau$Date, format = "%d/%m/%Y")
tabBGR_sau$system <- "sau"

tabBGR_sau_gath <- gather(tabBGR_sau, source, BGR, -Date, -system)
tabBGR_sau_gath$BGR <- as.numeric(tabBGR_sau_gath$BGR)

tabBGR_sau_gath$month <-  month(tabBGR_sau_gath$Date)
tabBGR_sau_gath$year <-  year(tabBGR_sau_gath$Date)


p_BGR_sau <- ggplot(tabBGR_sau_gath, aes(Date, BGR/1e+4, colour = source))+
  geom_abline(slope = 0, intercept = 1.3, alpha = 0.5)+
  geom_point()+geom_path()+
  theme_article()+
  xlab("")+
  ylab("BGR area*1e-4 [dl]")+theme(legend.position = c(0.05, 0.8), legend.title = element_blank())+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)

library(scales) 
p_inset2017 <- p_BGR_sau + scale_x_date(labels = date_format("%m"), breaks = "1 month", limits = as.Date(c("2017-01-01", "2017-12-31")))+
  theme(legend.position = 'none')+ylab("")

ggsave(filename = "BGR_sau_all_sources.png", plot = p_BGR_sau, path = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/", 
       scale = 0.5, width = 15, height = 8, units = 'in', dpi = 300)

ggsave(filename = "BGR_sau_all_sources_inset2017.png", plot = p_inset2017, path = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/", 
       scale = 0.5, width = 6, height = 4, units = 'in', dpi = 300)



p_BGR_sau <- p_BGR_sau + theme(legend.position = c(0.08, 0.8))
p_seas <- ggplot(tabBGR_sau_gath, aes(month, BGR*1e-4))+
  geom_abline(slope = 0, intercept = 1.3, alpha = 0.5)+
  geom_jitter(aes(colour = source))+
  geom_smooth(method = "loess", fill='grey', alpha=0.7, colour = "black")+
  theme_article()+theme(legend.position = 'none')+
  scale_x_continuous(breaks = seq(1,12))+
  ylab("")+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)

p_BGR_sau_trend_and_seas <- ggarrange(p_BGR_sau, p_seas, nrow = 1, widths = c(0.8,0.2))

ggsave(filename = "BGR_sau_all_sources_plus_seasonal.png", plot = p_BGR_sau_trend_and_seas, path = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/", 
       scale = 0.5, width = 20, height = 5, units = 'in', dpi = 300)




isF = T
for (y in unique(year(tabBGR_sau$Date))){
  slct <- tabBGR_sau[year(tabBGR_sau$Date) == y,]
  
  slct$meanBGR <- apply(as.matrix(slct[,seq(2,6)]), 1, mean, na.rm=T)
  slct <- slct[!is.na(slct$meanBGR),]
  
  df_n_whiting.temp <- data.frame(year = y,
                                  n = dim(slct)[1],
                                  n_above = dim(slct[slct$meanBGR>13000,])[1])
  
  df_n_whiting.temp$prop <- df_n_whiting.temp$n_above/df_n_whiting.temp$n*100
  if (isF){
    isF = F
    df_n_whiting <- df_n_whiting.temp
  } else {
    df_n_whiting <- rbind(df_n_whiting, df_n_whiting.temp)
  }
}
plot(df_n_whiting$year, df_n_whiting$prop, )
mean(df_n_whiting$prop, na.rm = T)

length(df_n_whiting$year[df_n_whiting$prop==0])



tabBGR_susqueda <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/BGR_susqueda_Landsat4578_Sentinel2_forR.csv",
                         header = T, sep = ",", na.strings = "NA")
tabBGR_susqueda$Date <- as.Date(tabBGR_susqueda$Date, format = "%d/%m/%Y")
tabBGR_susqueda$system <- "susqueda"

tabBGR_susqueda_gath <- gather(tabBGR_susqueda, source, BGR, -Date,-system)
tabBGR_susqueda_gath$BGR <- as.numeric(tabBGR_susqueda_gath$BGR)

tabBGR_susqueda_gath$year <-  year(tabBGR_susqueda_gath$Date)
tabBGR_susqueda_gath$month <-  month(tabBGR_susqueda_gath$Date)





ggplot()+
  geom_point(data = tabBGR_sau_gath, aes(Date, BGR, colour = "sau"))+
  geom_point(data = tabBGR_susqueda_gath, aes(Date, BGR, colour = "susqueda"))+
  theme_article()



ggplot()+
  geom_path(data = tabBGR_sau, aes(Date, Landsat7, colour = "sau"))+
  geom_path(data = tabBGR_susqueda, aes(Date, Landsat7, colour = "susqueda"))+
  theme_article()








spectra <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/rhown_spectra_sau_susqueda_20170622.csv",
                              header = T, sep = "\t", na.strings = "NA")

names(spectra) <- c("Wavelength", "Sau", "Susqueda")

spectra_gath <- gather(spectra, site, Rrs, -Wavelength)
p_spectra <- ggplot(spectra_gath, aes(Wavelength, Rrs, colour = site))+
  geom_path(size=1)+
  geom_point(size=2)+
  theme_article()+
  xlab("Wavelength [nm]")+
  ylab("Rrs [dl]")+theme(legend.position = c(0.12, 0.90), legend.title = element_blank())+
  scale_colour_viridis_d(option = "D", begin = 0.2, end = 0.8)


ggsave(filename = "spectra_Rrs_sau_Susqueda_20170622.png", plot = p_spectra, path = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/", 
       scale = 0.5, width = 10, height = 8, units = 'in', dpi = 300)


BRG_sau <- 0.5*abs(490 * spectra$Sau[3] + 560 * spectra$Sau[4] + 665 * spectra$Sau[2] - 560 * spectra$Sau[2] - 665 * spectra$Sau[3] - 490 * spectra$Sau[4])
BRG_susqueda <- 0.5*abs(490 * spectra$Susqueda[3] + 560 * spectra$Susqueda[4] + 665 * spectra$Susqueda[2] - 560 * spectra$Susqueda[2] - 665 * spectra$Susqueda[3] - 490 * spectra$Susqueda[4])


print(paste("BGR area in Sau = ",round(BRG_sau*100)/100,"; BGR area in Susqueda = ",round(BRG_susqueda*100)/100))









Rrs <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/remote_sensing/Rrs_sau_susqueda.csv",
                      header = T, sep = ",", na.strings = "NA")

Rrs$date <- as.Date(Rrs$date)
Rrs <- Rrs[order(Rrs$date),]
ggplot(Rrs, aes(date, ls_blue))+geom_path()+theme_article()

matplot(x = c(490, 560, 665, 860, 1600, 2200), y = t(Rrs[1:100, c("ls_blue", "ls_green", "ls_red", "ls_nir", "ls_swir1", "ls_swir2")])/1e+4, type = 'l')



Rrs_sau = Rrs[Rrs$plot_id<4, ]
Rrs_sau$BGRarea <- 0.5*1e-4*abs(490 * Rrs_sau$ls_green + 560 * Rrs_sau$ls_red + 665 * Rrs_sau$ls_blue - 560 * Rrs_sau$ls_blue - 665 * Rrs_sau$ls_green - 490 * Rrs_sau$ls_red)


Rrs_susqueda = Rrs[Rrs$plot_id>=4, ]
Rrs_susqueda$BGRarea <- 0.5*1e-4*abs(490 * Rrs_susqueda$ls_green + 560 * Rrs_susqueda$ls_red + 665 * Rrs_susqueda$ls_blue - 560 * Rrs_susqueda$ls_blue - 665 * Rrs_susqueda$ls_green - 490 * Rrs_susqueda$ls_red)


ggplot()+
  geom_path(data = Rrs_sau, aes(date, BGRarea, colour = satellite))+
  # geom_point(aes(colour = satellite), size = 0.5, alpha=0.2)+
  # geom_path(data = Rrs_susqueda, aes(date, BGRarea, colour = "susqueda"))+
  
  theme_article()




