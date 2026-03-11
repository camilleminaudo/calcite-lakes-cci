

library(sf)
library(readr)


setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
lakescci <- sf::read_sf("data/lakes_cci/lakescci_v2.1.0_alkalinity_rivernet2026.shp")
dim(lakescci)

st_crs(lakescci)

ggplot(lakescci, aes(alk_mean/1000))+geom_density()+theme_bw()


# check which lakes are the ones exceeding median value overall
median_alk <- median(lakescci$alk_mean, na.rm = T)

ggplot(lakescci, aes(alk_mean, fill = type))+
  geom_density(alpha=0.5)+
  geom_vline(xintercept = median_alk)+
  theme_bw()

lakescci_highestAlk <- lakescci[which(lakescci$alk_mean>median_alk),]

setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
sf::write_sf(lakescci_highestAlk, dsn = "PROJECTS/CALCYOM/lakes_cci_database/lakes_cci_highAlk.shp", append = F)



