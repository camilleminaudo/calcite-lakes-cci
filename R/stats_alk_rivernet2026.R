

library(sf)
library(readr)
library(ggplot2)

setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
lakescci <- sf::read_sf("data/lakes_cci/lakescci_v2.1.0_alkalinity_rivernet2026.shp")
dim(lakescci)

st_crs(lakescci)

ggplot(lakescci, aes(alk_mean/1000))+geom_density()+theme_bw()


# check which lakes are the ones exceeding median value overall
median_alk <- median(lakescci$alk_mean, na.rm = T)

ggplot(lakescci, aes(alk_mean/1000, fill = type))+
  geom_density(alpha=0.5)+
  geom_vline(xintercept = median_alk/1000)+
  theme_bw()

lakescci_highestAlk <- lakescci[which(lakescci$alk_mean>median_alk),]

setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
sf::write_sf(lakescci_highestAlk, dsn = "PROJECTS/CALCYOM/lakes_cci_database/lakes_cci_highAlk.shp", append = F)


dim(lakescci_highestAlk)

head(lakescci_highestAlk$short_name)

library("mclm")
setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
write_txt(x = unique(lakescci_highestAlk$short_name), file = "PROJECTS/CALCYOM/lakes_cci_database/list_lakes_high_alkalinity.txt", line_glue = "\n")





# order lake names in lakescci by levels of alkalinity
list_sorted <- lakescci$short_name[order(lakescci$alk_mean, decreasing = T)]
library("mclm")
setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
write_txt(x = unique(list_sorted), file = "PROJECTS/CALCYOM/lakes_cci_database/list_lakes_decreasing_alk.txt", line_glue = "\n")


# Parameters
chunk_size <- 200
file_prefix <- "list_lakes_decreasing_alk"

# Split into chunks
chunks <- split(list_sorted, ceiling(seq_along(list_sorted) / chunk_size))

# Write each chunk to a file
setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/CALCYOM/lakes_cci_database/")
for (i in seq_along(chunks)) {
  file_name <- paste0(file_prefix, "_", i, ".txt")
  writeLines(chunks[[i]], file_name)
}



