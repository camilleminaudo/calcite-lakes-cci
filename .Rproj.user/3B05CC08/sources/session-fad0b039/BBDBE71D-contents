



library(sf)
library(readr)


setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/")
lakescci <- sf::read_sf("data/lakes_cci/lakes_cci_v2.1.0_shp/shapefile/lakescci_v2.1.0_data-availability.shp")
dim(lakescci)

st_crs(lakescci)


list_sel <- read.csv(file = "PROJECTS/CALCYOM/lakes_cci_database/20260309_Xiaohan_candidate_lakes.txt")
list_sel <- as.character(list_sel[[1]])

lakescci_sel <- lakescci[which(lakescci$short_name %in% list_sel),]

dim(lakescci_sel)

lakescci_sel$name

sf::write_sf(lakescci_sel, dsn = "PROJECTS/CALCYOM/lakes_cci_database/lakes_cci_preselected.shp", append = F)
