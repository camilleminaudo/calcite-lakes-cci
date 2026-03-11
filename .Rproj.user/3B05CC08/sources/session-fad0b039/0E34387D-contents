# Extract data from C2RCC S2 data after Sencast extraction

# Camille Minaudo
# August 2022

# ------------------- Load packages -------------------------


cat("\014")
rm(list = ls())

library(ncdf4)
library(ncdf4.helpers)
library(cmocean)
library(ggplot2)
library(plot3D)
library(stringi)
library(stringr)
library(egg)
library(lubridate)


# ------------------- Load data -------------------------

isF = T
data_paths <- c("D:/MariaZambrano_UB2022/DIAS/output_data/parameters_sau_S2_sau_2015-01-01_2018-01-01/L2C2RCC/",
                "D:/MariaZambrano_UB2022/DIAS/output_data/parameters_sau_S2_sau_2018-01-01_2020-01-01/L2C2RCC/",
                "D:/MariaZambrano_UB2022/DIAS/output_data/parameters_sau_S2_sau_2020-01-01_2022-09-15/L2C2RCC/")

for (data_path in data_paths){
  print(".........................")
  print(data_path)
  lof = list.files(data_path, pattern = ".nc")
  for(my_f in lof){
    print(my_f)
    df_temp <- data.frame(product = unlist(str_split(my_f, "_"))[1],
                          satellite = unlist(str_split(my_f, "_"))[3],
                          instrument = unlist(str_split(my_f, "_"))[4],
                          datestamp =  unlist(str_split(my_f, "_"))[5])
    date_stamp_split <- unlist(str_split(df_temp$datestamp, pattern = "T"))
    df_temp$date <- as.POSIXct(paste(date_stamp_split[1],date_stamp_split[2],sep = " "), format = "%Y%m%d %H%M%S", tz = 'UTC')
    
    if(unlist(str_split(my_f, "_"))[1] == "Mosaic"){
      df_temp <- data.frame(product = unlist(str_split(my_f, "_"))[2],
                            satellite = unlist(str_split(my_f, "_"))[4],
                            instrument = unlist(str_split(my_f, "_"))[5],
                            datestamp =  unlist(str_split(my_f, "_"))[6])
      date_stamp_split <- unlist(str_split(df_temp$datestamp, pattern = ".nc"))
      df_temp$date <- as.POSIXct(date_stamp_split[1], format = "%Y%m%d", tz = 'UTC')
    }
    
    my_sat <- unlist(str_split(my_f, "_"))[3]
    my_sat <- unlist(str_split(my_f, "_"))[3]
    my_sat <- unlist(str_split(my_f, "_"))[3]
    my_sat <- unlist(str_split(my_f, "_"))[3]
    
    # open netCDF file
    setwd(data_path)
    ncin <- nc_open(my_f)
    
    nc.get.variable.list(ncin, min.dims = 1)
    
    lon = ncvar_get(ncin, "lon")
    lat = ncvar_get(ncin, "lat")
    
    pxl_classif = ncvar_get(ncin, "pixel_classif_flags")
    
    extract_raster_water <- function(ncin, varname, pxl_classif){
      rast = ncvar_get(ncin, varname)
      
      rast_water <- rast
      rast_water[pxl_classif != 49280] <- NA
      rast_water[rast_water==0] <- NA
      # image2D(lon, sort(lat), z = rast_water[,order(lat)], col=cmocean("algae")(56))
      
      return(rast_water)
    }
    
    chl_water <- extract_raster_water(ncin = ncin, varname = "conc_chl", pxl_classif = pxl_classif)
    apig_water = extract_raster_water(ncin = ncin, varname = "iop_apig", pxl_classif = pxl_classif)
    bwit_water = extract_raster_water(ncin = ncin, varname = "iop_bwit", pxl_classif = pxl_classif)
    atot_water = extract_raster_water(ncin = ncin, varname = "iop_atot", pxl_classif = pxl_classif)
    kd489_water = extract_raster_water(ncin = ncin, varname = "kd489", pxl_classif = pxl_classif)
    tsm_water = extract_raster_water(ncin = ncin, varname = "conc_tsm", pxl_classif = pxl_classif)
    r490_water = extract_raster_water(ncin = ncin, varname = "rhow_B2", pxl_classif = pxl_classif)
    r560_water = extract_raster_water(ncin = ncin, varname = "rhow_B3", pxl_classif = pxl_classif)
    r665_water = extract_raster_water(ncin = ncin, varname = "rhow_B4", pxl_classif = pxl_classif)
    
    BGR_area = 10000 * 0.5 * abs(490 * r560_water + 560 * r665_water + 665 * r490_water - 
                                   560 * r490_water - 665 * r560_water - 490 * r665_water)
    # image2D(lon, sort(lat), z = BGR_area[,order(lat)], col=cmocean("delta")(56))
    
    
    get_my_stats <- function(mat, varname, df_temp){
      df_temp[[paste("mean",varname, sep = "_")]] <- mean(mat, na.rm = T)
      df_temp[[paste("p10",varname, sep = "_")]] <- quantile(mat, na.rm = T, probs = 0.1)
      df_temp[[paste("med",varname, sep = "_")]] <- median(mat, na.rm = T)
      df_temp[[paste("p90",varname, sep = "_")]] <- quantile(mat, na.rm = T, probs = 0.9)
      df_temp[[paste("sd",varname, sep = "_")]] <- sd(mat, na.rm = T)
      return(df_temp)
    }
    
    
    df_temp$n_valid_pixels_sau = sum(!is.na(chl_water[lon < 2.45,]))
    df_temp$n_valid_pixels_susqueda = sum(!is.na(chl_water[lon >= 2.50,]))
    
    df_temp <- get_my_stats(mat = chl_water[lon < 2.45,], varname = "chl_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = apig_water[lon < 2.45,], varname = "apig_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = bwit_water[lon < 2.45,], varname = "bwit_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = atot_water[lon < 2.45,], varname = "atot_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = kd489_water[lon < 2.45,], varname = "kd489_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = tsm_water[lon < 2.45,], varname = "tsm_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = r560_water[lon < 2.45,], varname = "r560_sau", df_temp = df_temp)
    df_temp <- get_my_stats(mat = BGR_area[lon < 2.45,], varname = "BGR_area_sau", df_temp = df_temp)
    
    df_temp <- get_my_stats(mat = chl_water[lon >= 2.50,], varname = "chl_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = apig_water[lon >= 2.50,], varname = "apig_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = bwit_water[lon >= 2.50,], varname = "bwit_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = atot_water[lon >= 2.50,], varname = "atot_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = kd489_water[lon >= 2.50,], varname = "kd489_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = tsm_water[lon >= 2.50,], varname = "tsm_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = r560_water[lon >= 2.50,], varname = "r560_susqueda", df_temp = df_temp)
    df_temp <- get_my_stats(mat = BGR_area[lon >= 2.50,], varname = "BGR_area_susqueda", df_temp = df_temp)
    
    if(isF){
      isF = F
      df_out <- df_temp
    } else {
      df_out <- rbind(df_out, df_temp)
    }
  }
  
}


df_out <- df_out[order(df_out$date),]
df_out$doy <- as.numeric(strftime(df_out$date, format = "%j"))

write.table(x = df_out, file = "C:/Projects/data/sau/remote_sensing/sencast_sau_c2rcc.csv", quote = F, sep = ",")


# ------------------- Plots -------------------------


npixels_sau <- ggplot(df_out, aes(date, n_valid_pixels_sau, colour = satellite))+geom_point()+theme_article()
npixels_susqueda <- ggplot(df_out, aes(date, n_valid_pixels_susqueda, colour = satellite))+geom_point()+theme_article()
ggarrange(npixels_sau, npixels_susqueda, ncol = 1)


df_out_sel_sau <- df_out[df_out$n_valid_pixels_sau>400,]
# df_out_sel_sau <- df_out_sel_sau[df_out_sel_sau$med_tsm_sau<150,]
df_out_sel_susqueda <- df_out[df_out$n_valid_pixels_susqueda>300,]
# df_out_sel_susqueda <- df_out_sel_susqueda[df_out_sel_susqueda$med_tsm_susqueda<150,]

tsm_sau <- ggplot(df_out_sel_sau, aes(date, med_tsm_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
tsm_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_tsm_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(tsm_sau, tsm_susqueda, ncol = 1)

chl_sau <- ggplot(df_out_sel_sau, aes(date, p90_chl_sau ))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
chl_susqueda <- ggplot(df_out_sel_susqueda, aes(date, p90_chl_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(chl_sau, chl_susqueda, ncol = 1)


apig_sau <- ggplot(df_out_sel_sau, aes(date, med_apig_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
apig_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_apig_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(apig_sau, apig_susqueda, ncol = 1)


bwit_sau <- ggplot(df_out_sel_sau, aes(date, med_bwit_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
bwit_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_bwit_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(bwit_sau, bwit_susqueda, ncol = 1)


atot_sau <- ggplot(df_out_sel_sau, aes(date, med_atot_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
atot_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_atot_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(atot_sau, atot_susqueda, ncol = 1)


kd489_sau <- ggplot(df_out_sel_sau, aes(date, med_kd489_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
kd489_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_kd489_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(kd489_sau, kd489_susqueda, ncol = 1)

r560_sau <- ggplot(df_out_sel_sau, aes(date, med_r560_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
r560_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_r560_susqueda))+geom_point(aes(colour = satellite))+geom_path()+theme_article()
ggarrange(r560_sau, r560_susqueda, ncol = 1)

BGR_area_sau <- ggplot(df_out_sel_sau, aes(date, med_BGR_area_sau))+geom_point(aes(colour = satellite))+geom_path()+theme_article()+geom_abline(slope = 0, intercept = 13000)
BGR_area_susqueda <- ggplot(df_out_sel_susqueda, aes(date, med_BGR_area_susqueda ))+geom_point(aes(colour = satellite))+geom_path()+theme_article()+geom_abline(slope = 0, intercept = 13000)
ggarrange(BGR_area_sau, BGR_area_susqueda, ncol = 1)


df_out_sel_sau[which(df_out_sel_sau$med_BGR_area_sau > 20000),]





setwd("C:/Projects/data/Sau/Historical_data/")

data_chla <- read.table(file = "Chl-a.csv", header = T, sep = ",")
data_chla$date = as.POSIXct(data_chla$Fecha)

data_chla_0 <- data_chla[data_chla$Depth == 0,]


chl_sau + geom_point(data = data_chla_0[year(data_chla_0$date)>2014,], aes(date, Valor))


df_matchups <- df_out_sel_sau
df_matchups$inSitu_chla <- NA
for(ti in seq(1, dim(df_matchups)[1])){
  
  ind_closest <- which.min(abs(as.numeric(as.Date(df_matchups$date[ti])) - as.numeric(as.Date(data_chla_0$date))))
  my_diff_time <- difftime(as.Date(df_matchups$date[ti]), as.Date(data_chla_0$date[ind_closest]), units = "days")
  
  if(my_diff_time <= 2){
    df_matchups$inSitu_chla[ti] <- data_chla_0$Valor[ind_closest]
  }
  
}


ggplot(df_matchups)+
  geom_point(aes(date, inSitu_chla, colour = "in situ"))+
  geom_point(aes(date, mean_chl_sau, colour = "S2MSI"))+
  theme_article()+scale_colour_viridis_d()


ggplot(df_matchups)+
  geom_point(aes(inSitu_chla, mean_chl_sau))+
  theme_article()+scale_colour_viridis_d()




