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


# ------------------- Load data -------------------------

data_path <- "D:/MariaZambrano_UB2022/DIAS/output_data/parameters_pyramid_S2_pyramid_2020-01-01_2020-12-31/L2C2RCC/"
lof = list.files(data_path, pattern = ".nc")


mask <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/other_projects/MSCA_fellowship/remoteSensing/IDEPIX_WATER_Mask.txt", 
                   skip = 4, header = T)

plot(mask$Longitude, mask$Latitude)




my_f = lof[300]
setwd(data_path)
ncin <- nc_open(my_f)
nc.get.variable.list(ncin, min.dims = 1)
lon = ncvar_get(ncin, "lon")
lat = ncvar_get(ncin, "lat")
pxl_classif = ncvar_get(ncin, "pixel_classif_flags")
mypxl_water <- pxl_classif
mypxl_water[pxl_classif < 3e+4 | pxl_classif > 6e+4] <- NA
image2D(lon, sort(lat), z = mypxl_water[,order(lat)])

myMask <- !is.na(mypxl_water) # determined with lof[300], "Mosaic_L2C2RCC_NA_S2A_MSIL1C_20200204.nc"

nc_close(ncin)


isF = T
for(my_f in lof){
  print(my_f)
  df_temp <- data.frame(product = unlist(str_split(my_f, "_"))[1],
                        satellite = unlist(str_split(my_f, "_"))[3],
                        instrument = unlist(str_split(my_f, "_"))[4],
                        datestamp =  unlist(str_split(my_f, "_"))[5])
  date_stamp_split <- unlist(str_split(df_temp$datestamp, pattern = "T"))
  df_temp$date <- as.POSIXct(paste(date_stamp_split[1],date_stamp_split[2],sep = " "), format = "%Y%m%d %H%M%S", tz = 'UTC')
  
  # open netCDF file
  setwd(data_path)
  ncin <- nc_open(my_f)
  
  nc.get.variable.list(ncin, min.dims = 1)
  
  lon = ncvar_get(ncin, "lon")
  lat = ncvar_get(ncin, "lat")
  
  # mat_mask <- matrix(data = F, nrow = length(lon), ncol = length(lat))
  # for(i in dim(mat_mask)[1]){
  #   for (j in dim(mat_mask)[2]){
  #     lon[i]
  #     lat[j]
  #     nearest_valid <- which.min(abs())
  #     mat_mask[i,j]
  #   }
  # }
  
  
  pxl_classif = ncvar_get(ncin, "pixel_classif_flags")
  
  flags_content <- ncatt_get(ncin,"pixel_classif_flags")
  flag_values <- flags_content$flag_masks
  flag_descrpt <- str_split(flags_content$flag_descriptions, pattern = "\t")
  
  code_water <- flag_values[which(as.vector(unlist(flag_descrpt)) == "Water pixels" )]
  
  pxl_water <- pxl_classif
  pxl_water[pxl_classif!=code_water] <- NA
  # pxl_classif[!myMask] <- NA
  image2D(lon, sort(lat), z = pxl_water[,order(lat)])

  chl = ncvar_get(ncin, "unc_chl")
  apig = ncvar_get(ncin, "iop_apig")
  bwit = ncvar_get(ncin, "unc_bwit")
  r560 = ncvar_get(ncin, "rhow_B3")
  rtoa490 = ncvar_get(ncin, "rtoa_B2")
  rtoa560 = ncvar_get(ncin, "rtoa_B3")
  rtoa665 = ncvar_get(ncin, "rtoa_B4")
  
  
  nc_close(ncin)
  
  
  chl_water <- chl
  # chl_water[pxl_classif != 49280] <- NA
  chl_water[chl_water==0] <- NA
  # image2D(lon, sort(lat), z = chl_water[,order(lat)], col=cmocean("algae")(56))
  
  apig_water <- apig
  # apig_water[pxl_classif != 49280] <- NA
  apig_water[apig_water==0] <- NA
  # image(lon, sort(lat), apig_water[,order(lat)], col=cmocean("delta")(56))
  
  bwit_water <- bwit
  # bwit_water[pxl_classif != 49280] <- NA
  bwit_water[bwit_water==0] <- NA
  # image(lon, sort(lat), bwit_water[,order(lat)], col=cmocean("curl")(56), zlim = c(0,0.1))
  
  # ggplot(data = data.frame(bwit = as.vector(apig_water)), aes(bwit))+geom_density()
  
  r560_water <- r560
  r560_water[!myMask] <- NA
  r560_water[r560_water==0] <- NA
  # image2D(lon, sort(lat), z = r560_water[,order(lat)], col=cmocean("delta")(56))
  
  
  rtoa490_water <- rtoa490
  rtoa490_water[!myMask] <- NA
  rtoa490_water[rtoa490_water==0] <- NA
  
  rtoa560_water <- rtoa560
  rtoa560_water[as.matrix(!myMask)] <- NA
  # rtoa560_water[rtoa560_water==0] <- NA
  # image2D(lon, sort(lat), z = rtoa560_water[,order(lat)], col=cmocean("delta")(56))
  
  rtoa665_water <- rtoa665
  rtoa665_water[!myMask] <- NA
  rtoa665_water[rtoa665_water==0] <- NA
  
  # image2D(lon, sort(lat), z = rtoa560_water[,order(lat)], col=cmocean("delta")(56))
  
  BGR_area = 10000 * 0.5 * abs(490 * rtoa560_water + 560 * rtoa665_water + 665 * rtoa490_water - 
                                 560 * rtoa490_water - 665 * rtoa560_water - 490 * rtoa665_water)
  # image2D(lon, sort(lat), z = BGR_area[,order(lat)], col=cmocean("delta")(56))
  
  
  df_temp$chl_avg = mean(chl_water, na.rm = T)
  df_temp$chl_med = median(chl_water, na.rm = T)
  # df_temp$n_valid_pixels = sum(!is.na(chl_water))
  
  df_temp$apig_avg = mean(apig_water, na.rm = T)
  df_temp$apig_med = median(apig_water, na.rm = T)

  
  df_temp$bwit_avg = mean(bwit_water, na.rm = T)
  df_temp$bwit_med = median(bwit_water, na.rm = T)

  
  df_temp$r560_avg = mean(r560_water, na.rm = T)
  df_temp$r560_med = median(r560_water, na.rm = T)
  
  
  df_temp$rtoa560_avg = mean(rtoa560_water, na.rm = T)
  df_temp$rtoa560_med = median(rtoa560_water, na.rm = T)
  df_temp$rtoa560_sd = sd(rtoa560_water, na.rm = T)
  # df_temp$n_valid_pixels = sum(!is.na(mypxl_water))
  
  df_temp$rtoa560_pin1_med = median(rtoa560_water[lon > -119.55 & lon < -119.54, lat > 39.95 & lat < 39.96], na.rm = T)
  df_temp$rtoa560_pin2_med = median(rtoa560_water[lon > -119.57 & lon < -119.56, lat > 40.12 & lat < 40.13], na.rm = T)
  
  df_temp$BGR_area_avg = mean(BGR_area, na.rm = T)
  df_temp$BGR_area_med = median(BGR_area, na.rm = T)
  df_temp$BGR_area_sd = sd(BGR_area, na.rm = T)

  
  if(isF){
    isF = F
    df_out <- df_temp
  } else {
    df_out <- rbind(df_out, df_temp)
  }
}


# bwit <- ggplot(df_out, aes(date, bwit_med, colour = satellite))+geom_point()+theme_article()
# 
# chl <- ggplot(df_out, aes(date, chl_med, colour = satellite))+geom_point()+theme_article()
# 
# apig <- ggplot(df_out, aes(date, apig_med, colour = satellite))+geom_point()+theme_article()
# 
# r560 <- ggplot(df_out, aes(date, r560_med, colour = satellite))+geom_point()+theme_article()

rtoa560 <- ggplot(df_out, aes(date, rtoa560_med))+
  geom_smooth(method = "loess", colour = "white", span=0.2)+
  geom_point()+
  geom_point(aes(date, rtoa560_pin1_med, colour = "pin1"))+
  geom_point(aes(date, rtoa560_pin2_med, colour = "pin2"))+
  theme_article()+ylim(c(0,0.3))+
  ylab("Rrs at 560 nm (dl)")+
  xlab("")+
  geom_vline(xintercept = as.POSIXct("2020-08-09"), color = 'red')
rtoa560


BGR_area <- ggplot(df_out, aes(date, BGR_area_med))+geom_point()+
  # geom_ribbon(aes(x = date, ymin = BGR_area_med-BGR_area_sd, ymax = BGR_area_med+BGR_area_sd), alpha = 0.2)+
  # geom_smooth(method = "loess", colour = "white", )+
  theme_article()+ylim(c(0,NA))+
  ylab("BGR area (dl)")+
  xlab("")
BGR_area

ggplot(df_out, aes(date, n_valid_pixels, colour = satellite))+
  geom_point()+
  theme_article()
