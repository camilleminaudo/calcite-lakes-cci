

cat("/014")
rm(list = ls())

library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggpubr)
library(dplyr)


detect_whiting_events <- function(df, value_col = "p90",
                                  min_duration = 7,
                                  z_thresh = 5,
                                  max_gap_days = 15) {

  # Ensure sorted
  df <- df %>% arrange(date)

  values <- df[[value_col]]

  # ----------------------------
  # Robust stats
  # ----------------------------
  med <- median(values, na.rm = TRUE)
  mad_val <- mad(values, constant = 1, na.rm = TRUE)
  if (mad_val == 0) mad_val <- 1e-6

  df <- df %>%
    mutate(
      z = ( !!sym(value_col) - med ) / mad_val,
      event_flag = z > z_thresh
    )

  # ----------------------------
  # Keep only candidate points
  # ----------------------------
  df_events <- df %>%
    filter(event_flag) %>%
    arrange(date)

  if (nrow(df_events) == 0) {
    return(list(data = df, events = NULL))
  }

  # ----------------------------
  # Compute gaps between detections
  # ----------------------------
  df_events <- df_events %>%
    mutate(
      gap = as.numeric(date - lag(date)),
      new_event = ifelse(is.na(gap) | gap > max_gap_days, 1, 0),
      event_id = cumsum(new_event)
    )

  # ----------------------------
  # Summarise events
  # ----------------------------
  events <- df_events %>%
    group_by(event_id) %>%
    summarise(
      start_date = min(date),
      end_date   = max(date),
      duration   = n(),
      span_days  = as.numeric(max(date) - min(date)) + 1,
      peak_value = max(!!sym(value_col), na.rm = TRUE),
      mean_value = mean(!!sym(value_col), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration >= min_duration)

  # ----------------------------
  # Tag valid events back
  # ----------------------------
  df_events <- df_events %>%
    left_join(events %>% select(event_id) %>% mutate(valid_event = TRUE),
              by = "event_id") %>%
    mutate(valid_event = ifelse(is.na(valid_event), FALSE, TRUE))

  # Merge back into full dataset
  df <- df %>%
    left_join(df_events %>% select(date, valid_event), by = "date") %>%
    mutate(valid_event = ifelse(is.na(valid_event), FALSE, valid_event))

  return(list(data = df, events = events, median = med, mad = mad_val))
}


# Load your CSV

datapath <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/CALCYOM/lakes_cci_database/extracted_Rw_time_series/"
setwd(datapath)
myfile <- "Rw_L2_KIVU_GLWD00000067_20160101_20251231.csv"
df <- read.csv(myfile)
df$date <- as.Date(df$date)

# Keep one band (e.g., BGR)
df_bgr <- df %>% filter(band == "BGR")

res <- detect_whiting_events(df_bgr, value_col = "p95", z_thresh = 3)

df_out <- res$data
events <- res$events
events


ggplot(df_out)+
  geom_path(aes(date, p90), alpha=0.5, size=.5)+
  geom_point(data = df_out[df_out$valid_event,], aes(date, p95, colour = valid_event), alpha=0.5, size=3)+
  xlab("")+
  # ggtitle(paste0("OLCI L2 - lake ",lakename))+
  theme_bw()


