library(terra)
library(sf)
library(XML)

download_soilgrids_wcs <- function(lon, lat, buffer_m = 1000,
                                   variable = "clay",
                                   depth = "15-30cm",
                                   quantile = "mean",
                                   resolution = 250,
                                   outfile = paste0(variable, ".tif"),
                                   keep_xml = FALSE) {
  pt <- st_sfc(st_point(c(lon, lat)), crs = 4326)
  pt_igh <- st_transform(pt, "+proj=igh +datum=WGS84 +units=m +no_defs")
  coords <- st_coordinates(pt_igh)[1,]
  bb <- c(coords[1] - buffer_m, coords[2] - buffer_m,
          coords[1] + buffer_m, coords[2] + buffer_m)
  
  coverage <- paste0(variable, "_", depth, "_", quantile)
  wcs_path <- paste0("https://maps.isric.org/mapserv?map=/map/", variable, ".map")
  wcs_url <- paste(wcs_path, "SERVICE=WCS", "VERSION=2.0.1", sep = "&")
  
  xml_file <- paste0(variable, "_request.xml")
  l1 <- newXMLNode("WCS_GDAL")
  newXMLNode("ServiceURL", wcs_url, parent = l1)
  newXMLNode("CoverageName", coverage, parent = l1)
  saveXML(l1, file = xml_file)
  
  gdal_cmd <- sprintf(
    'gdal_translate "%s" "%s" -tr %d %d -projwin %.3f %.3f %.3f %.3f -projwin_srs \'+proj=igh +datum=WGS84 +units=m +no_defs\' -co "TILED=YES" -co "COMPRESS=DEFLATE" -co "PREDICTOR=2" -co "BIGTIFF=YES" -of GTiff -q',
    xml_file, outfile, resolution, resolution, bb[1], bb[4], bb[3], bb[2]
  )
  system(gdal_cmd)
  if (!keep_xml) unlink(xml_file)
  if (!file.exists(outfile)) stop("Download failed: ", outfile)
  return(outfile)
}

get_groundp_from_soilgrids <- function(lon, lat, buffer_m = 1000) {
  vars <- c("clay", "sand", "silt", "bdod")
  raster_list <- list()
  for (var in vars) {
    tif <- download_soilgrids_wcs(lon, lat, buffer_m, variable = var)
    rast_obj <- terra::rast(tif)
    pt <- terra::vect(matrix(c(lon, lat), ncol = 2), type = "points", crs = "EPSG:4326")
    pt_igh <- terra::project(pt, "+proj=igh +datum=WGS84 +units=m +no_defs")
    val <- terra::extract(rast_obj, pt_igh)[1, 2]
    raster_list[[var]] <- val
  }
  
  clay <- raster_list$clay / 100
  sand <- raster_list$sand / 100
  silt <- raster_list$silt / 100
  bdod <- raster_list$bdod / 1000
  
  groundp <- groundparams
  groundp$clay <- clay
  groundp$sand <- sand
  groundp$silt <- silt
  groundp$bd <- bdod
  return(groundp)
}


library(tidyverse)
library(lubridate)
library(zoo)

prepare_climate_blocks <- function(filepath) {
  rawdata <- read_csv(filepath, na = c("NULL", "", "NA", "NaN", "-9999"), show_col_types = TRUE) %>%
    filter(!is.na(datetime)) %>%
    mutate(datetime = parse_date_time(datetime, orders = c("ymd HMS", "ymd HM"), tz = "UTC")) %>%
    arrange(datetime)
  
  time_seq <- tibble(datetime = seq(from = min(rawdata$datetime),
                                    to = max(rawdata$datetime),
                                    by = "5 min"))
  data_full <- time_seq %>%
    left_join(rawdata, by = "datetime") %>%
    arrange(datetime)
  
  numeric_cols <- names(data_full)[sapply(data_full, is.numeric)]
  data_interp <- data_full %>%
    mutate(across(all_of(numeric_cols), ~ na.approx(., x = datetime, na.rm = FALSE)))
  
  data_interp <- data_interp %>%
    mutate(date = as.Date(datetime),
           hour = hour(datetime),
           minute = minute(datetime)) %>%
    filter(!(is.na(date) | is.na(hour) | is.na(minute)))
  
  daily_blocks <- data_interp %>%
    filter(hour * 60 + minute %% 1440 <= 1435) %>%
    group_by(date) %>%
    filter(n() == 288) %>%
    group_split()
  
  return(daily_blocks)
}


run_micropoint_comparison <- function(clim_block, vegp, soilp, lat, long) {
  climdata_user <- clim_block %>%
    select(obs_time = datetime, temp, relhum, pres, swdown, difrad, lwdown, windspeed, winddir, precip)
  
  paii_user <- PAIgeometry(PAI = vegp$pai, skew = 0.3, spread = 0.7, n = 20)
  
  mout_user <- runpointmodel(climdata = climdata_user, reqhgt = 1,
                             vegp = vegp, paii = paii_user,
                             groundp = soilp, lat = lat, long = long)
  
  mout_default <- runpointmodel(climdata = exampledata$climdata, reqhgt = 1,
                                vegp = forestparams, paii = NA,
                                groundp = groundparams, lat = lat, long = long)
  
  plot(climdata_user$obs_time, mout_user$tair, type = "l", col = "darkgreen",
       ylim = range(c(mout_user$tair, mout_default$tair), na.rm = TRUE),
       ylab = "Air Temperature (°C)", xlab = "Time", main = "Micropoint Modellvergleich")
  lines(exampledata$climdata$obs_time, mout_default$tair, col = "red", lty = 2)
  legend("topright", legend = c("Benutzer", "Standard"), col = c("darkgreen", "red"), lty = c(1, 2))
  
  return(list(user = mout_user, default = mout_default))
}

####
###Hauptskript
# source("scripts/01_soilgrids_download.R")
# source("scripts/02_climate_data_prepare.R")
# source("models/03_run_micropoint_model.R")

library(microclimf)
library(micropoint)

# Parameter
filepath <- "data/energie_bil_wiese.csv"
lon <- 8.77
lat <- 50.82

# Lade Bodenparameter
soilparams_user <- get_groundp_from_soilgrids(lon, lat)
soilparams_user$soildepth <- 0.3

# Lade Klimadaten
daily_blocks <- prepare_climate_blocks(filepath)
first_day <- daily_blocks[[1]]

# Veg-Parameter
vegparams_user <- forestparams
vegparams_user$canht <- 0.6
vegparams_user$pai   <- 1.2

# Starte Modellierung
results <- run_micropoint_comparison(first_day, vegparams_user, soilparams_user, lat, lon)

