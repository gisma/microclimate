library(terra)
library(sf)
library(XML)
library(microclimf)
library(micropoint)

# Funktion zum Herunterladen eines SoilGrids-Layers via WCS
download_soilgrids_wcs <- function(
    lon, lat, buffer_m = 1000,
    variable = "clay", 
    depth = "15-30cm", 
    quantile = "mean",  # options: Q0.05, Q0.5, Q0.95, mean
    resolution = 250,
    outfile = paste0(variable, ".tif"),
    keep_xml = FALSE
) {
  # Reprojektionsvorbereitung
  pt <- st_sfc(st_point(c(lon, lat)), crs = 4326)
  pt_igh <- st_transform(pt, "+proj=igh +datum=WGS84 +units=m +no_defs")
  coords <- st_coordinates(pt_igh)[1,]
  
  bb <- c(
    coords[1] - buffer_m,
    coords[2] - buffer_m,
    coords[1] + buffer_m,
    coords[2] + buffer_m
  )
  
  coverage <- paste0(variable, "_", depth, "_", quantile)
  wcs_path <- paste0("https://maps.isric.org/mapserv?map=/map/", variable, ".map")
  wcs_url <- paste(wcs_path, "SERVICE=WCS", "VERSION=2.0.1", sep = "&")
  
  # XML schreiben
  xml_file <- paste0(variable, "_request.xml")
  l1 <- newXMLNode("WCS_GDAL")
  newXMLNode("ServiceURL", wcs_url, parent = l1)
  newXMLNode("CoverageName", coverage, parent = l1)
  saveXML(l1, file = xml_file)
  
  # GDAL-Aufruf
  gdal_cmd <- sprintf(
    'gdal_translate "%s" "%s" -tr %d %d -projwin %.3f %.3f %.3f %.3f -projwin_srs \'+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs\' -co "TILED=YES" -co "COMPRESS=DEFLATE" -co "PREDICTOR=2" -co "BIGTIFF=YES" -of GTiff -q',
    xml_file, outfile, resolution, resolution, bb[1], bb[4], bb[3], bb[2]
  )
  
  message(sprintf("→ Downloading: %s", coverage))
  system(gdal_cmd)
  if (!keep_xml) unlink(xml_file)
  if (!file.exists(outfile)) stop("Download failed: ", outfile)
  return(outfile)
}

# Funktion: Alle Parameter laden und extrahieren
get_groundp_from_soilgrids <- function(lon, lat, buffer_m = 1000) {
  vars <- c("clay", "sand", "silt", "bdod")
  raster_list <- list()
  
  for (var in vars) {
    tif <- download_soilgrids_wcs(lon, lat, buffer_m, variable = var)
    rast_obj <- terra::rast(tif)
    pt <- terra::vect(matrix(c(lon, lat), ncol = 2), type = "points", crs = "EPSG:4326")
    pt_igh <- terra::project(pt, "+proj=igh +datum=WGS84 +units=m +no_defs")
    val <- terra::extract(rast_obj, pt_igh)[1,2]
    raster_list[[var]] <- val
  }
  
  # Werte zur Sicherheit normalisieren
  clay <- raster_list$clay / 100
  sand <- raster_list$sand / 100
  silt <- raster_list$silt / 100
  bdod <- raster_list$bdod / 1000  # kg/m³ → g/cm³
  
  # groundp bauen
  groundp <- groundparams
  groundp$clay <- clay
  groundp$sand <- sand
  groundp$silt <- silt
  groundp$bd <- bdod
  
  return(groundp)
}

lon <- 8.682639
lat <- 50.837055

groundp_custom <- get_groundp_from_soilgrids(lon, lat)
print(groundp_custom)

