library(terra)
library(sf)
library(XML)

download_soilgrids_wcs <- function(
    lon, lat,
    buffer_m = 500,
    variable = "clay",
    depth = "15-30cm",
    quantile = "mean",
    resolution = 0.00025,
    outfile = "soilgrids_output.tif",
    keep_xml = FALSE
) {
  deg_per_km_lat <- 1 / 111.32
  deg_per_km_lon <- 1 / (40075 * cos(lat * pi / 180) / 360)
  buffer_deg_lat <- buffer_m * deg_per_km_lat / 1000
  buffer_deg_lon <- buffer_m * deg_per_km_lon / 1000
  
  bb <- c(
    xmin = lon - buffer_deg_lon,
    ymin = lat - buffer_deg_lat,
    xmax = lon + buffer_deg_lon,
    ymax = lat + buffer_deg_lat
  )
  
  coverage <- paste0(variable, "_", depth, "_", quantile)
  wcs_url <- paste0("https://maps.isric.org/mapserv?map=/map/", variable, ".map&SERVICE=WCS&VERSION=2.0.1")
  
  xml_file <- paste0("request_", variable, ".xml")
  l1 <- XML::newXMLNode("WCS_GDAL")
  XML::newXMLNode("ServiceURL", wcs_url, parent = l1)
  XML::newXMLNode("CoverageName", coverage, parent = l1)
  XML::saveXML(l1, file = xml_file)
  
  gdal_cmd <- sprintf(
    'gdal_translate "%s" "%s" -tr %.6f %.6f -projwin %.6f %.6f %.6f %.6f -projwin_srs EPSG:4326 -co "TILED=YES" -co "COMPRESS=DEFLATE" -co "PREDICTOR=2" -co "BIGTIFF=YES" -of GTiff -q',
    xml_file, outfile, resolution, resolution,
    bb["xmin"], bb["ymax"], bb["xmax"], bb["ymin"]
  )
  
  message("Running GDAL WCS download for ", variable, "...")
  status <- system(gdal_cmd)
  
  if (!keep_xml) unlink(xml_file)
  if (!file.exists(outfile)) stop("GeoTIFF download failed for: ", variable)
  return(outfile)
}

# Extraktion + Umrechnung
extract_groundp_from_soilgrids <- function(lon, lat, buffer = 500) {
  vars <- c("clay", "sand", "silt", "bdod")
  vals <- list()
  
  for (v in vars) {
    tif <- download_soilgrids_wcs(
      lon = lon, lat = lat,
      variable = v,
      buffer_m = buffer,
      depth = "15-30cm", quantile = "mean",
      outfile = paste0("tmp_", v, ".tif")
    )
    rast <- terra::rast(tif)
    vals[[v]] <- terra::extract(rast, terra::vect(sf::st_sfc(sf::st_point(c(lon, lat)), crs = "EPSG:4326")))[1,2]
  }
  
  soilgrids_to_groundp(
    clay = vals$clay,
    sand = vals$sand,
    silt = vals$silt,
    bdod = vals$bdod
  )
}

# Berechnung groundp
soilgrids_to_groundp <- function(clay, sand, silt, bdod, ocd = NA, ksat = NA, depth_m = 0.3) {
  bd <- bdod / 1000
  theta_s <- 1 - bd / 2.65
  theta_fc <- 0.2576 + 0.0033 * clay + 0.0299 * silt
  theta_wp <- 0.0260 + 0.0050 * clay + 0.0158 * silt
  theta_fc <- pmin(theta_fc, theta_s * 0.95)
  theta_wp <- pmin(theta_wp, theta_fc * 0.9)
  
  if (is.na(ksat)) {
    ksat <- exp(12.012 - 0.0755 * clay + 0.2673 * sand^2 / 100 - 0.058 * bdod / 100)
    ksat <- ksat / 86400  # m/d → m/s
  }
  
  return(list(
    soildepth = depth_m,
    thetas = round(theta_s, 3),
    thetawp = round(theta_wp, 3),
    thetafc = round(theta_fc, 3),
    ksat = round(ksat, 6)
  ))
}
lon <- 8.7500
lat <- 50.8500
groundp <- extract_groundp_from_soilgrids(lon, lat)
print(groundp)