

#' @title Prepare Full-Day 5-Minute Climate Data Blocks
#'
#' @description
#' This function reads, interpolates, and extracts complete 24-hour blocks of 5-minute resolution
#' climate data from a CSV file. It ensures that only full days with exactly 288 valid time steps
#' (i.e. 24 hours × 12 values/hour) are returned. Missing values in numeric columns are linearly
#' interpolated, assuming equally spaced time steps. 
#'
#' The input CSV file must contain a `datetime` column in one of the formats "YYYY-MM-DD HH:MM[:SS]".
#' Missing values (`NULL`, `NA`, `NaN`, `-9999`) are treated as NA.
#'
#' This is essential for models like `micropoint` that require continuous full-day input.
#'
#' @param filepath Path to a CSV file with a datetime column and climate variables at 5-minute intervals.
#'
#' @return A list of data frames, each representing a complete day with 288 records.
#'
#' @examples
#' \dontrun{
#' blocks <- prepare_climate_blocks("data/energie_bilanz.csv")
#' plot(blocks[[1]]$datetime, blocks[[1]]$Ta_2m, type = "l")
#' }
#'
#' @seealso [micropoint::runpointmodel()]
#' @author Chris Reudenbach
library(readr)
library(dplyr)
library(lubridate)
library(zoo)

prepare_climate_blocks <- function(filepath) {
  # --- Load CSV and parse datetime correctly ---
  rawdata <- read_csv(filepath, na = c("NULL", "", "NA", "NaN", "-9999"), show_col_types = TRUE) 
  climdata_user <- rawdata %>%
    transmute(
      obs_time  = datetime,
      
      # Measured or placeholder values
      temp      = coalesce(Ta_2m, NA_real_),        # °C
      relhum    = coalesce(Huma_2m, NA_real_),      # %
      pres      = 101300,                           # Pa (default if not present)
      
      swdown    = coalesce(rad_bil, NA_real_),      # Shortwave incoming
      difrad    = coalesce(rad_bil * 0.3, NA_real_),# Estimate diffuse as 30% of SW
      lwdown    = coalesce(rad_lw, NA_real_),       # Longwave incoming
      
      windspeed = coalesce(Windspeed_2m, NA_real_), # m/s
      winddir   = 180,                              # Placeholder or default
      precip    = 0                                 # No precipitation assumed
    )
  
  # Parse datetime string to POSIXct
  climdata_user <- climdata_user %>%
    mutate(obs_time = dmy_hm(obs_time, tz = "UTC"),
           date = as.Date(obs_time))
  
  # Split by unique date (returns a list of tibbles, one per day)
  daily_blocks <- split(climdata_user, climdata_user$date)
  
  # Ergänze ggf. fehlenden 00:00-Zeitpunkt mit Kopie der ersten Zeile
  daily_blocks <- lapply(daily_blocks, function(day_df) {
    first_time <- format(day_df$obs_time[1], "%H:%M")
    if (first_time != "00:00") {
      new_row <- day_df[1, ]
      new_row$obs_time <- as.POSIXct(paste0(as.Date(new_row$obs_time), " 00:00"), tz = "UTC")
      day_df <- bind_rows(new_row, day_df) %>% arrange(obs_time)
    }
    return(day_df)
  })
}


#' Extract SoilGrids Soil Profile and Estimate Rooting Depth
#'
#' This function downloads and extracts point-specific soil property values
#' from the ISRIC SoilGrids Web Coverage Service (WCS) for a given geographic location.
#' It supports downloading variables like clay, sand, silt, and bulk density (bdod)
#' across six standard depth intervals. From these, a soil profile is constructed,
#' and the effective soil depth (rooting depth) is estimated using changepoint detection
#' on clay content and bulk density.
#'
#' ## Data Source
#' - ISRIC SoilGrids [https://soilgrids.org](https://soilgrids.org)
#' - Resolution: 250 m (default)
#' - Variables: `"clay"`, `"sand"`, `"silt"` (in g/kg), `"bdod"` (bulk density in kg/m³)
#'
#' ## Units and Scaling
#' - `clay`, `sand`, `silt`: Scaled from g/kg to volume fraction → divided by 1000
#' - `bdod`: Converted from kg/m³ to g/cm³ → divided by 1000
#'
#' ## Soil Depth Estimation
#' - Changepoint detection (`changepoint::cpt.mean()`) is applied to both clay and bdod profiles
#' - The shallower changepoint is interpreted as the limiting soil depth
#' - This approach avoids relying on hard thresholds and adapts to local profile structure
#'
#' ## Weighted Means
#' - Only profile layers above the estimated depth are used
#' - Weighting is based on thickness of each horizon
#'
#' @param lon Numeric. Longitude in decimal degrees (WGS84)
#' @param lat Numeric. Latitude in decimal degrees (WGS84)
#' @param variables Character vector of SoilGrids variables to retrieve (default: clay, sand, silt, bdod)
#' @param resolution Numeric. Target resolution in meters (default: 250)
#' @param outdir Character. Local directory for caching WCS GeoTIFFs (default: "data/tmp_soilgrids")
#'
#' @return A list with:
#' - `profile_df`: full profile dataframe (clay, sand, silt, bdod, thickness, z_mid, include)
#' - `weighted_means`: average values weighted by valid depth
#' - `estimated_soildepth`: inferred effective rooting depth
#' - `groundparams`: ready-to-use `groundparams` list for `micropoint`
#' @examples
#' # Example: Grubenwiese climate flux station near Caldern
#' # Coordinates: lat = 50.840461, lon = 8.683241
#' 
#' library(microclimf)
#' library(micropoint)
#' 
#' # Extract soil profile and estimate soil depth
#' soil_result <- get_soil_profile_from_soilgrids(
#'   lon = 8.683241,
#'   lat = 50.840461
#' )
#' 
#' # Inspect profile data
#' print(soil_result$profile_df)
#' 
#' # Show estimated rooting depth
#' print(soil_result$estimated_soildepth)
#' 
#' # Use in micropoint model
#' soilparams_user <- soil_result$groundparams
#' 
#' # Modify default forest parameters
#' vegparams_user <- forestparams
#' vegparams_user$canht <- 0.6
#' vegparams_user$pai <- 1.2
#' 
#' # Load daily weather block
#' clim_blocks <- prepare_climate_blocks("data/energie_bil_wiese.csv")
#' first_day <- clim_blocks[[1]]
#' 
#' # Run model comparison (user soil vs default)
#' results <- run_micropoint_comparison(
#'   weather = first_day,
#'   vegp = vegparams_user,
#'   groundp = soilparams_user,
#'   lat = 50.840461,
#'   long = 8.683241
#' )

#' @export
get_soil_profile_from_soilgrids <- function(lon, lat,
                                            variables = c("clay", "sand", "silt", "bdod"),
                                            resolution = 250,
                                            outdir = here::here("data/tmp_soilgrids")) {
  library(terra)
  library(XML)
  library(here)
  library(changepoint)
  
  make_safe_name <- function(var, depth, lon, lat) {
    lon_str <- gsub("\\.", "_", sprintf("%.4f", lon))
    lat_str <- gsub("\\.", "_", sprintf("%.4f", lat))
    gsub("[^a-zA-Z0-9_\\-]", "", sprintf("%s_%s_lon%s_lat%s", var, depth, lon_str, lat_str))
  }
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")
  profile <- data.frame(depth = depths)
  
  for (var in variables) {
    message(sprintf("→ Processing variable: %s", var))
    vals <- numeric(length(depths))
    
    for (i in seq_along(depths)) {
      coverage <- paste0(var, "_", depths[i], "_mean")
      base <- make_safe_name(var, depths[i], lon, lat)
      xml_file <- file.path(outdir, paste0(base, ".xml"))
      tif_file <- file.path(outdir, paste0(base, ".tif"))
      
      message(sprintf("   ↪ Depth: %s", depths[i]))
      
      if (!file.exists(xml_file)) {
        message(sprintf("     → Writing WCS XML: %s", basename(xml_file)))
        node <- newXMLNode("WCS_GDAL")
        newXMLNode("ServiceURL", paste0("https://maps.isric.org/mapserv?map=/map/", var, ".map&SERVICE=WCS&VERSION=2.0.1"), parent = node)
        newXMLNode("CoverageName", coverage, parent = node)
        saveXML(node, file = xml_file)
      } else {
        message("     → XML already exists")
      }
      
      abs_xml <- normalizePath(xml_file, winslash = "/", mustWork = FALSE)
      abs_tif <- normalizePath(tif_file, winslash = "/", mustWork = FALSE)
      
      if (!file.exists(abs_tif)) {
        message(sprintf("     → Downloading GeoTIFF to: %s", basename(abs_tif)))
        gdal_cmd <- paste(
          "gdalwarp",
          "-overwrite",
          "-te",
          sprintf("%.6f", lon - 0.001),
          sprintf("%.6f", lat - 0.001),
          sprintf("%.6f", lon + 0.001),
          sprintf("%.6f", lat + 0.001),
          "-ts", 1, 1,
          "-t_srs", "EPSG:4326",
          "-of", "GTiff",
          shQuote(abs_xml),
          shQuote(abs_tif)
        )
        system(gdal_cmd)
      } else {
        message("     → GeoTIFF already exists, skipping download")
      }
      
      if (!file.exists(abs_tif)) {
        warning(sprintf("⚠️ Missing data for %s at %s", var, depths[i]))
        vals[i] <- NA
        next
      }
      
      rast <- terra::rast(abs_tif)
      pt <- terra::vect(matrix(c(lon, lat), ncol = 2), type = "points", crs = "EPSG:4326")
      val <- terra::extract(rast, pt)[1, 2]
      vals[i] <- val
      message(sprintf("     → Extracted value: %s = %.4f", var, val))
    }
    
    profile[[var]] <- vals
  }
  
  # Unit conversion (raw SoilGrids in g/kg and kg/m³)
  profile$clay <- profile$clay / 1000
  profile$sand <- profile$sand / 1000
  profile$silt <- profile$silt / 1000
  profile$bdod <- profile$bdod / 1000
  profile$z_mid <- c(0.025, 0.10, 0.225, 0.45, 0.80, 1.5)
  
  # Rooting depth estimation based on changepoint in clay and bdod
  estimate_depth_changepoint <- function(vec, z_mid, method = "PELT", pen.value = "SIC") {
    if (all(is.na(vec))) return(NA)
    vec_interp <- approx(seq_along(vec), vec, xout = seq_along(vec), rule = 2)$y
    cpt <- cpt.mean(vec_interp, method = method, penalty = pen.value)
    cp_index <- cpts(cpt)
    return(if (length(cp_index) == 0) max(z_mid) else z_mid[min(cp_index)])
  }
  
  depth_bdod <- estimate_depth_changepoint(profile$bdod, profile$z_mid)
  depth_clay <- estimate_depth_changepoint(profile$clay, profile$z_mid)
  estimated_soildepth <- min(depth_bdod, depth_clay, na.rm = TRUE)
  
  # Restrict profile and compute weighted means
  thickness <- c(0.05, 0.10, 0.15, 0.30, 0.40, 1.00)
  profile$thickness <- thickness
  profile$include <- profile$z_mid <= estimated_soildepth
  
  weighted_means <- sapply(c("clay", "sand", "silt", "bdod"), function(var) {
    use <- profile$include & !is.na(profile[[var]])
    if (all(!use)) return(NA)
    sum(profile[[var]][use] * profile$thickness[use]) / sum(profile$thickness[use])
  })
  
  # Build groundparams object
  groundp <- groundparams
  groundp$clay <- weighted_means["clay"]
  groundp$sand <- weighted_means["sand"]
  groundp$silt <- weighted_means["silt"]
  groundp$bd   <- weighted_means["bdod"]
  groundp$soildepth <- estimated_soildepth
  
  message(sprintf("✅ Estimated soildepth: %.2f m", estimated_soildepth))
  
  return(list(
    profile_df = profile,
    weighted_means = weighted_means,
    estimated_soildepth = estimated_soildepth,
    groundparams = groundp
  ))
}


#### ============================================================================
#### Hauptskript: Mikropoint-Modellierung mit benutzerdefiniertem Bodenprofil
#### ============================================================================

# --- Bibliotheken laden ---
library(tidyverse)
library(lubridate)
library(zoo)
library(terra)
library(XML)
library(here)
library(changepoint)
library(micropoint)
library(knitr)

# ------------------------------------------------------------------------------
# 1. Parameter setzen
# ------------------------------------------------------------------------------

filepath <- "data/energie_bil_wiese.csv"
lon <- 8.683241
lat <- 50.840461
reqhgt <- 15

# ------------------------------------------------------------------------------
# 2. Bodendaten von SoilGrids laden und verarbeiten
# ------------------------------------------------------------------------------

soil <- get_soil_profile_from_soilgrids(lon = lon, lat = lat)
print(soil$profile_df)
print(soil$estimated_soildepth)
print(soil$groundparams)
soilp <- soil$groundparams  # als groundp für Modell verwenden

# ------------------------------------------------------------------------------
# 3. Klimadaten einlesen und Tagesblöcke vorbereiten
# ------------------------------------------------------------------------------

daily_blocks <- prepare_climate_blocks(filepath)

# ------------------------------------------------------------------------------
# 4. Schleife über zwei Tage (day 1 und 2)
# ------------------------------------------------------------------------------

for (day_index in 1:2) {
  
  first_day <- daily_blocks[[1]]
  print(first_day)
  
  # ----------------------------------------------------------------------------
  # 5. Vegetationsparameter definieren
  # ----------------------------------------------------------------------------
  
  vegparams_user <- list(
    h     = 2,   # Measurement height above ground (m), e.g. where sensors are placed
    pai   = 1,    # Total plant area index (m²/m²), typically 1–5 (use 4 for forest)
    x     = 1.0,    # Leaf angle distribution (1 = spherical)
    clump = 0.1,    # Clumping index (0 = fully clumped, 1 = random)
    lref  = 0.4,    # Leaf reflectance (0.2–0.5 typical)
    ltra  = 0.2,    # Leaf transmittance (0.05–0.3 typical)
    leafd = 0.05,   # Leaf width or characteristic leaf dimension (m)
    em    = 0.97,   # Surface emissivity (0.95–0.99 for vegetation)
    gsmax = 0.33,   # Maximum stomatal conductance (mol/m²/s), typical 0.2–0.6
    q50   = 100,    # PPFD at which gs is 50% of gsmax (µmol/m²/s), typical 50–200
    hgt   = 15.0    # Vegetation height (m), e.g. 10–30 m for forest
  )
  
  # Assign S3 class
  class(vegparams_user) <- "vegparams"
  paii_user <- paii_user <- PAIgeometry(
    PAI   = vegparams_user$pai,
    skew  = 4.9,     # niedrig = Konzentration in oberen Schichten
    spread = 1.9,    # höherer spread → verteilt sich weiter nach oben
    n     = 5       # Schichtenanzahl
  )
  
  
  # ----------------------------------------------------------------------------
  # 6. Mikropoint-Modellläufe (benutzerdefiniert vs. Standard)
  # ----------------------------------------------------------------------------
  
  mout_user <- runpointmodel(
    climdata = first_day,
    reqhgt   = reqhgt,
    vegp     = vegparams_user,
    paii     = paii_user,
    groundp  = soilp,
    lat      = lat,
    long     = lon
  )
  
  mout_default <- runpointmodel(
    climdata = first_day,
    reqhgt   = reqhgt,
    vegp     = vegparams_user,
    paii     = paii_user,
    groundp  = groundparams,
    lat      = lat,
    long     = lon
  )
  
  # ----------------------------------------------------------------------------
  # 7. Visualisierung: Zeitverlauf und Differenzplot
  # ----------------------------------------------------------------------------
  
  par(mfrow = c(1, 2))  # zwei Plots nebeneinander
  day_label <- format(first_day$obs_time[1], "%Y-%m-%d")
  
  # Plot 1: Zeitreihe Lufttemperatur
  plot(mout_user$obs_time, mout_user$tair, type = "l",
       xlab = "Time", ylab = "Air Temperature (°C)",
       ylim = range(c(mout_user$tair, mout_default$tair), na.rm = TRUE),
       col = rgb(1, 0, 0, 0.6), lwd = 2, 
       main = paste("User vs Default –", day_label))
  lines(mout_default$obs_time, mout_default$tair,
        col = rgb(0, 0, 1, 0.6), lty = 2, lwd = 2)
  legend("topright", legend = c("User", "Default"), col = c("red", "blue"),
         lty = c(1, 2), bty = "n")
  
  # Plot 2: Differenz User - Default
  temp_diff <- mout_user$tair - mout_default$tair
  plot(mout_user$obs_time, temp_diff, type = "l",
       xlab = "Time", ylab = "Δ Air Temperature (°C)",
       main = paste("User vs Default –", day_label),
       col = "darkorange", lwd = 2)
  abline(h = 0, lty = 2, col = "gray")
  
  par(mfrow = c(1, 1))  # Layout zurücksetzen
  
  # ----------------------------------------------------------------------------
  # 8. Tabellenvergleich der Bodenparameter: Benutzer vs. Standard (Loam)
  # ----------------------------------------------------------------------------
  
  default_soil <- soilparams %>% filter(Soil.type == "Loam") %>% slice(1)
  
  custom_soil <- tibble(
    Parameter = c("Smax", "Smin", "alpha", "n", "Ksat", "Vq", "Vm", "Mc", "rho", "b", "Psie"),
    Custom = c(
      soilp$Smax, soilp$Smin, soilp$alpha, soilp$n, soilp$Ksat,
      soilp$Vq, soilp$Vm, soilp$Mc, soilp$rho, soilp$b, soilp$Psie
    ),
    Default = c(
      default_soil$Smax, default_soil$Smin, default_soil$alpha, default_soil$n, default_soil$Ksat,
      default_soil$Vq, default_soil$Vm, default_soil$Mc, default_soil$rho, default_soil$b, default_soil$psi_e
    )
  )
  
  print(kable(custom_soil, caption = paste("Comparison of Soil Parameters (Day", day_index, ")")))
  
}


# run model for height mid-way to canopy top (6 cm above ground)
mout1 <- runpointmodel(first_day, reqhgt = .0, vegparams_user, paii_user,  groundparams, lat = lat, long = lon)
mout2 <- runpointmodel(first_day, reqhgt = -0.1, vegparams_user, paii_user,  groundparams, lat = lat, long = lon)
mout3 <- runpointmodel(first_day, reqhgt = -0.5, vegparams_user, paii_user,  groundparams, lat = lat, long = lon)
tme <- as.POSIXct(first_day$obs_time)
tme <- as.POSIXct(first_day$obs_time)

plot(mout1$tground ~ tme, type = "l", ylim = c(-5, 40),
     col = rgb(1, 0, 0, 0.5), xlab = "Time", ylab = "Temperature (°C)")

par(new = TRUE)
plot(mout2$tground ~ tme, type = "l", ylim = c(-5, 40),
     col = rgb(0, 0, 0, 0.5), lwd = 2, xlab = "", ylab = "")

par(new = TRUE)
plot(mout3$tground ~ tme, type = "l", ylim = c(-5, 40),
     col = rgb(0, 0, 1, 0.5), lwd = 3, xlab = "", ylab = "")

legend("topright",
       legend = c("reqhgt = 0 m", "reqhgt = -0.05 m", "reqhgt = -0.5 m"),
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
       lwd = c(1, 2, 3),
       bty = "n")


# run model for height mid-way to canopy top (6 cm above ground)
# run model for height mid-way to canopy top (6 cm above ground)
mout1 <- runpointmodel(first_day, reqhgt = 2, vegparams_user, paii = paii_user,  groundparams, lat = lat, long = lon)
mout2 <- runpointmodel(first_day, reqhgt = 7, vegparams_user, paii = paii_user,  groundparams, lat = lat, long = lon)
mout3 <- runpointmodel(first_day, reqhgt = 18., vegparams_user, paii = paii_user,  groundparams, lat = lat, long = lon)
tme <- as.POSIXct(first_day$obs_time)
# Compute global y-axis range from all model outputs
y_range <- range(
  mout1$tground,
  mout2$tair,
  mout3$tair,
  na.rm = TRUE
)
# Define time axis
tme <- as.POSIXct(first_day$obs_time)

# Compute dynamic y-axis range across all outputs
y_range <- range(
  mout1$tground,
  mout2$tair,
  mout3$tair,
  na.rm = TRUE
)

# Liste der Modelle
models <- list(
  "reqhgt = 0 m"  = mout1,
  "reqhgt = 5 m"  = mout2,
  "reqhgt = 10 m" = mout3
)

# Zeitachse
tme <- as.POSIXct(first_day$obs_time)

# Set up 2 rows × 3 columns layout
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))  # engeres margin

# Schleife über Modelle (jeweils Vergleich mit first_day$temp)
for (i in seq_along(models)) {
  model <- models[[i]]
  label <- names(models)[i]
  
  # Gemeinsame y-Achse für fairen Vergleich
  y_range <- range(first_day$temp, model$tair, na.rm = TRUE)
  
  # Plot gemessene Temperatur
  plot(tme, first_day$temp, type = "l", col = "black", lwd = 1.5,
       ylim = y_range,
       xlab = "Time", ylab = "Temperature (°C)",
       main = paste("Measured vs.", label))
  
  # Modellierte Temperatur überlagern
  lines(tme, model$tair, col = "red", lwd = 2)
  
  legend("topleft", legend = c("Measured", "Modelled"), col = c("black", "red"),
         lty = 1, lwd = 2, bty = "n", cex = 0.9)
}

# Falls du weitere Modelle (mout4 bis mout6) hast, hier erweitern:
# models[["reqhgt = 15 m"]] <- mout4
# ...





reqhgt = c(0.12, 1, 5, 10, 20, 28)  # Messhöhen von Boden bis Baumkrone
vegparams_user <- list(
  hgt   = 28,   # canopy height
  pai   = 4,    # plant area index
  x     = 1.0,
  clump = 0.15,
  lref  = 0.3,
  ltra  = 0.15,
  leafd = 0.05,
  em    = 0.97,
  gsmax = 0.3,
  q50   = 100
)
class(vegparams_user) <- "vegparams"

paii_user <- PAIgeometry(PAI = vegparams_user$pai, skew = 0.3, spread = 0.7, n = 30)

# ------------------------------------------------------------------------------
# 3. Define soil parameters
# ------------------------------------------------------------------------------
soilp <- default_groundparams()

# ------------------------------------------------------------------------------
# 4. Run model for multiple heights
# ------------------------------------------------------------------------------
reqhgt_seq <- c( 28,150,1500)
models <- lapply(reqhgt_seq, function(h) {
  runpointmodel(
    climdata = climdata,
    reqhgt = h,
    vegp = vegparams_user,
    paii = paii_user,
    groundp = soilp,
    lat = 49.96807,
    long = -5.215668
  )
})
names(models) <- paste0("h=", reqhgt_seq, "m")

# ------------------------------------------------------------------------------
# 5. Plot air temperature profiles for each height
# ------------------------------------------------------------------------------
colors <- rainbow(length(reqhgt_seq))
plot(models[[1]]$obs_time, models[[1]]$tair, type = "l", ylim = c(5, 35),
     xlab = "Time", ylab = "Air Temperature (°C)", col = colors[1], lwd = 2,
     main = "Tair at Different Heights in Canopy")
for (i in 2:length(models)) {
  lines(models[[i]]$obs_time, models[[i]]$tair, col = colors[i], lwd = 2)
}
legend("topright", legend = names(models), col = colors, lwd = 2, bty = "n")




