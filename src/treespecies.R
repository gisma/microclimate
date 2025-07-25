#------------------------------------------------------------------------------
# Script Type: Processing script
# Script Name: 30_filter_classified_species_map.R
# Author: Chris Reudenbach, creuden@gmail.com
#
# Description:
#   - Loads a classified raster map of tree species
#   - Applies spatial smoothing using OTB ClassificationMapRegularization
#   - Aggregates to 1 m resolution using median filtering
#   - Performs contextual correction: replaces isolated Douglas-fir pixels
#     with Beech or Oak if those are dominant in the local neighborhood
#
# Input:
#   - RDS and GeoTIFF classification of tree species (0.2 m resolution)
#
# Output:
#   - Cleaned and aggregated species raster (1 m resolution)
#
# Dependencies:
#   - OTB 9.1+ with PATH and OTB_APPLICATION_PATH correctly set
#
# Copyright: Chris Reudenbach 2021, GPL (>= 3)
# Git: https://github.com/gisma/envimetR.git
#------------------------------------------------------------------------------

# === Libraries ===
library(terra)           # raster handling
library(RColorBrewer)    # color palettes
library(link2GI)         # OTB integration
library(envimaR)         # environment management
library(tools)           # file name tools
library(mapview)         # interactive maps
library(dplyr)           # data manipulation

# === Environment and paths ===
root_folder <- find_rstudio_root_file()

# Set up OTB environment
otb <- link2GI::linkOTB(searchLocation = "~/apps/OTB-9.1.0-Linux/")
Sys.setenv(OTB_APPLICATION_PATH = file.path(dirname(as.character(otb$pathOTB)), "lib/otb/applications"))
Sys.setenv(PATH = paste(otb$pathOTB, Sys.getenv("PATH"), sep = ":"))

# === Parameters ===
target_res <- 1                # desired resolution in meters
min_tree_height <- 2           # (not used yet)
fn <- "5-25_MOF_rgb"           # image stem
epsg <- 25832                  # UTM32N
sapflow_ext <- raster::extent(477500, 478218, 5631730, 5632500)  # area of interest

# === Class ID legend ===
ts <- data.frame(
  ID = 1:12,
  value = c("agriculture", "alder", "ash", "beech", "douglas_fir", "larch",
            "oak", "pastures", "roads", "settlements", "spruce", "water")
)

#------------------------------------------------------------------------------
# FUNCTION: Replace isolated Douglas-fir with Beech or Oak if dominant around
#------------------------------------------------------------------------------
replace_douglas_in_buche_eiche <- function(rast_input,
                                           window_size = 5,
                                           douglas_value = 5,
                                           target_values = c(4, 7),
                                           target_res=1.0) {
  if (!inherits(rast_input, "SpatRaster")) {
    stop("Input must be a terra::SpatRaster object.")
  }
  
  if (window_size %% 2 == 0) stop("window_size must be odd")
  
  # Focal window matrix (square)
  w <- matrix(1, nrow = window_size, ncol = window_size)
  
  # Calculate local modal value (dominant class in window)
  r_mode <- focal(rast_input, w = w, fun = modal,
                  na.policy = "omit", na.rm = TRUE,
                  progress = "text")
  
  # Identify Douglas-fir pixels and surrounding Beech/Oak dominance
  is_douglas <- rast_input == douglas_value
  is_oak_beech_mode <- r_mode %in% target_values
  replace_mask <- is_douglas & is_oak_beech_mode
  
  # Replace Douglas-fir where Beech or Oak dominate
  r_new <- rast_input
  r_new[replace_mask] <- r_mode[replace_mask]
  

    writeRaster(r_new, sprintf("data/aerial/%s_%sm.tif", "agg_cleand", target_res), overwrite = TRUE)

  
  return(r_new)
}

#------------------------------------------------------------------------------
# STEP 1: Read tree species classification from RDS
#------------------------------------------------------------------------------
sapflow_species <- readRDS("data/aerial/sfprediction_ffs_5-25_MOF_rgb.rds")

# Write to GeoTIFF for further processing
raster::writeRaster(sapflow_species, "data/aerial/prediction_ffs.tif",
                    progress = "text", overwrite = TRUE)

# Crop to sapflow test area
sapflow_species <- raster::crop(sapflow_species, sapflow_ext)
raster::writeRaster(sapflow_species, "data/aerial/prediction_ffs_cut.tif",
                    progress = "text", overwrite = TRUE)

#------------------------------------------------------------------------------
# STEP 2: Run OTB ClassificationMapRegularization (majority filter)
#------------------------------------------------------------------------------
cmr <- parseOTBFunction("ClassificationMapRegularization", otb)
cmr$io.in <- "data/aerial/prediction_ffs.tif"
cmr$io.out <- "data/aerial/majority_out.tif"
cmr$progress <- "true"
cmr$ip.radius <- "1"

filter_treespecies <- runOTB(cmr, gili = otb$pathOTB,
                             quiet = FALSE, retRaster = TRUE)

#------------------------------------------------------------------------------
# STEP 3: Aggregate to 1 m resolution using median
#------------------------------------------------------------------------------
r <- rast("data/aerial/majority_out.tif")
cur_res <- res(r)[1]
fact <- round(target_res / cur_res)

if (target_res <= cur_res) stop("Zielauflösung ist kleiner als aktuelle.")

r_agg <- aggregate(r, fact = fact, fun = median, na.rm = TRUE)

# Build automatic filename
outfile <- sprintf("data/aerial/%s_%sm.tif",
                   tools::file_path_sans_ext(basename("data/aerial/aggregate.tif")),
                   target_res)

# Save aggregated raster
writeRaster(r_agg, outfile, overwrite = TRUE)

#------------------------------------------------------------------------------
# STEP 4: Clean Douglas-fir patches contextually
#------------------------------------------------------------------------------
  species_cleaned <- replace_douglas_in_buche_eiche(
    rast_input = r_agg,
    window_size = 5
  )

#------------------------------------------------------------------------------
# STEP 5: Visualize intermediate steps (interactive)
#------------------------------------------------------------------------------
# m1 <- mapview(sapflow_species, col.regions = brewer.pal(12, "Paired",),at = ts$ID,
#               maxpixels = 13821500)

m2 <- mapview(filter_treespecies, col.regions = brewer.pal(12, "Paired"),at = ts$ID,maxpixels = 13821500)

# m3 <- mapview(r_agg, col.regions = brewer.pal(12, "Paired"),at = ts$ID,
#               maxpixels = 13821500)

m4 <- mapview(species_cleaned, col.regions = brewer.pal(12, "Paired"),at = ts$ID, maxpixels = 13821500)

# Combine all visualizations

m2 + m4
