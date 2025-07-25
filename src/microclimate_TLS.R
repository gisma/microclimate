# =============================================================================
# Title:    TLS LAD Voxel Processing and Profile Preparation
# Author:   creu@gmail.com
# License:  MIT License or similar open license
# Date:     2025
# Purpose:  Normalize and preprocess TLS point cloud data for LAD profile
# =============================================================================

# --- Load required libraries -------------------------------------------------
library(lidR)         # Tools for processing LiDAR data
library(terra)        # Spatial raster and vector data handling
library(sf)           # Simple Features for vector data
library(data.table)   # Fast data manipulation
library(zoo)          # Rolling functions and smoothing
library(here)         # Project-root relative file paths
library(rstudioapi)

# --- Define parameters and file paths ----------------------------------------
zmax <- 40                            # Maximum height threshold for clipping [m]
grain.size <- 1                       # Voxel resolution in all directions [m]
project_root <- here::here()         # Root path of the current project

# Choose LAD method: "linear" or "beer"
# Beer–Lambert Notes:
# - Avoids log(0) and 1 by clipping near-extreme values
# - Use when cumulative light absorption or occlusion is relevant
# - Suitable if extinction coefficient is known or estimated from prior studies
lad_method <- "beer"  # Set to "linear" or "beer"

# Optional: extinction coefficient (used only for Beer–Lambert)
k_extinction <- 0.25


las_file <- file.path(project_root, "data/TLS/tree_08.laz")  
output_voxels <- file.path(project_root, "data/TLS/LAD_voxDF.rds")  
output_array <- file.path(project_root, "data/TLS/lad_array_m2m3.rds")  
output_profile_plot <- file.path(project_root, "data/TLS/lad_vertical_profile.pdf")  
envimet_out <- file.path(project_root,"data/envimet/treeo8_tls_envimet.pld")

# --- Load and normalize TLS point cloud --------------------------------------
las <- lidR::readLAS(las_file)  # Read the LAS/LAZ file (point cloud data)

# Normalize height so ground is at 0
las@data$Z <- las@data$Z - min(las@data$Z, na.rm = TRUE)  

# Limit height range to zmax
maxZ <- min(floor(max(las@data$Z, na.rm = TRUE)), zmax)  
las@data$Z[las@data$Z > maxZ] <- maxZ  # Clip height values above maxZ


### functions
  
  #' Export TLS-based LAD profile as Envi-met 3D Tree (.pld)
  #'
  #' Converts a voxel-based LAD profile into an Envi-met PLANT3D-compatible structure.
  #' LAD values are scaled (if needed), XML structure is wrapped in <ENVI-MET_Datafile>,
  #' and a complete .pld file is generated for direct use in Albero.
  #'
  #' @param lad_df A data.frame with `x`, `y`, `Height_bin`, `LAD_median`
  #' @param ID Unique Plant ID (e.g. "010101")
  #' @param Description Plant description (e.g. "TLS Oak")
  #' @param AlternativeName Latin name (e.g. "Quercus robur")
  #' @param Albedo Albedo (0–1), e.g. 0.18
  #' @param Width, Depth Crown dimensions in m (optional, auto-calculated from LAD spread)
  #' @param RootDiameter Root zone diameter in m
  #' @param cellsize Size of one voxel (default: 1 m)
  #' @param Transmittance Fraction of transmitted radiation (default: 0.3)
  #' @param SeasonProfile Vector of 12 monthly values (0–1), default: oak/beech leaf season
  #' @param BlossomProfile Vector of 12 monthly values (0–1), default: short spring blossom
  #' @param LSystem Logical, include L-System info (default: TRUE)
  #' @param scale_factor Numeric LAD multiplier (default: 3)
  #' @param file_out Output .pld filename (default: "tls_envimet_tree.pld")
  #'
  #' @return Writes an Envi-met compatible .pld file to disk
  #' @export
  export_lad_to_envimet3d <- function(lad_df,
                                      ID = "010101",
                                      Description = "TLS-Derived Tree",
                                      AlternativeName = "Genericus tlsii",
                                      Albedo = 0.18,
                                      Width = NULL,
                                      Depth = NULL,
                                      RootDiameter = 4.5,
                                      cellsize = 1,
                                      Transmittance = 0.3,
                                      SeasonProfile = c(0.3, 0.3, 0.3, 0.4, 0.7, 1, 1, 1, 0.8, 0.6, 0.3, 0.3),
                                      BlossomProfile = c(0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                      LSystem = FALSE,
                                      scale_factor = 3,
                                      file_out = "tls_envimet_tree.pld") {
    require(XML)
    
    # ─── 1. Clean and prepare ─────────────────────────────
    lad_df <- lad_df[!is.na(lad_df$LAD_median), ]
    lad_df$i <- as.integer(factor(lad_df$x))
    lad_df$j <- as.integer(factor(lad_df$y))
    z_map <- setNames(seq_along(sort(unique(lad_df$Height_bin))), sort(unique(lad_df$Height_bin)))
    lad_df$k <- z_map[as.character(lad_df$Height_bin)]
    lad_df$lad_value <- round(lad_df$LAD_median * scale_factor, 5)
    
    dataI <- max(lad_df$i)
    dataJ <- max(lad_df$j)
    zlayers <- max(lad_df$k)
    Width  <- if (is.null(Width)) dataI else Width
    Depth  <- if (is.null(Depth)) dataJ else Depth
    Height <- zlayers * cellsize
    
    # ─── 2. Create <PLANT3D> ──────────────────────────────
    plant_node <- newXMLNode("PLANT3D")
    addChildren(plant_node, newXMLNode("ID", ID))
    addChildren(plant_node, newXMLNode("Description", Description))
    addChildren(plant_node, newXMLNode("AlternativeName", AlternativeName))
    addChildren(plant_node, newXMLNode("Planttype", "0"))
    addChildren(plant_node, newXMLNode("Leaftype", "1"))
    addChildren(plant_node, newXMLNode("Albedo", sprintf("%.5f", Albedo)))
    addChildren(plant_node, newXMLNode("Eps", "0.96000"))
    addChildren(plant_node, newXMLNode("Transmittance", sprintf("%.5f", Transmittance)))
    addChildren(plant_node, newXMLNode("Height", sprintf("%.5f", Height)))
    addChildren(plant_node, newXMLNode("Width", sprintf("%.5f", Width)))
    addChildren(plant_node, newXMLNode("Depth", sprintf("%.5f", Depth)))
    addChildren(plant_node, newXMLNode("RootDiameter", sprintf("%.5f", RootDiameter)))
    addChildren(plant_node, newXMLNode("cellsize", sprintf("%.5f", cellsize)))
    addChildren(plant_node, newXMLNode("xy_cells", dataI))
    addChildren(plant_node, newXMLNode("z_cells", zlayers))
    addChildren(plant_node, newXMLNode("scalefactor", "1.00000"))
    
    # ─── 3. LAD Profile ──────────────────────────────────
    lad_lines <- apply(lad_df[, c("i", "j", "k", "lad_value")], 1, function(r) {
      sprintf("%d,%d,%d,%.5f", r[1], r[2], r[3], r[4])
    })
    lad_node <- newXMLNode("LAD-Profile",
                           attrs = c(type = "sparematrix-3D",
                                     dataI = dataI,
                                     dataJ = dataJ,
                                     zlayers = zlayers,
                                     defaultValue = "0.00000"),
                           .children = paste(lad_lines, collapse = "\n"))
    addChildren(plant_node, lad_node)
    
    # ─── 4. Season + Blossom + L-System ──────────────────
    addChildren(plant_node, newXMLNode("Season-Profile",
                                       paste(sprintf("%.5f", SeasonProfile), collapse = ",")))
    addChildren(plant_node, newXMLNode("Blossom-Profile",
                                       paste(sprintf("%.5f", BlossomProfile), collapse = ",")))
    if (LSystem) {
      addChildren(plant_node, newXMLNode("L-SystemBased", "1"))
      addChildren(plant_node, newXMLNode("Axiom", "F(2)V\\V\\\\V/////B"))
      addChildren(plant_node, newXMLNode("IterationDepth", "3"))
      addChildren(plant_node, newXMLNode("TermLString", "L"))
      addChildren(plant_node, newXMLNode("ApplyTermLString", "1"))
    } else {
      addChildren(plant_node, newXMLNode("L-SystemBased", "0"))
    }
    
    # ─── 5. Wrap in <ENVI-MET_Datafile> ──────────────────
    now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    header_node <- newXMLNode("Header")
    addChildren(header_node, newXMLNode("filetype", "DATA"))
    addChildren(header_node, newXMLNode("version", "1"))
    addChildren(header_node, newXMLNode("revisiondate", now))
    addChildren(header_node, newXMLNode("remark", paste("TLS-based tree generated on", now)))
    addChildren(header_node, newXMLNode("fileInfo", "Generated from R export"))
    addChildren(header_node, newXMLNode("checksum", "32767"))
    addChildren(header_node, newXMLNode("encryptionlevel", "1699612"))
    
    root <- newXMLNode("ENVI-MET_Datafile")
    addChildren(root, header_node)
    addChildren(root, plant_node)
    
    # ─── 6. Write to disk ─────────────────────────────────
    saveXML(root, file = file_out, indent = TRUE, encoding = "UTF-8")
    message("✔ Envi-met PLANT3D (.pld) written to: ", normalizePath(file_out))
  }  

  #' Count LiDAR Returns by Height Slice
#'
#' This function bins LiDAR point heights (`Z`) into integer height slices and counts 
#' the number of returns in each vertical meter. The output is a named list of counts 
#' for each height slice, ensuring all heights from 0 to `maxZ` are represented.
#'
#' @param Z A numeric vector of LiDAR point heights (e.g., Z coordinates from a point cloud).
#' @param maxZ Integer. The maximum vertical extent (in meters) to include in the slicing.
#'
#' @return A named list where each entry corresponds to a 1 m height slice and contains the
#'         number of LiDAR returns within that slice. The names follow the pattern:
#'         `"ground_0_1m"`, `"pulses_1_2m"`, ..., up to `maxZ`.
#'
#' @details
#' - Heights are floored to the nearest lower integer (i.e., `floor(Z)`).
#' - If a height slice from `0` to `maxZ` contains no points, it is still included with value `0`.
#' - Naming convention is consistent with vertical LAD or pulse profile analyses.
#'
#' @examples
#' set.seed(42)
#' Z <- runif(1000, 0, 20)  # Simulated LiDAR point heights
#' result <- pointsByZSlice(Z, maxZ = 25)
#' names(result)
#' result$ground_0_1m
#'
#' @export
pointsByZSlice = function(Z, maxZ) {
  heightSlices = as.integer(Z) # Round down
  zSlice = data.table::data.table(Z=Z, heightSlices=heightSlices)
  sliceCount = stats::aggregate(list(V1=Z), list(heightSlices=heightSlices), length)

  # Ensure full range 0:maxZ is represented
  colRange = 0:maxZ
  addToList = setdiff(colRange, sliceCount$heightSlices)
  n = length(addToList)
  if (n > 0) {
    bindDt = data.frame(heightSlices = addToList, V1 = integer(n))
    sliceCount = rbind(sliceCount, bindDt)
    sliceCount = sliceCount[order(sliceCount$heightSlices),]
  }

  colNames = as.character(sliceCount$heightSlices)
  colNames[1] = "ground_0_1m"
  colNames[-1] = paste0("pulses_", colNames[-1], "_", sliceCount$heightSlices[-1]+1, "m")
  metrics = list()
  metrics[colNames] = sliceCount$V1

  return(metrics)
}
  


#' Preprocess LiDAR Voxel Metrics for LAD or Pulse Analysis
#'
#' Aggregates normalized LiDAR pulse counts in vertical slices using voxel-based spatial binning.
#' Returns both a `terra::SpatRaster` object and a `data.frame` for further analysis.
#'
#' @param normlas A normalized `LAS` or `LAScatalog` object from `lidR`. Must contain height (`Z`) values.
#' @param grain.size Numeric. Horizontal and vertical resolution (in meters) of the voxel grid. Default is `1`.
#' @param maxP Maximum vertical extent (Z) to include (e.g., tree height or canopy top). Default is `zmax`.
#' @param normalize Logical. If `TRUE`, voxel values are divided by voxel volume (m³). Default is `TRUE`.
#' @param as_raster Logical. If `TRUE`, includes a `SpatRaster` object in the output. Default is `TRUE`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{raster}{(optional) A `terra::SpatRaster` object containing vertical pulse metrics per voxel.}
#'   \item{df}{A `data.frame` with columns `X`, `Y`, and one column per vertical slice (e.g. `ground_0_1m`, `pulses_1_2m`, ...).}
#' }
#'
#' @details
#' - Uses the internal function `pointsByZSlice()` to count returns in 1 m bins.
#' - All returns below `0` and above `maxP` are removed before aggregation.
#' - If `normalize = TRUE`, values are interpreted as density (pulses/m³).
#'
#' @examples
#' \dontrun{
#'   las <- readLAS("tree.las")
#'   norm_las <- normalize_height(las)
#'   vox <- preprocess_voxels(norm_las, grain.size = 1, maxP = 30)
#'   head(vox$df)
#' }
#'
#'
#'
#' @importFrom lidR pixel_metrics is.empty filter_poi
#' @importFrom terra xyFromCell values ncell
#' @export
preprocess_voxels <- function(normlas, grain.size = 1, maxP =zmax, normalize = TRUE, as_raster = TRUE) {  
  las <- normlas  
  
  # Filter height range
  las <- filter_poi(las, Z >= 0 & Z <= maxP)  
  if (lidR::is.empty(las)) return(NULL)
  # Determine Z-slices
  maxZ <- floor(max(las@data$Z))  
  maxZ <- min(maxZ, maxP)  
  
  
  # Compute voxel metrics
  func <- formula(paste0("~pointsByZSlice(Z, ", maxZ, ")"))  
  voxels <- pixel_metrics(las, func, res = grain.size)  # Calculate metrics in each voxel (3D grid cell)
  
  # Optionally normalize values by voxel volume
  if (normalize) {
    vvol <- grain.size^3  
    voxels <- voxels / vvol  
  }
  
  # Return as both terra::SpatRaster and data.frame
  result <- list()  
  
  if (as_raster) {
    result$raster <- voxels  
  }
  
  # Convert to data.frame
  xy <- terra::xyFromCell(voxels, seq_len(ncell(voxels)))  
  vals <- terra::values(voxels)  
  df <- cbind(xy, vals)  
  colnames(df)[1:2] <- c("X", "Y")  
  result$df <- df  
  
  return(result)
}


#' Convert TLS voxel pulse data to LAD using Beer–Lambert law with post-scaling
#'
#' @param df A data.frame with pulse columns (from TLS voxelization)
#' @param grainsize Numeric, vertical voxel height (e.g., 1 m)
#' @param k Extinction coefficient (default: 0.3)
#' @param scale_factor Optional multiplicative scale factor (default: 1.2)
#' @param lad_max Optional maximum LAD clamp (e.g. 2.5); set to NULL to disable
#' @param lad_min Optional minimum LAD threshold (e.g. 0.05); set to NULL to disable
#' @param keep_pulses Logical, whether to retain pulse columns (default: FALSE)
#'
#' @return Data.frame with LAD columns added
#' @export
convert_to_LAD_beer <- function(df,
                                grainsize = 1,
                                k = 0.3,
                                scale_factor = 1.2,
                                lad_max = 2.5,
                                lad_min = 0.05,
                                keep_pulses = FALSE) {
  df_lad <- df
  pulse_cols <- grep("^pulses_", names(df_lad), value = TRUE)
  
  for (col in pulse_cols) {
    lad_col <- paste0("lad_", sub("pulses_", "", col))
    p_rel <- df_lad[[col]] / max(df_lad[[col]], na.rm = TRUE)
    
    # Avoid log(0) and 1
    p_rel[p_rel >= 1] <- 0.9999
    p_rel[p_rel <= 0] <- 1e-5
    
    # Apply Beer–Lambert
    lad_vals <- -log(1 - p_rel) / (k * grainsize)
    
    # Apply scaling
    lad_vals <- lad_vals * scale_factor
    
    # Clamp LAD values if needed
    if (!is.null(lad_max)) {
      lad_vals <- pmin(lad_vals, lad_max)
    }
    if (!is.null(lad_min)) {
      lad_vals <- pmax(lad_vals, lad_min)
    }
    
    df_lad[[lad_col]] <- lad_vals
    
    if (!keep_pulses) {
      df_lad[[col]] <- NULL
    }
  }
  
  return(df_lad)
}


#' Convert TLS Pulse Counts to Leaf Area Density (LAD)
#'
#' Transforms vertically binned pulse counts (from voxelized TLS data) into Leaf Area Density (LAD, m²/m³)
#' by normalizing pulse values to a specified LAD maximum.
#'
#' @param df A `data.frame` containing voxelized TLS pulse data. Must include columns starting with `"pulses_"`, 
#'           each representing pulse returns per vertical layer (e.g. `pulses_1_2m`, `pulses_2_3m`, ...).
#' @param grainsize Numeric. The voxel edge length in meters (assumed cubic). Default is `1`.
#' @param LADmax Numeric. The maximum LAD value in m²/m³ for relative scaling. Common values: `4.0`–`6.0`. Default is `5.0`.
#' @param keep_pulses Logical. If `FALSE` (default), the original pulse columns are removed from the output. If `TRUE`, they are retained alongside the LAD columns.
#'
#' @return A modified `data.frame` with new LAD columns (`lad_1_2m`, `lad_2_3m`, ...) in m²/m³, scaled relatively to `LADmax`.
#'
#' @details
#' - Each `pulses_*` column is linearly normalized by the overall maximum value across all vertical bins and locations.
#' - The result is a relative LAD estimate, useful for ecological modeling, input to microclimate simulations (e.g., ENVI-met), or structural analysis.
#' - Voxel volume is implicitly considered constant due to cubic assumption (via `grainsize`) but is not explicitly used here.
#'
#' @examples
#' \dontrun{
#'   df_vox <- readRDS("TLS/voxel_metrics.rds")
#'   lad_df <- convert_to_LAD(df_vox, grainsize = 1, LADmax = 5)
#'   head(names(lad_df))  # Should show lad_* columns
#' }
#'
#' @export
convert_to_LAD <- function(df, grainsize = 1, LADmax = 5.0, keep_pulses = FALSE) {  
  # df: Data frame mit voxelisierten TLS-Daten
# grainsize: Voxelgröße in m (würfelförmig angenommen)
# LADmax: maximaler LAD-Wert (Literaturbasiert, z. B. 5.0 m²/m³)
  df_lad <- df  
  pulse_cols <- grep("^pulses_", names(df_lad), value = TRUE)  
  
  # Schichtanzahl = Anzahl Pulse-Spalten
  n_layers <- length(pulse_cols)  
  
  # Optional: originales Maximum zur linearen Skalierung (relativ)
  max_pulse <- max(df_lad[, pulse_cols], na.rm = TRUE)  
  
  # Umwandlung in LAD (m²/m³) – Skaliert auf LADmax oder absolut (siehe Kommentar)
  for (col in pulse_cols) {
    lad_col <- paste0("lad_", sub("pulses_", "", col))  
    
    # Hier wird RELATIV zu max_pulse skaliert → einfache Normalisierung
    df_lad[[lad_col]] <- (df_lad[[col]] / max_pulse) * LADmax  
    
    # Optional: löschen der Pulse-Spalten
    if (!keep_pulses) {
      df_lad[[col]] <- NULL  
    }
  }
  
  return(df_lad)
}

convert_matrix_to_df <- function(mat) {  
  df <- as.data.frame(mat)  
  colnames(df) <- attr(mat, "dimnames")[[2]]  
  return(df)
}

plot_lad_profiles <- function(lad_df, plotstyle = c("each_median", "all_median", "single_profile"),  
                              single_coords = c(NA, NA)) {
  plotstyle <- match.arg(plotstyle)  
  
  # Combine x and y coordinates into a unique column ID
  lad_df$col_id <- paste(lad_df$x, lad_df$y, sep = "_")  
  x_levels <- sort(unique(lad_df$x))  
  y_levels <- sort(unique(lad_df$y))  
  # Convert x/y coordinates to factor variables for matrix layout
  lad_df$x_f <- factor(lad_df$x, levels = x_levels)  
  lad_df$y_f <- factor(lad_df$y, levels = y_levels)  
  n_x <- length(x_levels)  
  n_y <- length(y_levels)  
  
  # Determine the maximum LAD value for relative normalization
  lad_max <- max(lad_df$LAD_median, na.rm = TRUE)  
  height_range <- range(lad_df$Height_bin, na.rm = TRUE)  
  dx <- 0.8  
  dy <- 0.8  
  
  par(mar = c(5, 5, 4, 5), xpd = TRUE)
  
  
  
  
  # Differentiate by plot type: all profiles, overall profile, or single profile
  if (plotstyle == "each_median") {
    # Load PNG legend
    legend_img <- png::readPNG("doc/output.png")
    
    # Define aspect-preserving image placement
    img_height_units <- 20
    img_width_units <- img_height_units * dim(legend_img)[2] / dim(legend_img)[1]  # preserve ratio
    
    # Define position
    img_x_left <- n_x + 1.5
    img_x_right <- img_x_left + img_width_units
    img_y_bottom <- 0
    img_y_top <- img_y_bottom + img_height_units
    
    # Begin plot
    plot(NA, xlim = c(1, n_x + img_width_units + 4), ylim = c(1, n_y),
         type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Vertical LAD Profiles in XY Matrix", asp = 1.2)
    
    
    # Draw all LAD profiles
    for (i in seq_along(x_levels)) {
      for (j in seq_along(y_levels)) {
        profile <- subset(lad_df, x == x_levels[i] & y == y_levels[j])
        if (nrow(profile) == 0) next
        lad_scaled <- profile$LAD_median / lad_max
        height_scaled <- (profile$Height_bin - min(height_range)) / diff(height_range)
        lines(x = lad_scaled * dx + i,
              y = height_scaled * dy + j,
              col = "darkgreen", lwd = 1)
      }
    }
    
    # Axis labels for ground position
    axis(1, at = 1:n_x, labels = round(x_levels, 1), las = 2)
    axis(2, at = 1:n_y, labels = round(y_levels, 1), las = 2)
    
    # Add the image
    rasterImage(legend_img,
                xleft = img_x_left,
                xright = img_x_right,
                ybottom = img_y_bottom,
                ytop = img_y_top)
    
    
    
  } else if (plotstyle == "all_median") {
    unique_heights <- sort(unique(lad_df$Height_bin))  
    lad_median <- numeric(length(unique_heights))  
    for (i in seq_along(unique_heights)) {
      h <- unique_heights[i]  
      lad_median[i] <- median(lad_df$LAD[lad_df$Height_bin == h], na.rm = TRUE)  
    }
    lad_smooth <- stats::filter(lad_median, rep(1/3, 3), sides = 2)  
    
    plot(
      lad_smooth, unique_heights,
      type = "l",
      col = "darkgreen",
      lwd = 2,
      xlab = "Leaf Area Density (m²/m³)",
      ylab = "Height (m)",
      main = "Vertical LAD Profile (smoothed)",
      xlim = c(0, max(lad_smooth, na.rm = TRUE)),
      ylim = range(unique_heights)
    )
    
    text(
      x = as.numeric(lad_smooth),
      y = unique_heights,
      labels = round(as.numeric(lad_smooth), 1),
      pos = 4,
      cex = 0.7,
      col = "black"
    )
    grid()
    
    
  } else if (plotstyle == "single_profile") {
    x_target <- single_coords[1]  
    y_target <- single_coords[2]  
    tol <- 1e-6  
    
    profile <- subset(lad_df, abs(x - x_target) < tol & abs(y - y_target) < tol)  
    
    if (nrow(profile) == 0) {
      # Show warning if no profile exists for selected coordinates
      warning("No data for the selected coordinates.")
      plot.new()
      title(main = paste("No profile at", x_target, "/", y_target))
      return(invisible(NULL))
    }
    
    # Normalize height and LAD
    height_range <- range(profile$Height_bin, na.rm = TRUE)  
    # Determine the maximum LAD value for relative normalization
    lad_max <- max(profile$LAD_median, na.rm = TRUE)  
    
    height_scaled <- (profile$Height_bin - min(height_range)) / diff(height_range)  
    height_unscaled <- profile$Height_bin
    # Determine the maximum LAD value for relative normalization
    lad_scaled <- profile$LAD_median / lad_max  
    
    plot(
      x = lad_scaled,
      y = height_unscaled, #height_scaled,
      type = "l",
      lwd = 2,
      col = "darkgreen",
      xlab = "LAD (normalized)",
      ylab = "Height (m)",
      main = paste("Profile at", x_target, "/", y_target)
    )
  }
}

# =============================================================================
# Process voxelized LAD data and export as ENVI-met-compatible tree profile
# =============================================================================
message("✔ Process voxelized LAD...")
# --- Preprocess LiDAR data into voxel metrics -------------------------------
vox_out <- preprocess_voxels(las, grain.size = 1, maxP = zmax)  # Calculate vertical pulse metrics
vox_df <- convert_matrix_to_df(vox_out$df)                      # Convert voxel array to data.frame

# --- Convert pulse counts to LAD (Leaf Area Density) ------------------------

# Intelligent method selection
if (lad_method == "beer") {
  message("✔ Using Beer–Lambert LAD conversion...")
  df_lad <- convert_to_LAD_beer(
    vox_df,
    grainsize = 1,
    k = k_extinction,
    scale_factor = 0.4,
    lad_max = 2.5,
    lad_min = 0.0
  )
} else if (lad_method == "linear") {
  message("Using linear LAD conversion...")
  df_lad <- convert_to_LAD(
    vox_df,
    grainsize = 1,
    LADmax = 5.0
  )
} else {
  stop("Unknown LAD conversion method: choose 'linear' or 'beer'")
}


# --- Convert LAD data.frame to SpatRaster -----------------------------------
xy <- df_lad[, c("X", "Y")]                                       # Extract XY coordinates
lad_vals <- df_lad[, grep("^lad_", names(df_lad), value = TRUE)] # Select all LAD layers

lad_raster <- rast(cbind(xy, lad_vals), type = "xyz")            # Create raster from LAD values
plot(lad_raster)                                                 # Visualize raster

# --- Reshape LAD data to long format ----------------------------------------

lad_df <- as.data.frame(lad_raster, xy = TRUE, na.rm = TRUE)     # Convert raster to data.frame

# 1. Extract LAD columns and XY coordinates
pulse_cols <- grep("^lad_", names(lad_df), value = TRUE)
xy_cols <- c("x", "y")  # Adjust to "X", "Y" if needed

# 2. Reshape to long format (one row per LAD layer)
lad_df <- reshape(
  data = lad_df[, c(xy_cols, pulse_cols)],
  varying = pulse_cols,
  v.names = "LAD",
  timevar = "layer",
  times = pulse_cols,
  direction = "long"
)

# 3. Extract z-layer information from column names
lad_df$z_low  <- as.numeric(sub("lad_(\\d+)_.*", "\\1", lad_df$layer))  
lad_df$z_high <- as.numeric(sub("lad_\\d+_(\\d+)m", "\\1", lad_df$layer))  

# 4. Compute mid-point height of each voxel layer
lad_df$Height <- (lad_df$z_low + lad_df$z_high) / 2  

# 5. Round to whole meters to create height classes
lad_df$Height_bin <- round(lad_df$Height)  

# --- Aggregate median LAD per 0.5 × 0.5 m column ----------------------------
setDT(lad_df)  # Use data.table for efficient aggregation

lad_by_column <- lad_df[  
  , .(LAD_median = median(LAD, na.rm = TRUE)), 
  by = .(x, y, Height_bin)
]

# Convert back to regular data.frame
lad_df <- as.data.frame(lad_by_column)

# --- Visualize LAD profiles -------------------------------------------------

# Option 1: Profile in each column
plot_lad_profiles(lad_df, plotstyle = "each_median")

# Option 2: Overall vertical LAD profile (median of all)
plot_lad_profiles(lad_df, plotstyle = "all_median")

# Option 3: Single profile at specified coordinates
plot_lad_profiles(lad_df, plotstyle = "single_profile", single_coords = c(57.5, -94.5))

# --- Export final profile as Envi-met PLANT3D tree --------------------------
export_lad_to_envimet3d(
  lad_df = lad_df,
  Description = "TLS Fagus Sylvatica No.8",
  AlternativeName = "TLS Fagus Sylvatica",
  ID = "TLSFS08",
  Albedo = 0.18,
  Transmittance = 0.25,
  RootDiameter = 6,
  LSystem = FALSE,
  scale_factor = 1,                  # Keep LAD values realistic (no artificial boost)
  file_out = envimet_out  # Output file name
)
rstudioapi::navigateToFile(envimet_out)

