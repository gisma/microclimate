#' --- ENVI-met 3DPLANT column generator from voxelized ALS data ---
#'
#' @author Chris Reudenbach
#'
#' This script performs the full pipeline from voxelization of ALS data to the export
#' of ENVI-met compatible 3DPLANT profiles. It includes LAD computation, clustering
#' of similar profiles, and export of:
#' - Point geometries with shared ENVIMET_IDs
#' - XML-based 3D plant profile database (.pld)
#'
#' Main steps:
#' 1. Voxelization of normalized ALS point cloud
#' 2. LAD computation per voxel
#' 3. Reshaping and clustering of LAD profiles
#' 4. Export of clustered ENVIMET_IDs to GIS and XML format

# Load required libraries
library(lidR)
library(terra)
library(dplyr)
library(sf)
library(here)
library(XML)
library(stats)

# Parameters
las_file <- here("data/ALS/las_mof.las")  # Path to ALS point cloud file
res_xy <- 2                                # Horizontal resolution of voxels (meters)
res_z  <- 2                                # Vertical resolution of voxels (meters)
k      <- 0.3                               # Light extinction coefficient for LAD
scale_factor <- 1.2                        # Optional scaling factor for LAD values
crs_code <- 25832                          # EPSG code of target CRS
output_gpkg <- "output/envimet_tree_points.gpkg"      # Output vector layer with tree positions
xml_output_file <- "output/tls_envimet_trees.pld"      # Output PLANT3D XML file
species_name <- "Fagus_sylvatica"         # Tree species name used in export
n_clusters <- 100                          # Number of LAD profile clusters

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Read and normalize LAS
las <- readLAS(las_file)
crs(las) <- "EPSG:25832"
las <- normalize_height(las, knnidw(k = 6, p = 2))
las <- filter_poi(las, Z > 0)  # Remove ground points

# Voxelize point cloud (counts per voxel)
voxels <- voxel_metrics(las, ~length(Z), res = res_xy, dz = res_z)

#' Convert voxel metrics to LAD profiles per column
#'
#' @param df voxel_metrics result
#' @param res_z vertical resolution (meters)
#' @param k extinction coefficient
#' @param scale_factor LAD scaling multiplier
#' @return data.frame with LAD values per (x,y,z)
convert_voxel_lad_long <- function(df, res_z = 2, k = 0.3, scale_factor = 1.2) {
  if (!all(c("X", "Y", "Z", "V1") %in% names(df))) stop("Missing required columns (X, Y, Z, V1)")
  df$xy_id <- paste(df$X, df$Y, sep = "_")
  df_split <- split(df, df$xy_id)
  lad_df_list <- lapply(df_split, function(group) {
    max_p <- max(group$V1, na.rm = TRUE)
    rel_p <- pmin(pmax(group$V1 / max_p, 1e-5), 0.9999)
    lad <- -log(1 - rel_p) / (k * res_z)
    lad <- lad * scale_factor
    data.frame(x = group$X, y = group$Y, z = group$Z, lad = lad)
  })
  do.call(rbind, lad_df_list)
}

# Compute LAD values
lad_df <- convert_voxel_lad_long(voxels, res_z = res_z, k = k, scale_factor = scale_factor)

# Prepare clustering
lad_df$xy_key <- paste(lad_df$x, lad_df$y)
lad_matrix <- lad_df %>% 
  tidyr::pivot_wider(names_from = z, values_from = lad, values_fill = 0) %>%
  column_to_rownames("xy_key") %>%
  as.matrix()

# Cluster LAD profiles using k-means
clustering <- kmeans(lad_matrix, centers = n_clusters, nstart = 10)
lad_df$cluster <- clustering$cluster[match(lad_df$xy_key, rownames(lad_matrix))]

#' Convert integer cluster to 6-character alphanumeric ENVIMET ID
#'
#' @param n integer to convert
#' @param width number of characters after the S-prefix
#' @return character string of format "Sxxxxx"
int_to_base36 <- function(n, width = 5) {
  chars <- c(0:9, LETTERS)
  base <- length(chars)
  result <- character()
  while (n > 0) {
    result <- c(chars[(n %% base) + 1], result)
    n <- n %/% base
  }
  result <- paste(result, collapse = "")
  padded <- sprintf(paste0("%0", width, "s"), result)
  paste0("S", substr(gsub(" ", "0", padded), 1, width))
}

# Assign ENVIMET_IDs per cluster
cluster_ids <- unique(lad_df$cluster)
cluster_mapping <- data.frame(
  cluster = cluster_ids,
  ENVIMET_ID = sapply(cluster_ids, int_to_base36)
)
lad_df <- left_join(lad_df, cluster_mapping, by = "cluster")

# Export point layer with one record per LAD profile
point_df <- lad_df[!duplicated(lad_df$xy_key), c("x", "y", "ENVIMET_ID")]
sf_points <- st_as_sf(point_df, coords = c("x", "y"), crs = crs_code)
st_write(sf_points, output_gpkg, delete_layer = TRUE)

#' Export clustered LAD profiles to ENVI-met PLANT3D XML
#'
#' @param lad_df data.frame with z, lad, and ENVIMET_ID
#' @param file_out output file path for .pld file
export_lad_to_envimet3d <- function(lad_df, file_out = "tls_envimet_tree.pld") {
  lad_df <- lad_df[!is.na(lad_df$lad), ]
  lad_df$i <- as.integer(factor(lad_df$x))
  lad_df$j <- as.integer(factor(lad_df$y))
  z_map <- setNames(seq_along(sort(unique(lad_df$z))), sort(unique(lad_df$z)))
  lad_df$k <- z_map[as.character(lad_df$z)]
  lad_df$lad_value <- round(lad_df$lad * scale_factor, 5)
  
  tree_ids <- unique(lad_df$ENVIMET_ID)
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  root <- newXMLNode("ENVI-MET_Datafile")
  header_node <- newXMLNode("Header")
  addChildren(header_node, newXMLNode("filetype", "DATA"))
  addChildren(header_node, newXMLNode("version", "1"))
  addChildren(header_node, newXMLNode("revisiondate", now))
  addChildren(header_node, newXMLNode("remark", "Clustered TLS-based trees"))
  addChildren(header_node, newXMLNode("fileInfo", "Clustered LAD Trees"))
  addChildren(header_node, newXMLNode("checksum", "32767"))
  addChildren(header_node, newXMLNode("encryptionlevel", "1699612"))
  addChildren(root, header_node)
  
  for (id in tree_ids) {
    tree_df <- lad_df[lad_df$ENVIMET_ID == id, ]
    profile <- tree_df %>% group_by(z = k) %>% summarise(lad_value = mean(lad_value, na.rm = TRUE))
    zlayers <- max(profile$z)
    dataI <- 1
    dataJ <- 1
    Height <- zlayers * res_z
    
    plant_node <- newXMLNode("PLANT3D")
    addChildren(plant_node, newXMLNode("ID", id))
    addChildren(plant_node, newXMLNode("Description", "Clustered TLS Tree"))
    addChildren(plant_node, newXMLNode("AlternativeName", species_name))
    addChildren(plant_node, newXMLNode("Planttype", "0"))
    addChildren(plant_node, newXMLNode("Leaftype", "1"))
    addChildren(plant_node, newXMLNode("Albedo", sprintf("%.5f", 0.18)))
    addChildren(plant_node, newXMLNode("Eps", "0.96000"))
    addChildren(plant_node, newXMLNode("Transmittance", sprintf("%.5f", 0.3)))
    addChildren(plant_node, newXMLNode("Height", sprintf("%.5f", Height)))
    addChildren(plant_node, newXMLNode("Width", sprintf("%.5f", 1)))
    addChildren(plant_node, newXMLNode("Depth", sprintf("%.5f", 1)))
    addChildren(plant_node, newXMLNode("RootDiameter", sprintf("%.5f", 4.5)))
    addChildren(plant_node, newXMLNode("cellsize", sprintf("%.5f", res_z)))
    addChildren(plant_node, newXMLNode("xy_cells", dataI))
    addChildren(plant_node, newXMLNode("z_cells", zlayers))
    addChildren(plant_node, newXMLNode("scalefactor", "1.00000"))
    
    lad_lines <- apply(profile, 1, function(r) {
      sprintf("%d,%d,%d,%.5f", 1, 1, r[1], r[2])
    })
    lad_node <- newXMLNode("LAD-Profile",
                           attrs = c(type = "sparematrix-3D",
                                     dataI = dataI,
                                     dataJ = dataJ,
                                     zlayers = zlayers,
                                     defaultValue = "0.00000"),
                           .children = paste(lad_lines, collapse = "\n"))
    addChildren(plant_node, lad_node)
    addChildren(plant_node, newXMLNode("Season-Profile",
                                       paste(sprintf("%.5f", rep(1, 12)), collapse = ",")))
    addChildren(plant_node, newXMLNode("Blossom-Profile",
                                       paste(sprintf("%.5f", rep(0, 12)), collapse = ",")))
    addChildren(plant_node, newXMLNode("L-SystemBased", "0"))
    
    addChildren(root, plant_node)
  }
  
  saveXML(root, file = file_out, indent = TRUE, encoding = "UTF-8")
  message("✔ Envi-met PLANT3D (.pld) written to: ", normalizePath(file_out))
}

# Run export
export_lad_to_envimet3d(lad_df, file_out = xml_output_file)
