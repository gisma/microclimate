#------------------------------------------------------------------------------
# Type: control script
# Name: 30_make_enviMet_simple_plants.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  derives tree hulls and the corresponding values for
#               tree species, height, LAD , albedo, transmissivity, root profile,
#               lucc classes and seasonality
# Data: point cloudand dsm and dtm as derived by 10_CHM_Catalog.R
#       sentinel 2 bands 2,3,4,8, classification of tree species  40_RS_high_resolution data
# Output: Envimet simple plant database file
# Copyright: Chris Reudenbach 2021, GPL (>= 3)
# git clone https://github.com/gisma/envimetR.git
#------------------------------------------------------------------------------

library(envimaR)
library(link2GI)
library(mapview)
library(terra)
library(rprojroot)
library(RColorBrewer)
root_folder = find_rstudio_root_file()




library(tidyr)
library(dplyr)
# 2 - define variables
#---------------------
fn = "5-25_MOF_rgb"
## ETRS89 / UTM zone 32N
epsg = 25832
# get viridris color palette
pal<-mapview::mapviewPalette("mapviewTopoColors")

otb  = link2GI::linkOTB(searchLocation = "~/apps/OTB-9.1.0-Linux/")
Sys.setenv(OTB_APPLICATION_PATH = paste0(dirname("/home/creu/apps/OTB-9.1.0-Linux/bin/"),"/lib/otb/applications"))

Sys.setenv(PATH = paste(otb$pathOTB, Sys.getenv("PATH"), sep = ":"))
           
# test area so called "sap flow halfmoon"
sapflow_ext=raster::extent(477500, 478218, 5631730, 5632500)
ts=data.frame(
  ID = 1:12,
  value = c(
    "agriculture",   # 1
    "alder",         # 2
    "ash",           # 3
    "beech",         # 4
    "douglas_fir",   # 5
    "larch",         # 6
    "oak",           # 7
    "pastures",      # 8
    "roads",         # 9
    "settlements",   #10
    "spruce",        #11
    "water"          #12
  )
)
min_tree_height = 2


library(terra)

replace_douglas_in_buche_eiche <- function(rast_input,
                                           window_size = 5,
                                           douglas_value = 5,
                                           target_values = c(4, 7),
                                           output_file = NULL) {
  # Sicherheitscheck
  if (!inherits(rast_input, "SpatRaster")) {
    stop("Input must be a terra::SpatRaster object.")
  }
  
  # Fokales Fenster (quadratisch)
  if (window_size %% 2 == 0) stop("window_size must be an odd number.")
  w <- matrix(1, nrow = window_size, ncol = window_size)
  
  # Lokaler Modus
  r_mode <- focal(rast_input, w = w, fun = modal, na.policy = "omit", na.rm = TRUE,progress = "text")
  
  # Maske: Douglasie
  is_douglas <- rast_input == douglas_value
  
  # Maske: Modus = Buche oder Eiche
  is_oak_beech_mode <- r_mode %in% target_values
  
  # Kombinierte Bedingung: Douglasie mit Eiche/Buche-Umgebung
  replace_mask <- is_douglas & is_oak_beech_mode
  
  # Ersetzen
  r_new <- rast_input
  r_new[replace_mask] <- r_mode[replace_mask]
  
  # Speichern (optional)
  if (!is.null(output_file)) {
    writeRaster(r_new, output_file, overwrite = TRUE)
  }
  
  return(r_new)
}

# get and filter the classification map
sapflow_species=readRDS("data/aerial/prediction_ffs_course2020_MOF_rgb.rds")
sapflow_species = raster::crop(sapflow_species,sapflow_ext)
raster::writeRaster(sapflow_species,"data/aerial/prediction_ffs.tif",progress="text",overwrite=TRUE)
#ClassificationMapRegularization majority filter
m1=mapview(sapflow_species, col.regions=brewer.pal( 12,"Paired"),maxpixels = 13821500)
cmr = parseOTBFunction("ClassificationMapRegularization",otb)
cmr$io.in = "data/aerial/prediction_ffs.tif"
cmr$io.out = "data/aerial/majority_out.tif"
cmr$progress = "true"
cmr$ip.radius = "1"



filter_treespecies = runOTB(cmr,gili = otb$pathOTB,quiet = FALSE,retRaster = TRUE)
#filter_treespecies=raster(paste0(envrmt$path_aerial,fn,"_majority_out.tif"))
f_tree_species=rast("data/aerial/majority_out.tif")
m2=mapview(f_tree_species, col.regions=brewer.pal( 12,"Paired"),maxpixels = 13821500)


# Eingabedatei laden (klassifiziertes Raster)
r <- rast("data/majority_in.tif")

# Zielauflösung in Meter (z. B. 1 m statt 0.2 m → Aggregationsfaktor = 5)
fact <- round(1 / res(r)[1])  # für 0.2 m Pixel ergibt das 5

# Aggregation per Median
r_agg <- aggregate(r, fact = fact, fun = median, na.rm = TRUE)
m3=mapview(r_agg, col.regions=brewer.pal( 12,"Paired"),maxpixels = 13821500)

# Ergebnis speichern
writeRaster(r_agg, "data/majority_in_1m.tif", overwrite = TRUE)
m1 + m2 +m3 


species_cleaned <- replace_douglas_in_buche_eiche(
  rast_input = r,
  window_size = 5,
  output_file = "data/aerial/majority_cleaned.tif"
)
## extract the values
species_ex = exactextractr::exact_extract(filter_treespecies, hulls_sf,  force_df = TRUE,
                                         include_cols = "treeID")
species_ex = dplyr::bind_rows(species_ex)

# calulate mode per tree
species_hulls = species_ex %>% group_by(treeID) %>%
  dplyr::reframe(species_median = median(value, na.rm=TRUE),
                   species_mode = modeest::mlv(value, method='mfv'))
species_hulls = inner_join(hulls_sf,species_hulls)
species_sf=species_hulls[,c("treeID","ZTOP","zmax","zmean","zsd","zskew","species_mode","area")]
st_write(species_hulls,file.path(envrmt$path_level1,"sapflow_tree_segments_multichm_dalponte2016_species.gpkg"),append=FALSE)
#species_hulls = st_read(file.path(envrmt$path_sapflow ,"sapflow_tree_segments_multichm_dalponte2016_species.gpkg"))

# calculate the LAD metrics for the derived trees
# review vertical LAI ditribution http://dx.doi.org/10.3390/rs12203457
lidR::writeLAS(sapflow_dalponte,file.path(envrmt$path_level1,"sapflow_dalponte.laz"))

VOXELS_LAD = lad.voxels(file.path(envrmt$path_level1,"sapflow_dalponte.laz"), grain.size = 2)
lad_profile = lad.profile(VOXELS_LAD,relative = F)


sapflow_dalponte@data$Z[sapflow_dalponte@data$Z < 0] = 0

maxZ = floor(max(sapflow_dalponte@data$Z))
func = formula(paste0("~pointsByZSlice(Z, ", maxZ, ")"))
t.binneds    = lidR::tree_metrics(sapflow_dalponte, func, res = 2,
                                  start = c(min(sapflow_dalponte@data$X), max(sapflow_dalponte@data$Y)))
lad_tree = tree_metrics(sapflow_dalponte, func = ~LAD(Z))
lad_tree = crown_metrics(sapflow_dalponte, ~metrics_lad(z = Z))
Fehler: Duplicated elements found. At least one of the metrics was not a number. Each metric should be a single number.


# Convert the result to a data frame for easier viewing
aggregated_lad_df <- data.frame(
  z_group = names(sum_lad_by_group),
  sum_lad = as.vector(sum_lad_by_group)
)

# GAP and transmittance metrics
# https://www.isprs.org/proceedings/xxxvi/3-W52/final_papers/Hopkinson_2007.pdf
gap_tree = tree_metrics(trees, func = ~lidR::gap_fraction_profile(Z))
gap_tree$gf = gap_tree$gf * 0.8
plot(gap_tree)
lt=st_drop_geometry(as(lad_tree,"sf"))
lad_trees = inner_join(lt,species_sf)

tmp = lt %>%
  group_by(treeID,z) %>%
  spread(z,lad,fill = 0,sep = "_")
lad_trees =  inner_join(tmp,lad_trees)


# read the setntinel data
albedo = raster(paste0(envrmt$path_sapflow,"S2B2A_20210613_108_sapflow_BOA_10_albedo.tif"))
lai = raster(paste0(envrmt$path_sapflow,"2021-06-13-00:00_2021-06-13-23:59_Sentinel-2_L1C_Custom_script.tiff"))

## extract the values
lai_ex = exactextractr::exact_extract(lai, hulls_sf,  force_df = TRUE,
                                         include_cols = "treeID")
lai_ex = dplyr::bind_rows(lai_ex)
lai_ex$coverage_fraction=NULL
names(lai_ex)=c("treeID","lai")
alb_ex = exactextractr::exact_extract(albedo, hulls_sf,  force_df = TRUE,
                                      include_cols = "treeID",)
alb_ex = dplyr::bind_rows(alb_ex)
alb_ex$coverage_fraction=NULL
names(alb_ex)=c("treeID","albedo")
l_trees =  inner_join(lai_ex,lad_trees)
a_trees = inner_join(alb_ex,l_trees)

# make mean of all unique treeIds
tree = a_trees %>% group_by(treeID) %>%
  mutate_all(.funs = mean) %>%
  distinct(.keep_all = TRUE)
tree$geom =NULL
t_lad=tree
tree_clust = tree[,c("treeID","albedo","lai","zmax","zmean","zskew","species_mode")]


tree_tmp = hulls_sf[ , names(hulls_sf) %in% c("treeID","geom")]
t=inner_join(tree_clust,tree_tmp)
tree_clust_sf=st_as_sf(t)
mapview(tree_clust_sf,zcol="albedo")
#t=st_drop_geometry(trees_sf)
header_height=seq(2.5,44.5,1)
#header_height=c(2.5, 3.5,   4.5,   5.5,   6.5,   7.5,   8.5,   9.5,   10.5,  11.5,  12.5,  13.5,  14.5,  15.5,  16.5,17.5,  18.5,  19.5,  20.5,  21.5,  22.5,  23.5,  24.5,  25.5,  26.5,  27.5,  28.5,  29.5,  30.5,  31.5,  32.5,33.5,  34.5,  35.5,  36.5,  37.5,  38.5,  39.5,  40.5,  41.5,  42.5,  43.5,  44.5)
lad= tree %>% distinct(across(contains(".5")))
norm1 = (lad[][2:length(lad)] - min(lad[][2:length(lad)] ))/ (max(lad[][2:length(lad)])-min(lad[][2:length(lad)]))
lad_spl = by(norm1, seq_len(nrow(norm1)), function(row) {try(smooth.spline(header_height, row,spar = 0.1))
  })

# for (i in seq(9000:9050)){
#   plot(lad_spl[[i]]$data$x,lad_spl[[i]]$data$y)
#   lines(lad_spl[[i]], lty = 2, col = "red")
#   lines( predict(lad_spl[[i]], xx), col = "green")
# }

xx = unique(sort(c(seq(4, 40, by = 4))))
lad_pred = lapply(lad_spl, function(x){predict(x, xx)[[2]]})
lad_pred = abs(as.data.frame(do.call(rbind, lad_pred)))
names(lad_pred)= paste0("lev_",xx)
lad_pred$treeID = lad$treeID


trees_all = inner_join(lad_pred,tree_clust)
trees_all =inner_join(trees_all,tree_tmp )
trees_all_sf=st_as_sf(trees_all)
st_write(trees_all_sf,file.path(envrmt$path_sapflow,"sapflow_tree_all_sf.gpkg"),append = FALSE)
trees_all_sf=st_read(file.path(envrmt$path_sapflow,"sapflow_tree_all_sf.gpkg"))
