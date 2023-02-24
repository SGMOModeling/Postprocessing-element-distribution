# General Comments --------------------------------------------------------
#
# Read Groundwater zonal budget HDF file directly for element level statistics and visualization
# The script is intended for C2VSimFG QA/QC purpose.
#

#' Read Groundwater Zonal Budget HDF File
#'
#' Reads a Groundwater zonal budget HDF file directly for element level statistics and visualization. This script is intended for C2VSimFG QA/QC purpose.
#'
#' @param WYstart The beginning water year of simulation.
#' @param WYend The ending water year of simulation.
#' @param model_dir A character string giving the path of the model directory.
#' @param model_run A character string giving the name of the model run.
#' @param hdf_file_name A character string giving the name of the HDF file.
#' @param shp_file1 A character string giving the path of the C2VSimFG StreamReaches shapefile.
#' @param shp_element A character string giving the path of the C2VSimFG Elements shapefile.
#' @param cbPalette A character vector giving the color palette with grey.
#' @param cbbPalette A character vector giving the color palette with black.
#' @return For each layer, the script creates an HTML
#' 

# User Inputs -------------------------------------------------------------


{
  # Temporal Information
  WYstart <- 1974 # Beginning water year of simulation
  WYend <- 2015 # Ending water year of simulation

  model_dir <- "C:/Users/ghuang/Documents/c2vsimfg_version1.01/Results/"
  
  # model_dir <- "C:/Users/ghuang/Documents/GitHub/Postprocessing-element-distribution/data/"

  model_run <- "his_v1.01"

  # model_run <- "CG_v1.0"

  hdf_file_name <- "C2VSimFG_GW_ZBudget.hdf"
  # hdf_file_name <- "C2VSimCG_GW_ZBudget.hdf"
  
  shp_file1 <- "/data/GIS/C2VSimFG_StreamReaches.shp"
  
  shp_element <- "/data/GIS/C2VSimFG_Elements.shp"
  # shp_element <- "/data/GIS/C2VSimCG_Elements.shp"
  
  # These are color-blind-friendly palettes, one with gray, and one with black for plotting.

  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}
#
# Script Initialization ---------------------------------------------------


library(dplyr)
library(hydroTSM)
library(hydroGOF)
library(lubridate)
library(sf)
library(mapview)

#  The package rhdf5 must be installed.

# install.packages("BiocManager")
# BiocManager::install("rhdf5")

library(rhdf5)

# close all active hdf files in R work space to avoid reading error.
h5closeAll()

# hdf file name
hdf_file <- paste0(model_dir, hdf_file_name)

hdf_list <- h5ls(hdf_file) # List the content of an HDF5 file

hdf <- H5Fopen(hdf_file) # Open an object in an HDF5 file

# hydroTSM: Sequence of monthly dates between "1961-01-01" and "1961-12-31" ##

time_m <- mip(paste0(WYstart - 1, "-10-01"), paste0(WYend, "-09-30"))

name_list <- hdf_list$name

# from hdf file, "Layer_1"  check budget columns

col_nam <- name_list[23:51]

# check file attributes

h5validObjects()
a1 <- h5readAttributes(hdf, hdf_list$name[1])
a1
#
#
# all element values
#
column_index <- h5read(hdf, name = paste0("Attributes/", "FullDataNames"))
column_index

# There are 26 zonal budget columns excluding face flow and storage columns

column_loc <- h5read(hdf, name = paste0("Attributes/", "cLocationNames"))

# map view visualization





dsn <- paste0(getwd(), shp_element)
dsn_streaM <- paste0(getwd(), shp_file1)


nc_streams <- st_read(dsn_streaM, quiet = TRUE)
nc_element <- st_read(dsn, quiet = TRUE)

#  Loop over all 4 layers and 26 water budget columns

for (k in 1:26)
{
  # Layer 1
  layer_id <- 1
  for (layer_id in 1:4)
  {
    b1 <- trimws(column_loc[k + (layer_id - 1) * 26])

    # from ft3 to acre-ft

    element_value <- h5read(hdf, name = b1) / 43560.0 # compressed storage in different lengths

    el_index <- h5read(hdf, name = paste0("Attributes/", "Layer", layer_id, "_ElemDataColumns"))

    m1 <- length(element_value[, k])
    m2 <- length(el_index[, k])

    # total time steps number
    NTime <- (WYend - WYstart + 1) * 12

    # element values Initialization ---------------------------------------------------
    el_dx <- as.matrix(array(1:m2, c(m2, NTime))) * 0.0 #

    # Loop over all elements
    if (m1 == 0) m1 <- 1 # no non-zero values

    for (m in 1:m2)
    {{ if (el_index[m, k] > 0.0) {
      for (j in 1:NTime) {
        el_dx[m, j] <- element_value[el_index[m, k], j]
      }
    } }}

    {
      data1 <- as.data.frame(rowSums(el_dx))

      nyears <- (WYend - WYstart + 1)

      nc_element2 <- nc_element %>% mutate(z_af = data1$`rowSums(el_dx)` / nyears)



      # Mapview visualization

      x1 <- summary(nc_element2$z_af)
      x1
      title <- column_index[k]
      if (max(abs(x1)) > 0.01) {
        p1 <- mapview(nc_element2,
          zcol = "z_af", color = "white", color.region = cbPalette[3],
          alpha.regions = 0.5,
          at = seq(0.0 * x1[1] + 0.05, x1[6] * 1.2, (x1[6] - x1[1]) / 10),
          layer.name = paste0("Map for ", model_run, "_", title)
        )

        p2 <- p1 + mapview(nc_streams, zol = "Name")
        p2
        ## create standalone .html

        mapshot(p2, url = paste0(getwd(), "//output//groundwater//Map_", k, "_Layer", layer_id, ".htm"))

        #  rename map files for better read
        from_name <- paste0(getwd(), "//output//groundwater//Map_", k, "_Layer", layer_id, ".htm")
        to_name <- paste0(getwd(), "/output/groundwater/", model_run, "_'", trimws(column_index[k]), "_Layer", layer_id, ".htm")
        file.rename(from_name, to_name)
        # for csv ouput

        y <- cbind(c(1:m2), nc_element2$z_af)
        colnames(y) <- c("Element ID", paste0(b1, "(AF/year)"))

        write.csv(y, file = paste0(getwd(), "/output/groundwater/", model_run, "_'", trimws(column_index[k]), "_Layer", layer_id, ".csv"))
      }
    }
  }
}
write.csv(column_index, file = paste0(getwd(), "/output/groundwater/", model_run, "_map_list", ".csv"))
