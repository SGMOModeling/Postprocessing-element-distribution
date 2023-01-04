#
# Read Groundwater zonal budget HDF file directly for element level statistics and visualization 
#
#
# User Inputs -------------------------------------------------------------

# Temporal Information

WYstart <- 1974    # Beginning water year of simulation
WYend <- 2015     # Ending water year of simulation
model_dir <- "C:/Users/ghuang/Documents/c2vsimfg_version1.01/Results/"

model_run <- "his_v1.01"
hdf_file_name <- "C2VSimFG_GW_ZBudget.hdf"

#These are color-blind-friendly palettes, one with gray, and one with black for plotting.

#The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#
# Script Initialization ---------------------------------------------------


library(dplyr)
library(hydroTSM)
library(hydroGOF)
library(lubridate)
library(sf)
library(mapview)

#  The package rhdf5 must be installed.

#install.packages("BiocManager")
#BiocManager::install("rhdf5")

library(rhdf5)

# close all active hdf files in R work space to avoid reading error.
h5closeAll()

# hdf file name
hdf_file=paste0(model_dir,hdf_file_name)

hdf_list <- h5ls(hdf_file)  #List the content of an HDF5 file

hdf = H5Fopen(hdf_file)  #Open an object in an HDF5 file

# hydroTSM: Sequence of monthly dates between "1961-01-01" and "1961-12-31" ##

time_m <- mip(paste0(WYstart-1,"-10-01"), paste0(WYend,"-09-30"))

name_list <- hdf_list$name

# from hdf file, "Layer_1"  budget columns

col_nam <- name_list[23:51]

# check file attributes

h5listIdentifier()
h5validObjects()
a1 <- h5readAttributes(hdf,hdf_list$name[1])
a1
#
#
# all element values
#

 
 column_index <- h5read(hdf, name =paste0("Attributes/","FullDataNames"))
 column_index 
 
 k=9
 k=17
 k=22
   
 {
 column_loc <- h5read(hdf, name =paste0("Attributes/","cLocationNames")) 
 
 b1 <- trimws(column_loc[k])
 b2 <- trimws(column_loc[k+2])

 # for GW storage
  if(k==1) b2 <-  trimws(column_loc[k+1])
 
 # from ft3 to acre-ft  
 
 element_value <- h5read(hdf, name =b1)/43560.0 # compressed storage in different lengths
 element_value1 <- h5read(hdf, name =b2)/43560.0  # compressed storage in different lengths
 
 
 
 el_index <- h5read(hdf, name =paste0("Attributes/","Layer1_ElemDataColumns"))
 
 m1 <- length(element_value[,k])
 m2 <- length(el_index[,k-1])
 
 # total time steps number
 NTime <- (WYend - WYstart+1)*12
 
 # element values Initialization ---------------------------------------------------
 el_dx <-as.matrix(array(1:m2,c(m2,NTime)))*0.0  # 
 
 # Loop over all elements
 
 for (m in 1:32537)
{
 
  {
    if (el_index[m, k] > 0.0) {
      for (j in 1:NTime) {
        el_dx[m, j] <- element_value[el_index[m, k], j] 
      }
    }
  }
}


# Loop over all elements
 
 sgn=1.0
 k1=k+2
if(k==1) {
  k1=2
  sgn=-1.0}
  
for (m in 1:32537)
{
  #  for(k in 1:26)
  {
    if (el_index[m, k1] > 0.0) {
      for (j in 1:504) {
        el_dx[m, j] <- el_dx[m, j] + element_value1[el_index[m, k1], j]*sgn
      }
    }
  }
}
 
 # spot check single element data
 
 j=14860  #Sub09
 j=4876   #Sub05
 j=9024   #sub05
 j=11470  #sub06
 j0=0
 


  # time series plot    
    
    {  
      j0=j0+1  
      plot(time_m, el_dx[j,], pch=16, col=cbPalette[2],xlab="Month",ylab=column_index[k], main=paste0("Figure ",j0," ", column_index[k], "  Element ID: ", j))
      #   
      lines(time_m, el_dx[j,])
    }
   
 
  }

# map view visualization


shp_file1 <- "C:/Users/ghuang/Documents/ArcGIS/C2VSimFG-V1_0_GIS/C2VSimFG-V1_0_GIS/Shapefiles/C2VSimFG_StreamReaches.shp"

shp_file2 <- "C:/Users/ghuang/Documents/ArcGIS/C2VSimFG-V1_0_GIS/C2VSimFG-V1_0_GIS/Shapefiles/C2VSim_Stream_Nodes.shp"

shp_element <- "C2VSimFG_Elements.shp"
shp_Node <- "C2VSimFG_Nodes.shp"

dsn <- paste0("C:/Users/ghuang/Documents/ArcGIS/C2VSimFG-V1_0_GIS/C2VSimFG-V1_0_GIS/Shapefiles/", shp_element)
dsn_node <- paste0("C:/Users/ghuang/Documents/ArcGIS/C2VSimFG-V1_0_GIS/C2VSimFG-V1_0_GIS/Shapefiles/", shp_Node)

nc_c2v <- st_read(shp_file1, quiet = TRUE)
nc_streamnodes <- st_read(shp_file2, quiet = TRUE)
nc_element <- st_read(dsn, quiet = TRUE)

nc_node <- st_read(dsn_node, quiet = FALSE)

{
data1 <- as.data.frame(rowSums(el_dx))


#value1 <-   data1$`rowSums(element_value)`/94.0/nc_element$Acres
#initializing 

nyears <- (WYend - WYstart+1)

nc_element2 <- nc_element %>% mutate(z_af=data1$`rowSums(el_dx)`/nyears)  #/nc_element$Acres)

write.csv(cbind(nc_element2$ElementID, nc_element2$z_af),file=paste0(getwd(),"/output/" ,column_index[k],".csv"))


#x1 <- data1 %>% summarise(mean)
# water year, oct to sep
#nc_element2$Acres <- data1$`rowSums(element_value)`/94.0/1000.0
x1 <- summary(nc_element2$z_af)
x1
title=column_index[k]
#title="small_watersheds_inflow"
p1 <- mapview(nc_element2,zcol="z_af", color = "white", color.region=cbPalette, 
              alpha.regions =0.5,
              at = seq(0.0*x1[1]+0.05, x1[6]*1.2, (x1[6]-x1[1])/10),
              layer.name=paste0("Map for ",model_run,"_",title))

p2 <- p1 + mapview(nc_c2v, zol="Name")
p2
## create standalone .html

mapshot(p2, url = paste0(getwd(),"/output/Map1.htm"))

#mapshot(p2, url = paste0(getwd(), "/map2",model_run,"_deep_perc",".htm"))
 
#plot(nc_element2)

}

