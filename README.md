# Postprocessing-element-distribution
 R script for IWFM output distribution of water budget by element 

R script which can read and process IWFM model (e.g. C2VSimFG or C2VSimCG, SVSim) element level water budget output in zonal budget HDF files, output to CSV files and visualize in web html pages as leaflet maps. 

The R script has been tested on: 

RStudio 2022.12.0+353 "Elsbeth Geranium" Release (7d165dcfc1b6d300eb247738db2c7076234f6ef0, 2022-12-03) for Windows
Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) RStudio/2022.12.0+353 Chrome/102.0.5005.167 Electron/19.1.3 Safari/537.36


sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rhdf5_2.42.0     mapview_2.11.0   sf_1.0-9         lubridate_1.9.0  timechange_0.1.1 hydroGOF_0.4-0   hydroTSM_0.6-0   xts_0.12.2      
 [9] zoo_1.8-11       dplyr_1.0.10    

loaded via a namespace (and not attached):
 [1] satellite_1.0.4         httr_1.4.4              webshot_0.5.4           RColorBrewer_1.1-3      rprojroot_2.0.3         R.cache_0.16.0         
 [7] tools_4.2.2             utf8_1.2.2              R6_2.5.1                KernSmooth_2.23-20      DBI_1.1.3               lazyeval_0.2.2         
[13] colorspace_2.0-3        rhdf5filters_1.10.0     raster_3.6-11           withr_2.5.0             sp_1.5-1                tidyselect_1.2.0       
[19] processx_3.8.0          leaflet_2.1.1           compiler_4.2.2          automap_1.0-16          leafem_0.2.0            cli_3.4.1              
[25] gstat_2.1-0             xml2_1.3.3              desc_1.4.2              plotly_4.10.1           scales_1.2.1            classInt_0.4-8         
[31] callr_3.7.3             proxy_0.4-27            stringr_1.5.0           digest_0.6.30           foreign_0.8-83          R.utils_2.12.2         
[37] rmarkdown_2.19          lintr_3.0.2             base64enc_0.1-3         pkgconfig_2.0.3         htmltools_0.5.4         styler_1.8.1           
[43] fastmap_1.1.0           htmlwidgets_1.6.0       rlang_1.0.6             ggthemes_4.2.4          rstudioapi_0.14         FNN_1.1.3.1            
[49] farver_2.1.1            generics_0.1.3          jsonlite_1.8.4          crosstalk_1.2.0         R.oo_1.25.0             magrittr_2.0.3         
[55] Rcpp_1.0.9              munsell_0.5.0           Rhdf5lib_1.20.0         fansi_1.0.3             R.methodsS3_1.8.2       lifecycle_1.0.3        
[61] terra_1.6-47            stringi_1.7.8           yaml_2.3.6              plyr_1.8.8              grid_4.2.2              maptools_1.1-5         
[67] crayon_1.5.2            lattice_0.20-45         cowplot_1.1.1           knitr_1.41              ps_1.7.2                pillar_1.8.1           
[73] spacetime_1.2-8         codetools_0.2-18        stats4_4.2.2            glue_1.6.2              evaluate_0.19           rex_1.2.1              
[79] leaflet.providers_1.9.0 data.table_1.14.6       remotes_2.4.2           png_0.1-8               vctrs_0.5.1             gtable_0.3.1           
[85] purrr_0.3.5             tidyr_1.2.1             reshape_0.8.9           assertthat_0.2.1        ggplot2_3.4.0           xfun_0.36              
[91] e1071_1.7-12            cyclocomp_1.1.0         viridisLite_0.4.1       class_7.3-20            tibble_3.1.8            intervals_0.15.2       
[97] units_0.8-1             ellipsis_0.3.2          xmlparsedata_1.0.5  