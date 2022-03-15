# Species Portfolio Effects Dominate Seasonal Zooplankton Stabilization Within a Large Temperate Lake

Data for this analysis can be found online (https://data.ontario.ca/dataset/lake-simcoe-monitoring). 

# Data Prepping
Download the 1980-2019 csv files and open up the file named "Simcoe_Zooplankton&Bythotrephes.csv". Copy and paste the zooplankton sheet into a new excel file called "Simcoe_Zooplankton.csv". The column ZDATE in the excel file needs to be switched from GENERAL to SHORT DATE so that the dates in the ZDATE column read "YYYY/MM/DD". Save as .csv in appropriate R folder.

# R Scripts

R scripts are numbered (1-10) in the order that they are supposed to be used.

R scripts 1-4 represent the analysis of the data when 5 stations are included as well as code for statistical tests, final figures, and supplmemental figures.

R scripts 5-7 represent the analysis of the data when 3 stations are included as well as code for statistical tests and supplmemental figures.

R scripts 8-10 represent the analysis of the data when 7 stations are included as well as code for statistical tests and supplmemental figures.

The folder labelled "Figures" contains the raw figures as they are generated by the R scripts and may not have axes labelled the same as how they appear in the article.



# R Studio Session Info - sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] mefa_3.2-7      vegan_2.5-7     lattice_0.20-41 permute_0.9-5   reshape2_1.4.4  forcats_0.5.1  
 [7] stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4     readr_2.0.1     tidyr_1.1.3     tibble_3.1.4   
[13] ggplot2_3.3.5   tidyverse_1.3.1

loaded via a namespace (and not attached):
 [1] tinytex_0.26      tidyselect_1.1.0  xfun_0.25         splines_4.0.2     haven_2.3.1      
 [6] colorspace_1.4-1  vctrs_0.3.8       generics_0.1.0    mgcv_1.8-31       utf8_1.1.4       
[11] rlang_0.4.11      pillar_1.6.2      glue_1.4.2        withr_2.4.2       DBI_1.1.0        
[16] dbplyr_2.1.1      sessioninfo_1.1.1 modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.0  
[21] plyr_1.8.6        munsell_0.5.0     gtable_0.3.0      cellranger_1.1.0  rvest_1.0.1      
[26] tzdb_0.1.2        parallel_4.0.2    fansi_0.4.1       broom_0.7.9       Rcpp_1.0.5       
[31] backports_1.1.10  scales_1.1.1      jsonlite_1.7.2    fs_1.5.0          hms_1.1.0        
[36] stringi_1.5.3     grid_4.0.2        cli_3.0.1         tools_4.0.2       magrittr_2.0.1   
[41] cluster_2.1.0     crayon_1.4.1      pkgconfig_2.0.3   Matrix_1.2-18     ellipsis_0.3.2   
[46] MASS_7.3-51.6     xml2_1.3.2        reprex_2.0.1      lubridate_1.7.10  assertthat_0.2.1 
[51] httr_1.4.2        rstudioapi_0.13   R6_2.4.1          nlme_3.1-148      compiler_4.0.2   
