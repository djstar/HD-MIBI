# HD-MIBI
Code and Data accompanying the HD-MIBI paper


### 1. System requirements

This R notebook was written and executed on a custom built server with the following specifications.
2 x Intel(R) Xeon(R) CPU E5-2650
12 x 32GB Ram

The R sessioninfo is as follows:
R version 3.5.3 (2019-03-11)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] fastcluster_1.1.25  spatstat_1.58-2     rpart_4.1-13       
 [4] nlme_3.1-137        spatstat.data_1.4-0 gridExtra_2.3      
 [7] RColorBrewer_1.1-2  Rtsne_0.15          gplots_3.0.1.1     
[10] ggplot2_3.1.0       bioimagetools_1.1.3

loaded via a namespace (and not attached):
 [1] gtools_3.8.1          tidyselect_0.2.5      locfit_1.5-9.1       
 [4] xfun_0.4              purrr_0.3.0           splines_3.5.3        
 [7] lattice_0.20-38       colorspace_1.4-0      spatstat.utils_1.13-0
[10] htmltools_0.3.6       mgcv_1.8-27           rlang_0.4.0          
[13] pillar_1.3.1          glue_1.3.0            withr_2.1.2          
[16] EBImage_4.24.0        BiocGenerics_0.28.0   jpeg_0.1-8           
[19] plyr_1.8.4            munsell_0.5.0         gtable_0.2.0         
[22] htmlwidgets_1.3       caTools_1.17.1.2      knitr_1.21           
[25] Rcpp_1.0.2            tensor_1.5            KernSmooth_2.23-15   
[28] scales_1.0.0          gdata_2.18.0          abind_1.4-5          
[31] deldir_0.1-16         png_0.1-7             digest_0.6.18        
[34] tiff_0.1-5            dplyr_0.8.3           polyclip_1.9-1       
[37] grid_3.5.3            tools_3.5.3           bitops_1.0-6         
[40] magrittr_1.5          goftest_1.1-1         lazyeval_0.2.1       
[43] RCurl_1.95-4.11       tibble_2.0.1          crayon_1.3.4         
[46] pkgconfig_2.0.2       Matrix_1.2-15         assertthat_0.2.0     
[49] httr_1.4.0            rstudioapi_0.9.0      R6_2.3.0             
[52] fftwtools_0.9-8       compiler_3.5.3 



### 2. Installation

The following packages will be needed to execute the R notebook. The exact versions are listed above.
	- bioimagetools
	- ggplot2
	- gplots
	- Rtsne
	- RColorBrewer
	- gridExtra
	- spatstat
	- parallel
	- fastcluster

Additional custom functions called upon are provided in "srIBI_functions.R"

### 3. Instructions for use

Step by step instructions, along with data related to this manuscript, are provided in the R notebook: "Submission_Extract_Features_5_Cells_largervoxels.Rmd"

