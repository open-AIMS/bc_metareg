AIMS <img src="logo.png" width = 180 alt="AIMS Logo" align="right" />
=========================================================================================

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![license](https://img.shields.io/badge/license-GPL--2-lightgrey.svg)](https://choosealicense.com/)
[![Ask Us Anything
!](https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg)](https://github.com/open-AIMS/bc_metareg/issues/new)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

This repository contains code and data needed to reproduce the article:

**Fulton CJ, Barneche DR, Davis K, Faubel C, Pascelli C, Vercelloni J, Wilson SK** (2025) Tracing blue carbon flows across diverse seascapes. *Global Change Biology*. doi: [10.1111/gcb.70420](https://doi.org/10.1111/gcb.70420)

## Instructions

**How to download this project for people not familiar with GitHub:**
on the project main page on GitHub, click on the green button `Code` and then
click on `Download ZIP`.

The data needs to be downloaded independently from another source (to added once archival completed), and stored in a folder called `data` inside the project root directory.

All analyses were done in [`R`](https://cran.r-project.org/) and uses [`Quarto`](https://quarto.org/) for rendering documents. Please use these links to have both software properly installed.

We use the [targets](https://github.com/ropensci/targets) R package to wrangle
the data, run the Bayesian models, and create all figures. First install
`targets`:

```r
install.packages("targets")
```

Next you need to open an R session with working directory set to the root of the project.

This routine uses some packages detailed under `_targets.R`, so make sure to install them before running anything.

Then, to reproduce all targets, run:

```r
source("run.R") # this may take 10--20 min depending on your machine specs.
```

The output figure will be saved to `output`.

## Reproducibility statement

```
R version 4.3.3 (2024-02-29)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

time zone: Australia/Perth
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] knitr_1.48         insight_0.19.10    viridisLite_0.4.2  ggdist_3.3.0       tidybayes_3.0.6    DHARMa_0.4.6       brms_2.21.0        Rcpp_1.0.13        rstan_2.32.6      
[10] StanHeaders_2.32.9 patchwork_1.2.0    readxl_1.4.3       lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.4       
[19] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1     svUnit_1.0.6         loo_2.7.0            TH.data_1.1-2        tensorA_0.36.2.1     timechange_0.2.0     estimability_1.4.1   lifecycle_1.0.4     
 [9] survival_3.5-8       magrittr_2.0.3       posterior_1.5.0      compiler_4.3.3       rlang_1.1.4          tools_4.3.3          utf8_1.2.4           bridgesampling_1.1-2
[17] pkgbuild_1.4.4       curl_5.2.1           abind_1.4-5          multcomp_1.4-25      withr_3.0.1          grid_4.3.3           stats4_4.3.3         fansi_1.0.6         
[25] xtable_1.8-4         colorspace_2.1-1     inline_0.3.19        emmeans_1.8.7        scales_1.3.0         MASS_7.3-60.0.1      cli_3.6.3            mvtnorm_1.2-6       
[33] generics_0.1.3       RcppParallel_5.1.8   tzdb_0.4.0           minqa_1.2.5          splines_4.3.3        bayesplot_1.11.1     parallel_4.3.3       cellranger_1.1.0    
[41] matrixStats_1.3.0    vctrs_0.6.5          V8_4.3.3             boot_1.3-29          Matrix_1.6-5         sandwich_3.0-2       jsonlite_1.8.8       hms_1.1.3           
[49] arrayhelpers_1.1-0   glue_1.7.0           nloptr_2.0.3         codetools_0.2-19     distributional_0.4.0 stringi_1.8.4        gtable_0.3.5         QuickJSR_1.3.1      
[57] lme4_1.1-35.3        munsell_0.5.1        pillar_1.9.0         Brobdingnag_1.2-9    R6_2.5.1             lattice_0.22-5       backports_1.5.0      rstantools_2.4.0    
[65] coda_0.19-4.1        gridExtra_2.3        nlme_3.1-164         checkmate_2.3.2      xfun_0.45            zoo_1.8-12           pkgconfig_2.0.3     
```

## License

This repository is provided by the authors under the MIT License
([MIT](https://opensource.org/license/mit/)).

## Bug reporting

Please [report any issues or bugs](https://github.com/open-AIMS/bc_metareg/issues).
