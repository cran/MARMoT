
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MARMoT

<!-- badges: start -->
<!-- badges: end -->

This R package contains the function to apply MARMoT balancing
technique; furthermore, it contains a function for computing the
Deloof’s approximation of the average rank (and also a parallelized
version) and a function to compute the Absolute Standardized Bias.

## Installation

You can install the development version of MARMoT like so:

``` r
devtools::install_github("AlbertoCalore/MARMoT")
```

## Example

This is a basic example of the function “MARMoT”:

``` r
library(MARMoT)
out = MARMoT(data = MARMoT_data, confounders = c("race", "age"), treatment = "hospital", n.cores = 1)
#> [1] "Absolute standardized bias - before balancing"
#>          Min 1st quartile       Median 3rd quartile          Max         Mean 
#>        0.010        2.875        6.859       12.027       50.700        8.667 
#>      Over 5%     Over 10% 
#>      214.000      118.000 
#> [1] "... Computing average rank ..."
#> [1] "... Balancing ..."
#> [1] "Absolute standardized bias - after balancing"
#>          Min 1st quartile       Median 3rd quartile          Max         Mean 
#>        0.132        0.138        0.412        0.570       43.302        1.055 
#>      Over 5%     Over 10% 
#>       23.000        4.000
```

### Foundings

This package was developed within the project funded by Next Generation
EU - “Age-It - Ageing well in an agening society” project (PE0000015),
National Recovery and Resilience Plan (NRRP) - PE8 - Mission 4, C2,
Intervention 1.3”. The views and opinions expressed are only those of
the authors and do not necessarily reflect those of the European Union
or the European Commission. Neither the European Union nor the European
Commission can be held responsible for them.
