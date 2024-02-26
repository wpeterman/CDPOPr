# CDPOPr

This package facilitates the running of basic CDPOP simulations from R and importing the results as `genind` objects. As inputs, it requires `SpatialPoints` and `raster` objects because the package `gdistance` is used to calculate effective distance between sample points.

## Installation

```{r}
devtools::install_github('wpeterman/CDPOPr',
                         dependencies = T,
                         build_vignettes = T)
```

The example in the vignette makes use of the `NLMR` package. If you plan to integrate `CDPOPr` with raster surface simulation, this can be a useful package. Follow instructions [HERE](https://github.com/ropensci/NLMR) to download `NLMR`.
