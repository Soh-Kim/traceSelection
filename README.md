
# What's traceSelection?

traceSelection is a R package specializes to calculate selection intensity for each cultivar, using pedigree and genotypic, and phenotypic information. This package is especially suitable for breeding data or germplasm collection data, where estimating selection intensity each individual underwent will provide with valueble information. The package has two main functions to calculate selection intensity and two helper functions, which is beneficial in handling pedigree file. 

### Functions provided
### main functions
**calcMarkerEff**
**calcSI**

### helper functions
**orderPed**
**selectPed**

## Installation

You can install the development version of traceSelection from [GitHub](https://github.com/) with:

**install using R package "pak"**
``` r
# if pak is not installed...
install.packages("pak")

pak::pak("Soh-Kim/selectTrace")
```

**install using R package "devtools"**
``` r
# if devtools is not installed...
install.packages("devtools")

devtools::install_github("Soh-Kim/selectTrace")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(traceSelection)
## basic example code
```

## Reference
