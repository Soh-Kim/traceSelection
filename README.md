
# What's traceSelection?

**traceSelection** is a R package specializes to calculate selection intensity for each individual, using pedigree and genotypic, and phenotypic information. This package is especially suitable for breeding data or germplasm collection data, where estimating selection intensity that  each individual underwent will provide with valueble information. The package has two main functions to calculate selection intensity and two helper functions, which is beneficial in handling pedigree file. 

### Functions provided
### main functions
**calcMarkerEff**  
Function to calculate marker effect using either Ridge, Lasso, BRR, BayesB, or Bayes C. The estimated marker effect is used in calcSI function to calculate selection intensity.    
**calcSI**:  
Function to calculate selection intensity using pedigree, phased-genotypic data as well as estimated marker effect.

### helper functions
**orderPed**   
Function to sort pedigree information in ancestral order. The pedigree file are searched from ancestors in base population (whose parents are unknown), and individuals that can not be traced from ancestors are to be excluded.  

**selectPed**    
Function to extract pedigree information of target individuals. It is usefull to exclude unintenresting individuals.


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
