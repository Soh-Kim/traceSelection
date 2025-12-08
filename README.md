
# What's traceSelection?

**traceSelection** is a R package specializes to calculate selection intensity for each individual, using pedigree and genotypic, and phenotypic information. This package is especially suitable for breeding data or germplasm collection data, where estimating selection intensity that  each individual underwent will provide with valueble information. The package has two main functions to calculate selection intensity and two helper functions, which is beneficial in handling pedigree file. 

### main functions
**・calcMarkerEff**  
Function to calculate marker effect using either Ridge, Lasso, BRR, BayesB, or Bayes C. The estimated marker effect is used in calcSI function to calculate selection intensity.    

**・calcSI**  
Function to calculate selection intensity using pedigree, phased-genotypic data as well as estimated marker effect.

### helper functions
**・orderPed**   
Function to sort pedigree information in ancestral order. The pedigree file are searched from ancestors in base population (whose parents are unknown), and individuals that can not be traced from ancestors are to be excluded.  

**・selectPed**    
Function to extract pedigree information of target individuals. It is usefull to exclude unintenresting individuals.


## Installation

You can install the development version of traceSelection from [GitHub](https://github.com/) by two different ways:

**install using R package "pak"**
``` r
# if pak is not installed...
install.packages("pak")

pak::pak("Soh-Kim/traceSeletion")
```

**install using R package "devtools"**
``` r
# if devtools is not installed...
install.packages("devtools")

devtools::install_github("Soh-Kim/traceSeletion")
```

## Example

### calcMarkerEff function
Suppose we have apply phenotypic data (color of fruit skin) and phased-genotypic data, labeled as "PhenoCol" and "genotype", respectively. 

``` r
> head(PhenoCol)
                FruitColor
GoldenDelicious       5.60
GrannySmith           6.24
Fuji                    NA
RedDelicious          7.26
PinkLady             13.32
Honeycrisp           18.56
```

```r
> genotype[1:6,1:6]
                m1001 m1003 m1021 m1042 m1103 m1256
GoldenDelicious "1|1" "0|1" "1|1" "1|1" "1|0" "0|1"
GrannySmith     "1|0" "0|1" "1|1" "1|0" "1|0" "1|0"
Fuji            "0|0" "0|1" "0|1" "1|0" "1|0" "0|1"
RedDelicious    "0|0" "0|1" "1|0" "0|0" "0|1" "1|0"
PinkLady        "0|0" "0|0" "1|1" "1|0" "1|1" "0|0"
Honeycrisp      "1|1" "0|1" "0|1" "1|0" "1|0" "1|1"
```

Now, marker effect is easily calculated by calcMarkerEff function as follows;

```r
> result <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "Ridge", Effect = "A" )

> result$coefDetermin
> result$rFit
> result$mseFit
> result$mEffect$Add
> result$mEffect$Dom
```

The different model can be easily implemented by changing the parameter of Model

```r
> result.Lasso <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "Lasso", Effect = "A" )
> result.BRR <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "BRR", Effect = "A" )
```

Not only additive effect but also dominance effect is calculated if parameter of Effect is replaced by "AD" (abbreviation of "Additive and Dominance"").
```r
> result <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "Ridge", Effect = "AD" )
> result$mEffect$Add
> result$mEffect$Dom
```

### calcSI function

### Visualization (Example)

## Reference
