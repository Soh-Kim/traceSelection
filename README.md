
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

### <ins>calcMarkerEff function</ins>
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

> result$model
[1] "Ridge"
> result$coefDetermin
[1] 0.5799931
> result$rFit
[1] 0.7793479
> result$mseFit
[1] 0.4980829

> head( result$yEst )
    GoldenDelicious         GrannySmith                Fuji        RedDelicious            PinkLady          Honeycrisp 
         -0.4366007           0.0471181          -0.2775965          -0.4318490          -0.6626601          -0.3674068

> head( result$mEffect$Add )
                      s0
m1001       0.0082345564
m1003      -0.0019489414
m1021      -0.0005759275
m1042       0.0084794079
m1103       0.0099916246
m1256       0.0090663859

> head( result$mEffect$Dom )
            s0
m1001       0
m1003       0
m1021       0
m1042       0
m1103       0
m1256       0
```
When "A" is applied to argument Effect, all the marker effects for dominance (result$mEffect$Dom) are estimated to be 0.

The different model can be easily implemented by changing the parameter of Model
```r
> result.Lasso <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "Lasso", Effect = "A" )
> result.BRR <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "BRR", Effect = "A" )
```

Not only additive effect but also dominance effect is calculated if parameter of Effect is replaced by "AD" (abbreviation of "Additive and Dominance"").
```r
> result <- calcMarkerEff( PhenoCol, genotype, phased = TRUE, Model = "Ridge", Effect = "AD" )

> head( result$mEffect$Add )
                      s0
m1001       0.002644945
m1003      -0.002283451
m1021      -0.001762202
m1042       0.005706774
m1103       0.006405869
m1256       0.005625770

> head( result$mEffect$Dom )
                     s0
m1001       0.009962970
m1003       0.005466636
m1021       0.006012601
m1042       0.001665153
m1103       0.001018533
m1256      -0.001178075
```

### <ins>calcSI function</ins>

In order to calculate selection intensity, we have to prepare tree different objects.
One is the phased genotype matrix, which is exactly same as above one (genotype).
Another one is the matrix summarizing chromosome, position, additive effect, and dominance effect, which looks like this;
```r
> head(mark)
  MarkerID Chr   Map  MarkerEff.a  MarkerEff.d
1    m1001   1 0.000  0.004766140  0.013356718
2    m1003   1 0.000 -0.002847033  0.007843717
3    m1021   1 0.000 -0.001987517  0.008748732
4    m1042   1 0.000  0.006850184  0.002274927
5    m1103   1 2.782  0.007411546  0.001073315
6    m1256   1 4.568  0.006699810 -0.001434129
```

Finally, pedigree information is required. For the unknown parents, any coding is acceptable including NA, "0", "unknown", and so on. 
```r
> head(ped)
              IID        Seed        Pollen
1 GoldenDelicious        <NA>          <NA>
2     GrannySmith        <NA>          <NA>
3            Fuji        <NA>          <NA>
4    RedDelicious        <NA>          <NA>
5        PinkLady        <NA>          <NA>
6             JS5        Fuji   GrannySmith
```

Together with these object, the selection intensity is calculated with function "calcSI" as follows;
```r
result <- calcSI( Marker = mark, Pedigree = ped, genoPhased = genome, nCore = 1 )
```
Here, nCore indicate the number of cores used in the calculation, and can be specified as you want.

The output is slightly different whether dominance effect was considered or not.
When dominance effect was considered (MarkerEff.d!=0), the output is a list containing tree element; Add, Dom, and Tot, and when dominance effect was not considered, the output is a list containing one element; Tot.  
Regardless of which element, each element contains data frame with 6 column, like this; 
```r
> head( result$Tot )
                  gEffect_s gEffect_p gEffect_mid  gEffect_o       Vg_s       Vg_p         Vg Select_Intens
GoldenDelicious          NA        NA          NA         NA         NA         NA         NA            NA
    GrannySmith          NA        NA          NA         NA         NA         NA         NA            NA
           Fuji          NA        NA          NA         NA         NA         NA         NA            NA
   RedDelicious          NA        NA          NA         NA         NA         NA         NA            NA
       PinkLady          NA        NA          NA         NA         NA         NA         NA            NA
            JS5  0.09927577 -0.954009  -0.3704914 -0.3169922 0.04111576 0.05608076 0.02607191     0.3313301
```
Here, each row corresponds to the information of each individual, and gEffect_s, gEffect_p are the genetic effect of seed and pollen parents, respectively. gEffect_mid is the expected value of progeny distribution, and gEffect_o is the genetic effect of the target individual. In addition, Vg is the variance of progeny distribution, and Vg_s, and Vg_p are the variance of progeny distributions attributable to seed and pollen parents, respectively. Finally, Select_Intens is the selection intensity.  
Because, these quantities are calculated based on the genome information of both target individual, and its parents, the row corresponding to the individuals with unknown parent becomes NA.  
  
Basically, the Tot is the result of all effects included, and Add and Dom are the results attributable to additive and dominance effects, respectively. More specifically, Add is the result obtained by setting all dominance effects to be zero, and in the same way, Dom is the result obtained by setting all additive effects to be zero. Because the variance of pedigree distribution in Dom (i.e. Vg) can not be partitioned into maternal and paternal contribution, Vg_s and Vg_p in Dom are always zero. 

### <ins>Visualization (Example)</ins>

## Reference
