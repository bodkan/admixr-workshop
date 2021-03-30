Demo tutorial
================

## Introduction

Lets load all pre-requisites first. You should first try to install them
by running `install.packages(c("admixr", "tidyverse"))`. Note that this
make take a long time because there are a lots of C/C++ dependencies
which need to be compiled first (on a freshly installed system, this
could easily take 20 minutes).

``` r
library(admixr)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.0     ✓ dplyr   1.0.2
    ## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

## Demo \#1

#### EIGENSTRAT file format

Let’s direct *admixr* to the location of the EIGENSTRAT trio data files.
Note that this doesn’t load any data! It just creates “pointers” to
where the data files are located on your disk.s

``` r
prefix <- file.path(getwd(), "dataset")

data <- eigenstrat(prefix)

data
```

    ## EIGENSTRAT object
    ## =================
    ## components:
    ##   ind file: /Users/martin_petr/Google/postdoc/talks/2021-03-31_admixr/dataset.ind
    ##   snp file: /Users/martin_petr/Google/postdoc/talks/2021-03-31_admixr/dataset.snp
    ##   geno file: /Users/martin_petr/Google/postdoc/talks/2021-03-31_admixr/dataset.geno

Each EIGENSTRAT data set has three components: a trio of `ind/snp/geno`
files, which can be loaded with `read_*()` functions:

-   `ind` file:

``` r
ind <- read_ind(data)
ind
```

    ## # A tibble: 46 x 3
    ##    id           sex   label        
    ##    <chr>        <chr> <chr>        
    ##  1 Bichon       M     Bichon       
    ##  2 KK1          M     Kotias       
    ##  3 SATP         M     Satsurblia   
    ##  4 Motala12     M     Motala12     
    ##  5 I9030        M     Villabruna   
    ##  6 I0062        M     Vestonice16  
    ##  7 I0876        M     Kostenki14   
    ##  8 I0869.damage F     Ostuni1      
    ##  9 I0907.damage F     ElMiron      
    ## 10 I9050.damage F     AfontovaGora3
    ## # … with 36 more rows

-   `snp` file:

``` r
snp <- read_snp(data)
snp
```

    ## # A tibble: 2,144,502 x 6
    ##    id            chrom     gen    pos ref   alt  
    ##    <chr>         <chr>   <dbl>  <int> <chr> <chr>
    ##  1 Affx-34462167 1     0.00567 567137 T     C    
    ##  2 rs3094315     1     0.00753 752566 G     A    
    ##  3 rs12562034    1     0.00768 768448 A     G    
    ##  4 rs12124819    1     0.00776 776546 A     G    
    ##  5 1_777122      1     0.00777 777122 A     T    
    ##  6 1_832756      1     0.00833 832756 T     G    
    ##  7 rs28765502    1     0.00833 832918 T     C    
    ##  8 1_834832      1     0.00835 834832 G     C    
    ##  9 1_838555      1     0.00839 838555 C     A    
    ## 10 1_838931      1     0.00839 838931 A     C    
    ## # … with 2,144,492 more rows

``` r
nrow(snp) / 1e6 # how many millions of sites we have genotyped?
```

    ## [1] 2.144502

-   `geno` file:

``` r
geno <- read_geno(data)
geno
```

    ## # A tibble: 2,144,502 x 46
    ##    Bichon   KK1  SATP Motala12 I9030 I0062 I0876 I0869.damage I0907.damage
    ##     <int> <int> <int>    <int> <int> <int> <int>        <int>        <int>
    ##  1      0     0     0        0    NA    NA    NA           NA           NA
    ##  2      0     0     0        0     0     0     0            2            0
    ##  3     NA    NA    NA        0     0    NA     0           NA           NA
    ##  4      2     0     2        2    NA    NA     2           NA            2
    ##  5      0     0    NA        0    NA    NA    NA           NA           NA
    ##  6      2     0     2        0    NA     2    NA           NA           NA
    ##  7      2     2     2        0     2    NA     2           NA           NA
    ##  8      2     2     2        2    NA     0     2           NA            0
    ##  9      2     2     2        2    NA    NA     0           NA           NA
    ## 10      2     2    NA        2     2    NA     0           NA           NA
    ## # … with 2,144,492 more rows, and 37 more variables: I9050.damage <int>,
    ## #   Ranchot.damage <int>, Q116 <int>, Rochedane <int>, I1577 <int>, MA1 <int>,
    ## #   KO1 <int>, LaBrana1 <int>, I0061 <int>, Loschbour <int>, Ust_Ishim <int>,
    ## #   Stuttgart <int>, Estonian <int>, Mbuti <int>, Korean <int>, Yoruba <int>,
    ## #   Hezhen <int>, Orcadian <int>, Japanese <int>, Tujia <int>, Greek <int>,
    ## #   Russian <int>, Dai <int>, Han <int>, Dinka <int>, Czech <int>,
    ## #   Icelandic <int>, Uygur <int>, Burmese <int>, Tuscan <int>, Papuan <int>,
    ## #   Spanish <int>, Adygei <int>, French <int>, Chimp <int>, Vindija <int>,
    ## #   Altai <int>

We can see that the only genotypes present are reference and alternative
homozygotes. There are no heterozygous sites in our demo data because
all genomes have been pseudo-haplodized (at each heterozygous position,
a random allele was drawn from the two alleles present). In other cases,
allele types 0, 1, and 2 would be present.

``` r
table(geno$I0876) # Kostenki-14
```

    ## 
    ##       0       2 
    ##  461818 1285547

``` r
table(geno$Vindija)
```

    ## 
    ##       0       2 
    ##  597596 1361700

Every function in *admixr* works with the EIGENSTRAT trio and the
`eigenstrat()` function provided by the package makes this easier by
tracking where those files are and keeping them in one place in memory.

Note: in case you want to generate your own EIGENSTRAT data (from VCF
files, etc.) or if you want to somehow process or filter data you
already have, functions `write_snp/write_ind/write_geno` can be very
useful.

## Demo \#2

``` r
na.omit(geno[, c("Dinka", "Yoruba", "Vindija", "Chimp")])
```

    ## # A tibble: 1,897,864 x 4
    ##    Dinka Yoruba Vindija Chimp
    ##    <int>  <int>   <int> <int>
    ##  1     2      2       2     2
    ##  2     2      0       2     0
    ##  3     2      0       2     2
    ##  4     2      2       0     0
    ##  5     2      2       2     2
    ##  6     2      2       2     2
    ##  7     2      0       2     2
    ##  8     0      0       0     2
    ##  9     0      0       2     2
    ## 10     2      0       2     2
    ## # … with 1,897,854 more rows

``` r
f4_result1 <- f4(
  W = "Dinka", X = "Yoruba", Y = "Vindija", Z = "Chimp",
  data = data
)

f4_result1
```

    ## # A tibble: 1 x 10
    ##   W     X      Y       Z            f4   stderr Zscore  BABA  ABBA   nsnps
    ##   <chr> <chr>  <chr>   <chr>     <dbl>    <dbl>  <dbl> <dbl> <dbl>   <dbl>
    ## 1 Dinka Yoruba Vindija Chimp -0.000438 0.000281  -1.56 56211 57043 1897864

``` r
f4_result2 <- f4(
  W = "Dinka", X = "French", Y = "Vindija", Z = "Chimp",
  data = data
)

f4_result2
```

    ## # A tibble: 1 x 10
    ##   W     X      Y       Z           f4   stderr Zscore  BABA  ABBA   nsnps
    ##   <chr> <chr>  <chr>   <chr>    <dbl>    <dbl>  <dbl> <dbl> <dbl>   <dbl>
    ## 1 Dinka French Vindija Chimp -0.00220 0.000319  -6.87 55638 59805 1898330

### Homebrew *f*<sub>4</sub> statistic (and *D* statistic)

``` r
# subset to only four columns, for brevity
gt <- geno[, c("Dinka", "Yoruba", "Vindija", "Chimp")]

# detect two kinds of site patterns
baba <- with(gt, Dinka == Vindija & Yoruba == Chimp & Vindija != Chimp)
abba <- with(gt, Dinka == Chimp & Yoruba == Vindija & Vindija != Chimp)

# count the two site patterns
(total_baba <- sum(baba, na.rm = T))
```

    ## [1] 56211

``` r
(total_abba <- sum(abba, na.rm = T))
```

    ## [1] 57043

D statistic:

``` r
(total_baba - total_abba) / (total_baba + total_abba)
```

    ## [1] -0.007346319

*f*<sub>4</sub> statistic:

``` r
(total_baba - total_abba) / nrow(na.omit(gt))
```

    ## [1] -0.0004383876

Of course, this doesn’t give us Z-scores (p-values) and so we can’t
determine significance (although we could easily add that too) but it
servers as a useful demonstration that ADMIXTOOLS, *f*-statistics,
*admixr* and other tools are not black boxes built on high-level
techniques. In the end, they all boil down to this kind of calculation.

## Demo \#3 - *f*<sub>4</sub>-ratio statistic
