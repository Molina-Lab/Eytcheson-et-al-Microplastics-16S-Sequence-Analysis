Analysis
================
Sarah Brown
2023-04-28

## Load Packages and Import Data

Here, we import the .qza files exported from QIIME2 using the qqiime2R
package. The final object is a phyloseq object called MPphyseq.

``` r
#Load libraries
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.2.3

``` r
library(vegan)
```

    ## Warning: package 'vegan' was built under R version 4.2.3

    ## Loading required package: permute

    ## Warning: package 'permute' was built under R version 4.2.3

    ## Loading required package: lattice

    ## Warning: package 'lattice' was built under R version 4.2.3

    ## This is vegan 2.6-4

``` r
library(plyr)
```

    ## Warning: package 'plyr' was built under R version 4.2.3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.2.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(scales)
```

    ## Warning: package 'scales' was built under R version 4.2.3

``` r
library(grid)
library(reshape2)
```

    ## Warning: package 'reshape2' was built under R version 4.2.3

``` r
library(phyloseq)
library(picante)
```

    ## Warning: package 'picante' was built under R version 4.2.3

    ## Loading required package: ape

    ## Warning: package 'ape' was built under R version 4.2.3

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    ## Loading required package: nlme

    ## Warning: package 'nlme' was built under R version 4.2.3

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
library(tidyr)
```

    ## Warning: package 'tidyr' was built under R version 4.2.3

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:reshape2':
    ## 
    ##     smiths

``` r
library(viridis)
```

    ## Warning: package 'viridis' was built under R version 4.2.3

    ## Loading required package: viridisLite

    ## Warning: package 'viridisLite' was built under R version 4.2.3

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
library(qiime2R)
library(DESeq2)
```

    ## Warning: package 'DESeq2' was built under R version 4.2.2

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.2.2

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:plyr':
    ## 
    ##     rename

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:plyr':
    ## 
    ##     desc

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 4.2.2

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 4.2.2

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 4.2.3

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following object is masked from 'package:plyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

``` r
library(patchwork)
library(RColorBrewer)

set.seed(13289)

#Create a phyloseq object from the .qza files exported from qiime2 using 
#the qiime2R package
MPphyseq <- qza_to_phyloseq(
  features="Qiime output/dada2-table.qza",
  tree="Qiime output/rooted-tree.qza",
  taxonomy="Qiime output/MP-taxonomy.qza",
  metadata = "Qiime output/MP-metadata.tsv"
)
```

## Preparing the Data

``` r
#Check the rank names to make sure they are accurate
rank_names(MPphyseq)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

``` r
#Correct output:[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  
#                   "Genus"   "Species"

#Check sample variables
sample_variables(MPphyseq)
```

    ## [1] "effluent"      "week"          "polymer_type"  "bead_diameter"
    ## [5] "Channel"       "sample_type"   "particle_type"

``` r
#Remove chloroplast sequences and any contaminant sequences
MPphyseq <- subset_taxa(MPphyseq, Kingdom != "d__Eukaryota")
MPphyseq <- subset_taxa(MPphyseq, Kingdom != "d__Archaea")
MPphyseq <- subset_taxa(MPphyseq, Order != "o__Chloroplast")
MPphyseq <- subset_taxa(MPphyseq, Order != "Chloroplast")
MPphyseq <- subset_taxa(MPphyseq, Family != "f__Mitochondria")
MPphyseq <- subset_taxa(MPphyseq, Family != "Mitochondria")

#Check that contaminant sequences are removed (easiest to save as data frame and search)
taxtabl <- as.data.frame(tax_table(MPphyseq))

#At this point, blanks and positive controls should also be removed
MPphyseq = subset_samples(MPphyseq, effluent != "FB")
```

``` r
#Examining the total number of reads and ASV's
readsumsdf = data.frame(nreads = sort(taxa_sums(MPphyseq), TRUE), 
                        sorted = 1:ntaxa(MPphyseq), type = "ASVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(MPphyseq), 
                        TRUE), sorted = 1:nsamples(MPphyseq), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
```

![](Analysis_files/figure-gfm/asv%20plot-1.png)<!-- -->

``` r
#Rarefaction curve using vegan
#From: https://micca.readthedocs.io/en/latest/phyloseq.html
taxa_are_rows(MPphyseq)
```

    ## [1] TRUE

``` r
mat <- t(otu_table(MPphyseq))
class(mat) <- "matrix"
```

    ## Warning in class(mat) <- "matrix": Setting class(x) to "matrix" sets attribute
    ## to NULL; result will no longer be an S4 object

``` r
class(mat)
```

    ## [1] "matrix" "array"

``` r
mat <- as(t(otu_table(MPphyseq)), "matrix")
class(mat)
```

    ## [1] "matrix" "array"

``` r
raremax <- min(rowSums(mat))

system.time(rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE))
```

    ## Warning in rarecurve(mat, step = 100, sample = raremax, col = "blue", label =
    ## FALSE): most observed count data have counts 1, but smallest count is 2

![](Analysis_files/figure-gfm/rarefaction%20curve-1.png)<!-- -->

    ##    user  system elapsed 
    ##   56.04    4.19   61.44

``` r
#Transform to relative abundance. Save as new object.
MPphyseqRA = transform_sample_counts(MPphyseq, function(x){x / sum(x)})

#Re-label the variables; if you type in MPphyseq you will see that sample_data 
#is the matrix that holds the information that we want to change, so we need to 
#include sample_data in our code here

#Order factors
sample_data(MPphyseqRA)$effluent <- factor(sample_data(MPphyseqRA)$effluent, 
                                     levels = c("CON", "TWW"),
                                     labels = c("CON", "TWW"))

sample_data(MPphyseqRA)$polymer_type <- factor(sample_data(MPphyseqRA)$polymer_type, 
                                     levels = c("Glass", "HDPE", "LDPE", "PP", "PS", "water"),
                                     labels = c("Glass", "HDPE", "LDPE", "PP", "PS", "Water"))

sample_data(MPphyseqRA)$week <- factor(sample_data(MPphyseqRA)$week, 
                                               levels = c("0", "2", "6", "10"),
                                               labels = c("0", "2", "6", "10"))

sample_data(MPphyseqRA)$sample_type <- factor(sample_data(MPphyseqRA)$sample_type, 
                                       levels = c("Particle", "water"),
                                       labels = c("Particle", "Water"))

sample_data(MPphyseqRA)$particle_type <- factor(sample_data(MPphyseqRA)$particle_type, 
                                              levels = c("MP", "Glass", "water"),
                                              labels = c("MP", "Glass", "Water"))

#Subset samples into groups
sample_data(MPphyseqRA)
```

    ##                effluent week polymer_type bead_diameter Channel sample_type
    ## ESF21MP-100         CON   10         HDPE         3.175     7.1    Particle
    ## ESF21MP-101         TWW   10         HDPE         3.175     5.2    Particle
    ## ESF21MP-102         CON   10         LDPE         3.175     7.2    Particle
    ## ESF21MP-103         TWW   10           PP         3.175     5.2    Particle
    ## ESF21MP-104         CON   10           PP         3.175     7.2    Particle
    ## ESF21MP-105         TWW   10         HDPE         3.175     8.1    Particle
    ## ESF21MP-106         CON   10         LDPE         3.175     7.1    Particle
    ## ESF21MP-107         TWW   10         HDPE         3.175     5.1    Particle
    ## ESF21MP-108         CON   10         LDPE         3.175     8.2    Particle
    ## ESF21MP-109         CON   10           PP         3.175     8.2    Particle
    ## ESF21MP-110         TWW   10         HDPE         3.175     6.2    Particle
    ## ESF21MP-111         CON   10         HDPE         3.175     6.1    Particle
    ## ESF21MP-112         CON   10        Glass         3.175     8.2    Particle
    ## ESF21MP-113         TWW   10         LDPE         3.175     6.2    Particle
    ## ESF21MP-114         CON   10           PS         3.175     7.2    Particle
    ## ESF21MP-115         TWW   10        Glass         3.175     5.2    Particle
    ## ESF21MP-116         TWW   10           PP         3.175     6.2    Particle
    ## ESF21MP-117         CON   10        Glass         3.175     7.1    Particle
    ## ESF21MP-118         CON   10           PS         3.175     7.1    Particle
    ## ESF21MP-119         CON   10         HDPE         3.175     7.2    Particle
    ## ESF21MP-120         CON   10        Glass         3.175     6.1    Particle
    ## ESF21MP-121         CON   10         HDPE         3.175     8.2    Particle
    ## ESF21MP-122         CON   10         LDPE         3.175     6.1    Particle
    ## ESF21MP-123         TWW   10           PS         3.175     6.2    Particle
    ## ESF21MP-124         CON   10        Glass         3.175     7.2    Particle
    ## ESF21MP-125         CON   10           PS         3.175     8.2    Particle
    ## ESF21MP-126         TWW   10           PP         3.175     5.1    Particle
    ## ESF21MP-127         TWW   10           PS         3.175     5.2    Particle
    ## ESF21MP-128         TWW   10        Glass         3.175     6.2    Particle
    ## ESF21MP-129         CON   10           PP         3.175     6.1    Particle
    ## ESF21MP-13          TWW    2        Glass         3.175     5.1    Particle
    ## ESF21MP-130         TWW   10         LDPE         3.175     5.2    Particle
    ## ESF21MP-131         TWW   10         LDPE         3.175     5.1    Particle
    ## ESF21MP-132         TWW   10           PS         3.175     8.1    Particle
    ## ESF21MP-133         TWW   10        Glass         3.175     5.1    Particle
    ## ESF21MP-134         TWW   10        Glass         3.175     8.1    Particle
    ## ESF21MP-135         TWW   10         LDPE         3.175     8.1    Particle
    ## ESF21MP-14          CON    2           PP         3.175     8.2    Particle
    ## ESF21MP-15          CON    2           PS         3.175     7.2    Particle
    ## ESF21MP-16          TWW    2         LDPE         3.175     8.1    Particle
    ## ESF21MP-17          CON    2           PP         3.175     7.1    Particle
    ## ESF21MP-18          TWW    2         HDPE         3.175     6.2    Particle
    ## ESF21MP-19          CON    2           PP         3.175     7.2    Particle
    ## ESF21MP-20          CON    2           PS         3.175     8.2    Particle
    ## ESF21MP-21          TWW    2        Glass         3.175     6.2    Particle
    ## ESF21MP-22          TWW    2           PP         3.175     5.2    Particle
    ## ESF21MP-23          TWW    2           PP         3.175     5.1    Particle
    ## ESF21MP-24          CON    2         LDPE         3.175     7.1    Particle
    ## ESF21MP-25          CON    2           PP         3.175     6.1    Particle
    ## ESF21MP-26          TWW    2           PS         3.175     6.2    Particle
    ## ESF21MP-27          TWW    2         HDPE         3.175     5.1    Particle
    ## ESF21MP-28          CON    2           PS         3.175     7.1    Particle
    ## ESF21MP-29          TWW    2         LDPE         3.175     5.1    Particle
    ## ESF21MP-30          CON    2        Glass         3.175     7.1    Particle
    ## ESF21MP-31          CON    2         LDPE         3.175     8.2    Particle
    ## ESF21MP-33          CON    2         HDPE         3.175     6.1    Particle
    ## ESF21MP-34          CON    2         HDPE         3.175     7.2    Particle
    ## ESF21MP-35          TWW    2         LDPE         3.175     6.2    Particle
    ## ESF21MP-36          TWW    2           PP         3.175     8.1    Particle
    ## ESF21MP-37          TWW    2           PS         3.175     5.2    Particle
    ## ESF21MP-38          TWW    2         HDPE         3.175     8.1    Particle
    ## ESF21MP-40          CON    2         LDPE         3.175     7.2    Particle
    ## ESF21MP-41          TWW    2           PS         3.175     8.1    Particle
    ## ESF21MP-42          TWW    2         LDPE         3.175     5.2    Particle
    ## ESF21MP-43          TWW    2           PP         3.175     6.2    Particle
    ## ESF21MP-44          TWW    2        Glass         3.175     8.1    Particle
    ## ESF21MP-46          CON    2         HDPE         3.175     8.2    Particle
    ## ESF21MP-47          CON    2           PS         3.175     6.1    Particle
    ## ESF21MP-48          TWW    2        Glass         3.175     5.2    Particle
    ## ESF21MP-49          CON    2         LDPE         3.175     6.1    Particle
    ## ESF21MP-5-1-0       TWW    0        Water         0.000     5.1       Water
    ## ESF21MP-5-1-10      TWW   10        Water         0.000     5.1       Water
    ## ESF21MP-5-1-2       TWW    2        Water         0.000     5.1       Water
    ## ESF21MP-5-1-6       TWW    6        Water         0.000     5.1       Water
    ## ESF21MP-5-2-0       TWW    0        Water         0.000     5.2       Water
    ## ESF21MP-5-2-10      TWW   10        Water         0.000     5.2       Water
    ## ESF21MP-5-2-2       TWW    2        Water         0.000     5.2       Water
    ## ESF21MP-5-2-6       TWW    6        Water         0.000     5.2       Water
    ## ESF21MP-50          TWW    2         HDPE         3.175     5.2    Particle
    ## ESF21MP-51          CON    2         HDPE         3.175     7.1    Particle
    ## ESF21MP-52          CON    2        Glass         3.175     7.2    Particle
    ## ESF21MP-53          TWW    2           PS         3.175     5.1    Particle
    ## ESF21MP-54          CON    2        Glass         3.175     8.2    Particle
    ## ESF21MP-55          CON    2        Glass         3.175     6.1    Particle
    ## ESF21MP-56          TWW    6           PS         3.175     5.1    Particle
    ## ESF21MP-57          CON    6         LDPE         3.175     8.2    Particle
    ## ESF21MP-58          CON    6         HDPE         3.175     8.2    Particle
    ## ESF21MP-59          CON    6           PS         3.175     7.1    Particle
    ## ESF21MP-6-1-0       CON    0        Water         0.000     6.1       Water
    ## ESF21MP-6-1-10      CON   10        Water         0.000     6.1       Water
    ## ESF21MP-6-1-2       CON    2        Water         0.000     6.1       Water
    ## ESF21MP-6-1-6       CON    6        Water         0.000     6.1       Water
    ## ESF21MP-6-2-0       TWW    0        Water         0.000     6.2       Water
    ## ESF21MP-6-2-10      TWW   10        Water         0.000     6.2       Water
    ## ESF21MP-6-2-2       TWW    2        Water         0.000     6.2       Water
    ## ESF21MP-6-2-6       TWW    6        Water         0.000     6.2       Water
    ## ESF21MP-60          CON    6           PP         3.175     6.1    Particle
    ## ESF21MP-61          TWW    6        Glass         3.175     8.1    Particle
    ## ESF21MP-62          TWW    6           PP         3.175     8.1    Particle
    ## ESF21MP-63          TWW    6           PP         3.175     5.2    Particle
    ## ESF21MP-64          TWW    6        Glass         3.175     5.2    Particle
    ## ESF21MP-65          TWW    6        Glass         3.175     6.2    Particle
    ## ESF21MP-66          CON    6         HDPE         3.175     6.1    Particle
    ## ESF21MP-67          CON    6        Glass         3.175     8.2    Particle
    ## ESF21MP-68          TWW    6         LDPE         3.175     6.2    Particle
    ## ESF21MP-69          TWW    6        Glass         3.175     5.1    Particle
    ## ESF21MP-7-1-0       CON    0        Water         0.000     7.1       Water
    ## ESF21MP-7-1-10      CON   10        Water         0.000     7.1       Water
    ## ESF21MP-7-1-2       CON    2        Water         0.000     7.1       Water
    ## ESF21MP-7-1-6       CON    6        Water         0.000     7.1       Water
    ## ESF21MP-7-2-0       CON    0        Water         0.000     7.2       Water
    ## ESF21MP-7-2-10      CON   10        Water         0.000     7.2       Water
    ## ESF21MP-7-2-2       CON    2        Water         0.000     7.2       Water
    ## ESF21MP-7-2-6       CON    6        Water         0.000     7.2       Water
    ## ESF21MP-70          CON    6           PP         3.175     7.1    Particle
    ## ESF21MP-71          CON    6           PP         3.175     8.2    Particle
    ## ESF21MP-72          CON    6        Glass         3.175     7.2    Particle
    ## ESF21MP-73          TWW    6           PP         3.175     6.2    Particle
    ## ESF21MP-74          TWW    6         HDPE         3.175     5.2    Particle
    ## ESF21MP-75          CON    6         LDPE         3.175     7.2    Particle
    ## ESF21MP-76          CON    6         LDPE         3.175     6.1    Particle
    ## ESF21MP-77          TWW    6         LDPE         3.175     8.1    Particle
    ## ESF21MP-78          TWW    6         LDPE         3.175     5.2    Particle
    ## ESF21MP-79          CON    6         HDPE         3.175     7.2    Particle
    ## ESF21MP-8-1-0       TWW    0        Water         0.000     8.1       Water
    ## ESF21MP-8-1-10      TWW   10        Water         0.000     8.1       Water
    ## ESF21MP-8-1-2       TWW    2        Water         0.000     8.1       Water
    ## ESF21MP-8-1-6       TWW    6        Water         0.000     8.1       Water
    ## ESF21MP-8-2-0       CON    0        Water         0.000     8.2       Water
    ## ESF21MP-8-2-10      CON   10        Water         0.000     8.2       Water
    ## ESF21MP-8-2-2       CON    2        Water         0.000     8.2       Water
    ## ESF21MP-8-2-6       CON    6        Water         0.000     8.2       Water
    ## ESF21MP-80          TWW    6           PS         3.175     6.2    Particle
    ## ESF21MP-81          CON    6           PP         3.175     7.2    Particle
    ## ESF21MP-82          CON    6         LDPE         3.175     7.1    Particle
    ## ESF21MP-83          TWW    6         HDPE         3.175     6.2    Particle
    ## ESF21MP-84          CON    6           PS         3.175     6.1    Particle
    ## ESF21MP-86          TWW    6         HDPE         3.175     8.1    Particle
    ## ESF21MP-87          CON    6        Glass         3.175     6.1    Particle
    ## ESF21MP-88          TWW    6         HDPE         3.175     5.1    Particle
    ## ESF21MP-89          CON    6         HDPE         3.175     7.1    Particle
    ## ESF21MP-90          TWW    6           PS         3.175     5.2    Particle
    ## ESF21MP-91          TWW    6         LDPE         3.175     5.1    Particle
    ## ESF21MP-92          CON    6           PS         3.175     8.2    Particle
    ## ESF21MP-93          TWW    6           PP         3.175     5.1    Particle
    ## ESF21MP-94          TWW    6           PS         3.175     8.1    Particle
    ## ESF21MP-95          CON    6        Glass         3.175     7.1    Particle
    ## ESF21MP-96          CON   10           PP         3.175     7.1    Particle
    ## ESF21MP-97          TWW   10           PP         3.175     8.1    Particle
    ## ESF21MP-98          TWW   10           PS         3.175     5.1    Particle
    ## ESF21MP-99          CON   10           PS         3.175     6.1    Particle
    ##                particle_type
    ## ESF21MP-100               MP
    ## ESF21MP-101               MP
    ## ESF21MP-102               MP
    ## ESF21MP-103               MP
    ## ESF21MP-104               MP
    ## ESF21MP-105               MP
    ## ESF21MP-106               MP
    ## ESF21MP-107               MP
    ## ESF21MP-108               MP
    ## ESF21MP-109               MP
    ## ESF21MP-110               MP
    ## ESF21MP-111               MP
    ## ESF21MP-112            Glass
    ## ESF21MP-113               MP
    ## ESF21MP-114               MP
    ## ESF21MP-115            Glass
    ## ESF21MP-116               MP
    ## ESF21MP-117            Glass
    ## ESF21MP-118               MP
    ## ESF21MP-119               MP
    ## ESF21MP-120            Glass
    ## ESF21MP-121               MP
    ## ESF21MP-122               MP
    ## ESF21MP-123               MP
    ## ESF21MP-124            Glass
    ## ESF21MP-125               MP
    ## ESF21MP-126               MP
    ## ESF21MP-127               MP
    ## ESF21MP-128            Glass
    ## ESF21MP-129               MP
    ## ESF21MP-13             Glass
    ## ESF21MP-130               MP
    ## ESF21MP-131               MP
    ## ESF21MP-132               MP
    ## ESF21MP-133            Glass
    ## ESF21MP-134            Glass
    ## ESF21MP-135               MP
    ## ESF21MP-14                MP
    ## ESF21MP-15                MP
    ## ESF21MP-16                MP
    ## ESF21MP-17                MP
    ## ESF21MP-18                MP
    ## ESF21MP-19                MP
    ## ESF21MP-20                MP
    ## ESF21MP-21             Glass
    ## ESF21MP-22                MP
    ## ESF21MP-23                MP
    ## ESF21MP-24                MP
    ## ESF21MP-25                MP
    ## ESF21MP-26                MP
    ## ESF21MP-27                MP
    ## ESF21MP-28                MP
    ## ESF21MP-29                MP
    ## ESF21MP-30             Glass
    ## ESF21MP-31                MP
    ## ESF21MP-33                MP
    ## ESF21MP-34                MP
    ## ESF21MP-35                MP
    ## ESF21MP-36                MP
    ## ESF21MP-37                MP
    ## ESF21MP-38                MP
    ## ESF21MP-40                MP
    ## ESF21MP-41                MP
    ## ESF21MP-42                MP
    ## ESF21MP-43                MP
    ## ESF21MP-44             Glass
    ## ESF21MP-46                MP
    ## ESF21MP-47                MP
    ## ESF21MP-48             Glass
    ## ESF21MP-49                MP
    ## ESF21MP-5-1-0          Water
    ## ESF21MP-5-1-10         Water
    ## ESF21MP-5-1-2          Water
    ## ESF21MP-5-1-6          Water
    ## ESF21MP-5-2-0          Water
    ## ESF21MP-5-2-10         Water
    ## ESF21MP-5-2-2          Water
    ## ESF21MP-5-2-6          Water
    ## ESF21MP-50                MP
    ## ESF21MP-51                MP
    ## ESF21MP-52             Glass
    ## ESF21MP-53                MP
    ## ESF21MP-54             Glass
    ## ESF21MP-55             Glass
    ## ESF21MP-56                MP
    ## ESF21MP-57                MP
    ## ESF21MP-58                MP
    ## ESF21MP-59                MP
    ## ESF21MP-6-1-0          Water
    ## ESF21MP-6-1-10         Water
    ## ESF21MP-6-1-2          Water
    ## ESF21MP-6-1-6          Water
    ## ESF21MP-6-2-0          Water
    ## ESF21MP-6-2-10         Water
    ## ESF21MP-6-2-2          Water
    ## ESF21MP-6-2-6          Water
    ## ESF21MP-60                MP
    ## ESF21MP-61             Glass
    ## ESF21MP-62                MP
    ## ESF21MP-63                MP
    ## ESF21MP-64             Glass
    ## ESF21MP-65             Glass
    ## ESF21MP-66                MP
    ## ESF21MP-67             Glass
    ## ESF21MP-68             Glass
    ## ESF21MP-69             Glass
    ## ESF21MP-7-1-0          Water
    ## ESF21MP-7-1-10         Water
    ## ESF21MP-7-1-2          Water
    ## ESF21MP-7-1-6          Water
    ## ESF21MP-7-2-0          Water
    ## ESF21MP-7-2-10         Water
    ## ESF21MP-7-2-2          Water
    ## ESF21MP-7-2-6          Water
    ## ESF21MP-70                MP
    ## ESF21MP-71                MP
    ## ESF21MP-72             Glass
    ## ESF21MP-73                MP
    ## ESF21MP-74                MP
    ## ESF21MP-75                MP
    ## ESF21MP-76                MP
    ## ESF21MP-77                MP
    ## ESF21MP-78                MP
    ## ESF21MP-79                MP
    ## ESF21MP-8-1-0          Water
    ## ESF21MP-8-1-10         Water
    ## ESF21MP-8-1-2          Water
    ## ESF21MP-8-1-6          Water
    ## ESF21MP-8-2-0          Water
    ## ESF21MP-8-2-10         Water
    ## ESF21MP-8-2-2          Water
    ## ESF21MP-8-2-6          Water
    ## ESF21MP-80                MP
    ## ESF21MP-81                MP
    ## ESF21MP-82                MP
    ## ESF21MP-83                MP
    ## ESF21MP-84                MP
    ## ESF21MP-86                MP
    ## ESF21MP-87             Glass
    ## ESF21MP-88                MP
    ## ESF21MP-89                MP
    ## ESF21MP-90                MP
    ## ESF21MP-91                MP
    ## ESF21MP-92                MP
    ## ESF21MP-93                MP
    ## ESF21MP-94                MP
    ## ESF21MP-95             Glass
    ## ESF21MP-96                MP
    ## ESF21MP-97                MP
    ## ESF21MP-98                MP
    ## ESF21MP-99                MP

``` r
mp.con <- subset_samples(MPphyseqRA, effluent == "CON")
sample_data(mp.con)
```

    ##                effluent week polymer_type bead_diameter Channel sample_type
    ## ESF21MP-100         CON   10         HDPE         3.175     7.1    Particle
    ## ESF21MP-102         CON   10         LDPE         3.175     7.2    Particle
    ## ESF21MP-104         CON   10           PP         3.175     7.2    Particle
    ## ESF21MP-106         CON   10         LDPE         3.175     7.1    Particle
    ## ESF21MP-108         CON   10         LDPE         3.175     8.2    Particle
    ## ESF21MP-109         CON   10           PP         3.175     8.2    Particle
    ## ESF21MP-111         CON   10         HDPE         3.175     6.1    Particle
    ## ESF21MP-112         CON   10        Glass         3.175     8.2    Particle
    ## ESF21MP-114         CON   10           PS         3.175     7.2    Particle
    ## ESF21MP-117         CON   10        Glass         3.175     7.1    Particle
    ## ESF21MP-118         CON   10           PS         3.175     7.1    Particle
    ## ESF21MP-119         CON   10         HDPE         3.175     7.2    Particle
    ## ESF21MP-120         CON   10        Glass         3.175     6.1    Particle
    ## ESF21MP-121         CON   10         HDPE         3.175     8.2    Particle
    ## ESF21MP-122         CON   10         LDPE         3.175     6.1    Particle
    ## ESF21MP-124         CON   10        Glass         3.175     7.2    Particle
    ## ESF21MP-125         CON   10           PS         3.175     8.2    Particle
    ## ESF21MP-129         CON   10           PP         3.175     6.1    Particle
    ## ESF21MP-14          CON    2           PP         3.175     8.2    Particle
    ## ESF21MP-15          CON    2           PS         3.175     7.2    Particle
    ## ESF21MP-17          CON    2           PP         3.175     7.1    Particle
    ## ESF21MP-19          CON    2           PP         3.175     7.2    Particle
    ## ESF21MP-20          CON    2           PS         3.175     8.2    Particle
    ## ESF21MP-24          CON    2         LDPE         3.175     7.1    Particle
    ## ESF21MP-25          CON    2           PP         3.175     6.1    Particle
    ## ESF21MP-28          CON    2           PS         3.175     7.1    Particle
    ## ESF21MP-30          CON    2        Glass         3.175     7.1    Particle
    ## ESF21MP-31          CON    2         LDPE         3.175     8.2    Particle
    ## ESF21MP-33          CON    2         HDPE         3.175     6.1    Particle
    ## ESF21MP-34          CON    2         HDPE         3.175     7.2    Particle
    ## ESF21MP-40          CON    2         LDPE         3.175     7.2    Particle
    ## ESF21MP-46          CON    2         HDPE         3.175     8.2    Particle
    ## ESF21MP-47          CON    2           PS         3.175     6.1    Particle
    ## ESF21MP-49          CON    2         LDPE         3.175     6.1    Particle
    ## ESF21MP-51          CON    2         HDPE         3.175     7.1    Particle
    ## ESF21MP-52          CON    2        Glass         3.175     7.2    Particle
    ## ESF21MP-54          CON    2        Glass         3.175     8.2    Particle
    ## ESF21MP-55          CON    2        Glass         3.175     6.1    Particle
    ## ESF21MP-57          CON    6         LDPE         3.175     8.2    Particle
    ## ESF21MP-58          CON    6         HDPE         3.175     8.2    Particle
    ## ESF21MP-59          CON    6           PS         3.175     7.1    Particle
    ## ESF21MP-6-1-0       CON    0        Water         0.000     6.1       Water
    ## ESF21MP-6-1-10      CON   10        Water         0.000     6.1       Water
    ## ESF21MP-6-1-2       CON    2        Water         0.000     6.1       Water
    ## ESF21MP-6-1-6       CON    6        Water         0.000     6.1       Water
    ## ESF21MP-60          CON    6           PP         3.175     6.1    Particle
    ## ESF21MP-66          CON    6         HDPE         3.175     6.1    Particle
    ## ESF21MP-67          CON    6        Glass         3.175     8.2    Particle
    ## ESF21MP-7-1-0       CON    0        Water         0.000     7.1       Water
    ## ESF21MP-7-1-10      CON   10        Water         0.000     7.1       Water
    ## ESF21MP-7-1-2       CON    2        Water         0.000     7.1       Water
    ## ESF21MP-7-1-6       CON    6        Water         0.000     7.1       Water
    ## ESF21MP-7-2-0       CON    0        Water         0.000     7.2       Water
    ## ESF21MP-7-2-10      CON   10        Water         0.000     7.2       Water
    ## ESF21MP-7-2-2       CON    2        Water         0.000     7.2       Water
    ## ESF21MP-7-2-6       CON    6        Water         0.000     7.2       Water
    ## ESF21MP-70          CON    6           PP         3.175     7.1    Particle
    ## ESF21MP-71          CON    6           PP         3.175     8.2    Particle
    ## ESF21MP-72          CON    6        Glass         3.175     7.2    Particle
    ## ESF21MP-75          CON    6         LDPE         3.175     7.2    Particle
    ## ESF21MP-76          CON    6         LDPE         3.175     6.1    Particle
    ## ESF21MP-79          CON    6         HDPE         3.175     7.2    Particle
    ## ESF21MP-8-2-0       CON    0        Water         0.000     8.2       Water
    ## ESF21MP-8-2-10      CON   10        Water         0.000     8.2       Water
    ## ESF21MP-8-2-2       CON    2        Water         0.000     8.2       Water
    ## ESF21MP-8-2-6       CON    6        Water         0.000     8.2       Water
    ## ESF21MP-81          CON    6           PP         3.175     7.2    Particle
    ## ESF21MP-82          CON    6         LDPE         3.175     7.1    Particle
    ## ESF21MP-84          CON    6           PS         3.175     6.1    Particle
    ## ESF21MP-87          CON    6        Glass         3.175     6.1    Particle
    ## ESF21MP-89          CON    6         HDPE         3.175     7.1    Particle
    ## ESF21MP-92          CON    6           PS         3.175     8.2    Particle
    ## ESF21MP-95          CON    6        Glass         3.175     7.1    Particle
    ## ESF21MP-96          CON   10           PP         3.175     7.1    Particle
    ## ESF21MP-99          CON   10           PS         3.175     6.1    Particle
    ##                particle_type
    ## ESF21MP-100               MP
    ## ESF21MP-102               MP
    ## ESF21MP-104               MP
    ## ESF21MP-106               MP
    ## ESF21MP-108               MP
    ## ESF21MP-109               MP
    ## ESF21MP-111               MP
    ## ESF21MP-112            Glass
    ## ESF21MP-114               MP
    ## ESF21MP-117            Glass
    ## ESF21MP-118               MP
    ## ESF21MP-119               MP
    ## ESF21MP-120            Glass
    ## ESF21MP-121               MP
    ## ESF21MP-122               MP
    ## ESF21MP-124            Glass
    ## ESF21MP-125               MP
    ## ESF21MP-129               MP
    ## ESF21MP-14                MP
    ## ESF21MP-15                MP
    ## ESF21MP-17                MP
    ## ESF21MP-19                MP
    ## ESF21MP-20                MP
    ## ESF21MP-24                MP
    ## ESF21MP-25                MP
    ## ESF21MP-28                MP
    ## ESF21MP-30             Glass
    ## ESF21MP-31                MP
    ## ESF21MP-33                MP
    ## ESF21MP-34                MP
    ## ESF21MP-40                MP
    ## ESF21MP-46                MP
    ## ESF21MP-47                MP
    ## ESF21MP-49                MP
    ## ESF21MP-51                MP
    ## ESF21MP-52             Glass
    ## ESF21MP-54             Glass
    ## ESF21MP-55             Glass
    ## ESF21MP-57                MP
    ## ESF21MP-58                MP
    ## ESF21MP-59                MP
    ## ESF21MP-6-1-0          Water
    ## ESF21MP-6-1-10         Water
    ## ESF21MP-6-1-2          Water
    ## ESF21MP-6-1-6          Water
    ## ESF21MP-60                MP
    ## ESF21MP-66                MP
    ## ESF21MP-67             Glass
    ## ESF21MP-7-1-0          Water
    ## ESF21MP-7-1-10         Water
    ## ESF21MP-7-1-2          Water
    ## ESF21MP-7-1-6          Water
    ## ESF21MP-7-2-0          Water
    ## ESF21MP-7-2-10         Water
    ## ESF21MP-7-2-2          Water
    ## ESF21MP-7-2-6          Water
    ## ESF21MP-70                MP
    ## ESF21MP-71                MP
    ## ESF21MP-72             Glass
    ## ESF21MP-75                MP
    ## ESF21MP-76                MP
    ## ESF21MP-79                MP
    ## ESF21MP-8-2-0          Water
    ## ESF21MP-8-2-10         Water
    ## ESF21MP-8-2-2          Water
    ## ESF21MP-8-2-6          Water
    ## ESF21MP-81                MP
    ## ESF21MP-82                MP
    ## ESF21MP-84                MP
    ## ESF21MP-87             Glass
    ## ESF21MP-89                MP
    ## ESF21MP-92                MP
    ## ESF21MP-95             Glass
    ## ESF21MP-96                MP
    ## ESF21MP-99                MP

``` r
mp.tww <- subset_samples(MPphyseqRA, effluent == "TWW")
sample_data(mp.tww)
```

    ##                effluent week polymer_type bead_diameter Channel sample_type
    ## ESF21MP-101         TWW   10         HDPE         3.175     5.2    Particle
    ## ESF21MP-103         TWW   10           PP         3.175     5.2    Particle
    ## ESF21MP-105         TWW   10         HDPE         3.175     8.1    Particle
    ## ESF21MP-107         TWW   10         HDPE         3.175     5.1    Particle
    ## ESF21MP-110         TWW   10         HDPE         3.175     6.2    Particle
    ## ESF21MP-113         TWW   10         LDPE         3.175     6.2    Particle
    ## ESF21MP-115         TWW   10        Glass         3.175     5.2    Particle
    ## ESF21MP-116         TWW   10           PP         3.175     6.2    Particle
    ## ESF21MP-123         TWW   10           PS         3.175     6.2    Particle
    ## ESF21MP-126         TWW   10           PP         3.175     5.1    Particle
    ## ESF21MP-127         TWW   10           PS         3.175     5.2    Particle
    ## ESF21MP-128         TWW   10        Glass         3.175     6.2    Particle
    ## ESF21MP-13          TWW    2        Glass         3.175     5.1    Particle
    ## ESF21MP-130         TWW   10         LDPE         3.175     5.2    Particle
    ## ESF21MP-131         TWW   10         LDPE         3.175     5.1    Particle
    ## ESF21MP-132         TWW   10           PS         3.175     8.1    Particle
    ## ESF21MP-133         TWW   10        Glass         3.175     5.1    Particle
    ## ESF21MP-134         TWW   10        Glass         3.175     8.1    Particle
    ## ESF21MP-135         TWW   10         LDPE         3.175     8.1    Particle
    ## ESF21MP-16          TWW    2         LDPE         3.175     8.1    Particle
    ## ESF21MP-18          TWW    2         HDPE         3.175     6.2    Particle
    ## ESF21MP-21          TWW    2        Glass         3.175     6.2    Particle
    ## ESF21MP-22          TWW    2           PP         3.175     5.2    Particle
    ## ESF21MP-23          TWW    2           PP         3.175     5.1    Particle
    ## ESF21MP-26          TWW    2           PS         3.175     6.2    Particle
    ## ESF21MP-27          TWW    2         HDPE         3.175     5.1    Particle
    ## ESF21MP-29          TWW    2         LDPE         3.175     5.1    Particle
    ## ESF21MP-35          TWW    2         LDPE         3.175     6.2    Particle
    ## ESF21MP-36          TWW    2           PP         3.175     8.1    Particle
    ## ESF21MP-37          TWW    2           PS         3.175     5.2    Particle
    ## ESF21MP-38          TWW    2         HDPE         3.175     8.1    Particle
    ## ESF21MP-41          TWW    2           PS         3.175     8.1    Particle
    ## ESF21MP-42          TWW    2         LDPE         3.175     5.2    Particle
    ## ESF21MP-43          TWW    2           PP         3.175     6.2    Particle
    ## ESF21MP-44          TWW    2        Glass         3.175     8.1    Particle
    ## ESF21MP-48          TWW    2        Glass         3.175     5.2    Particle
    ## ESF21MP-5-1-0       TWW    0        Water         0.000     5.1       Water
    ## ESF21MP-5-1-10      TWW   10        Water         0.000     5.1       Water
    ## ESF21MP-5-1-2       TWW    2        Water         0.000     5.1       Water
    ## ESF21MP-5-1-6       TWW    6        Water         0.000     5.1       Water
    ## ESF21MP-5-2-0       TWW    0        Water         0.000     5.2       Water
    ## ESF21MP-5-2-10      TWW   10        Water         0.000     5.2       Water
    ## ESF21MP-5-2-2       TWW    2        Water         0.000     5.2       Water
    ## ESF21MP-5-2-6       TWW    6        Water         0.000     5.2       Water
    ## ESF21MP-50          TWW    2         HDPE         3.175     5.2    Particle
    ## ESF21MP-53          TWW    2           PS         3.175     5.1    Particle
    ## ESF21MP-56          TWW    6           PS         3.175     5.1    Particle
    ## ESF21MP-6-2-0       TWW    0        Water         0.000     6.2       Water
    ## ESF21MP-6-2-10      TWW   10        Water         0.000     6.2       Water
    ## ESF21MP-6-2-2       TWW    2        Water         0.000     6.2       Water
    ## ESF21MP-6-2-6       TWW    6        Water         0.000     6.2       Water
    ## ESF21MP-61          TWW    6        Glass         3.175     8.1    Particle
    ## ESF21MP-62          TWW    6           PP         3.175     8.1    Particle
    ## ESF21MP-63          TWW    6           PP         3.175     5.2    Particle
    ## ESF21MP-64          TWW    6        Glass         3.175     5.2    Particle
    ## ESF21MP-65          TWW    6        Glass         3.175     6.2    Particle
    ## ESF21MP-68          TWW    6         LDPE         3.175     6.2    Particle
    ## ESF21MP-69          TWW    6        Glass         3.175     5.1    Particle
    ## ESF21MP-73          TWW    6           PP         3.175     6.2    Particle
    ## ESF21MP-74          TWW    6         HDPE         3.175     5.2    Particle
    ## ESF21MP-77          TWW    6         LDPE         3.175     8.1    Particle
    ## ESF21MP-78          TWW    6         LDPE         3.175     5.2    Particle
    ## ESF21MP-8-1-0       TWW    0        Water         0.000     8.1       Water
    ## ESF21MP-8-1-10      TWW   10        Water         0.000     8.1       Water
    ## ESF21MP-8-1-2       TWW    2        Water         0.000     8.1       Water
    ## ESF21MP-8-1-6       TWW    6        Water         0.000     8.1       Water
    ## ESF21MP-80          TWW    6           PS         3.175     6.2    Particle
    ## ESF21MP-83          TWW    6         HDPE         3.175     6.2    Particle
    ## ESF21MP-86          TWW    6         HDPE         3.175     8.1    Particle
    ## ESF21MP-88          TWW    6         HDPE         3.175     5.1    Particle
    ## ESF21MP-90          TWW    6           PS         3.175     5.2    Particle
    ## ESF21MP-91          TWW    6         LDPE         3.175     5.1    Particle
    ## ESF21MP-93          TWW    6           PP         3.175     5.1    Particle
    ## ESF21MP-94          TWW    6           PS         3.175     8.1    Particle
    ## ESF21MP-97          TWW   10           PP         3.175     8.1    Particle
    ## ESF21MP-98          TWW   10           PS         3.175     5.1    Particle
    ##                particle_type
    ## ESF21MP-101               MP
    ## ESF21MP-103               MP
    ## ESF21MP-105               MP
    ## ESF21MP-107               MP
    ## ESF21MP-110               MP
    ## ESF21MP-113               MP
    ## ESF21MP-115            Glass
    ## ESF21MP-116               MP
    ## ESF21MP-123               MP
    ## ESF21MP-126               MP
    ## ESF21MP-127               MP
    ## ESF21MP-128            Glass
    ## ESF21MP-13             Glass
    ## ESF21MP-130               MP
    ## ESF21MP-131               MP
    ## ESF21MP-132               MP
    ## ESF21MP-133            Glass
    ## ESF21MP-134            Glass
    ## ESF21MP-135               MP
    ## ESF21MP-16                MP
    ## ESF21MP-18                MP
    ## ESF21MP-21             Glass
    ## ESF21MP-22                MP
    ## ESF21MP-23                MP
    ## ESF21MP-26                MP
    ## ESF21MP-27                MP
    ## ESF21MP-29                MP
    ## ESF21MP-35                MP
    ## ESF21MP-36                MP
    ## ESF21MP-37                MP
    ## ESF21MP-38                MP
    ## ESF21MP-41                MP
    ## ESF21MP-42                MP
    ## ESF21MP-43                MP
    ## ESF21MP-44             Glass
    ## ESF21MP-48             Glass
    ## ESF21MP-5-1-0          Water
    ## ESF21MP-5-1-10         Water
    ## ESF21MP-5-1-2          Water
    ## ESF21MP-5-1-6          Water
    ## ESF21MP-5-2-0          Water
    ## ESF21MP-5-2-10         Water
    ## ESF21MP-5-2-2          Water
    ## ESF21MP-5-2-6          Water
    ## ESF21MP-50                MP
    ## ESF21MP-53                MP
    ## ESF21MP-56                MP
    ## ESF21MP-6-2-0          Water
    ## ESF21MP-6-2-10         Water
    ## ESF21MP-6-2-2          Water
    ## ESF21MP-6-2-6          Water
    ## ESF21MP-61             Glass
    ## ESF21MP-62                MP
    ## ESF21MP-63                MP
    ## ESF21MP-64             Glass
    ## ESF21MP-65             Glass
    ## ESF21MP-68             Glass
    ## ESF21MP-69             Glass
    ## ESF21MP-73                MP
    ## ESF21MP-74                MP
    ## ESF21MP-77                MP
    ## ESF21MP-78                MP
    ## ESF21MP-8-1-0          Water
    ## ESF21MP-8-1-10         Water
    ## ESF21MP-8-1-2          Water
    ## ESF21MP-8-1-6          Water
    ## ESF21MP-80                MP
    ## ESF21MP-83                MP
    ## ESF21MP-86                MP
    ## ESF21MP-88                MP
    ## ESF21MP-90                MP
    ## ESF21MP-91                MP
    ## ESF21MP-93                MP
    ## ESF21MP-94                MP
    ## ESF21MP-97                MP
    ## ESF21MP-98                MP

## Alpha Diversity

Note that these data are not rarefied. We’ll use the Shannon diversity
index (based on richness AND evenness; examines how many different taxa
are present and how evenly they’re distributed within a sample) to
analyze alpha diversity between variable types. This means it considers
both the number of species and the inequality between species
abundances.

### Calculate Shannon Diversity per Sample

Using the estimate_richness function in the phyloseq package to
calculate Shannon diversity. The estimate_richness function can also
take measures “Chao1” “ACE” “Simpson” and “Fisher”.

``` r
#Calculating richness - Shannon diversity in a new dataframe
richness <- data.frame(estimate_richness(MPphyseqRA, measures = c("Shannon")))
```

    ## Warning in estimate_richness(MPphyseqRA, measures = c("Shannon")): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
richness <- setNames(cbind(rownames(richness), richness, row.names = NULL), 
                     c("sample-id", "Shannon"))

#Add the sample metadata to the dataframe
s <- data.frame(sample_data(MPphyseqRA))
s <- setNames(cbind(rownames(s), s, row.names = NULL), 
              c("sample-id", "effluent", "week", "polymer_type", 
                "bead_diameter", "Channel", "sample_type", "particle_type"))

alphadiv <- merge(s, richness, by = "sample-id")

#Order factors
alphadiv$polymer_type <- factor(alphadiv$polymer_type, 
                        levels = c("Glass", "HDPE", "LDPE", "PP", "PS", "Water"),
                        labels = c("Glass", "HDPE", "LDPE", "PP", "PS", "Water")) 
#Shows the calculated indices
knitr::kable(head(alphadiv)) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div
style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
sample-id
</th>
<th style="text-align:left;">
effluent
</th>
<th style="text-align:left;">
week
</th>
<th style="text-align:left;">
polymer_type
</th>
<th style="text-align:right;">
bead_diameter
</th>
<th style="text-align:left;">
Channel
</th>
<th style="text-align:left;">
sample_type
</th>
<th style="text-align:left;">
particle_type
</th>
<th style="text-align:right;">
Shannon
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ESF21MP-100
</td>
<td style="text-align:left;">
CON
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
HDPE
</td>
<td style="text-align:right;">
3.175
</td>
<td style="text-align:left;">
7.1
</td>
<td style="text-align:left;">
Particle
</td>
<td style="text-align:left;">
MP
</td>
<td style="text-align:right;">
5.500475
</td>
</tr>
<tr>
<td style="text-align:left;">
ESF21MP-101
</td>
<td style="text-align:left;">
TWW
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
HDPE
</td>
<td style="text-align:right;">
3.175
</td>
<td style="text-align:left;">
5.2
</td>
<td style="text-align:left;">
Particle
</td>
<td style="text-align:left;">
MP
</td>
<td style="text-align:right;">
5.285657
</td>
</tr>
<tr>
<td style="text-align:left;">
ESF21MP-102
</td>
<td style="text-align:left;">
CON
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
LDPE
</td>
<td style="text-align:right;">
3.175
</td>
<td style="text-align:left;">
7.2
</td>
<td style="text-align:left;">
Particle
</td>
<td style="text-align:left;">
MP
</td>
<td style="text-align:right;">
5.045258
</td>
</tr>
<tr>
<td style="text-align:left;">
ESF21MP-103
</td>
<td style="text-align:left;">
TWW
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
PP
</td>
<td style="text-align:right;">
3.175
</td>
<td style="text-align:left;">
5.2
</td>
<td style="text-align:left;">
Particle
</td>
<td style="text-align:left;">
MP
</td>
<td style="text-align:right;">
5.439320
</td>
</tr>
<tr>
<td style="text-align:left;">
ESF21MP-104
</td>
<td style="text-align:left;">
CON
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
PP
</td>
<td style="text-align:right;">
3.175
</td>
<td style="text-align:left;">
7.2
</td>
<td style="text-align:left;">
Particle
</td>
<td style="text-align:left;">
MP
</td>
<td style="text-align:right;">
5.174548
</td>
</tr>
<tr>
<td style="text-align:left;">
ESF21MP-105
</td>
<td style="text-align:left;">
TWW
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
HDPE
</td>
<td style="text-align:right;">
3.175
</td>
<td style="text-align:left;">
8.1
</td>
<td style="text-align:left;">
Particle
</td>
<td style="text-align:left;">
MP
</td>
<td style="text-align:right;">
3.987267
</td>
</tr>
</tbody>
</table>

</div>

## Session Info

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.1 (2022-06-23 ucrt)
    ##  os       Windows 10 x64 (build 22000)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       America/New_York
    ##  date     2023-04-28
    ##  pandoc   2.18 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package              * version    date (UTC) lib source
    ##    ade4                   1.7-22     2023-02-06 [1] CRAN (R 4.2.3)
    ##    annotate               1.76.0     2022-11-03 [1] Bioconductor
    ##    AnnotationDbi          1.60.2     2023-03-10 [1] Bioconductor
    ##    ape                  * 5.7-1      2023-03-13 [1] CRAN (R 4.2.3)
    ##    backports              1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
    ##    base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
    ##    Biobase              * 2.58.0     2022-11-01 [1] Bioconductor
    ##    BiocGenerics         * 0.44.0     2022-11-01 [1] Bioconductor
    ##    BiocParallel           1.32.6     2023-03-17 [1] Bioconductor
    ##    biomformat             1.26.0     2022-11-01 [1] Bioconductor
    ##    Biostrings             2.66.0     2022-11-01 [1] Bioconductor
    ##    bit                    4.0.5      2022-11-15 [1] CRAN (R 4.2.3)
    ##    bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.2.3)
    ##    bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
    ##    blob                   1.2.4      2023-03-17 [1] CRAN (R 4.2.3)
    ##    cachem                 1.0.7      2023-02-24 [1] CRAN (R 4.2.3)
    ##    callr                  3.7.3      2022-11-02 [1] CRAN (R 4.2.3)
    ##    checkmate              2.1.0      2022-04-21 [1] CRAN (R 4.2.3)
    ##    cli                    3.6.1      2023-03-23 [1] CRAN (R 4.2.3)
    ##    cluster                2.1.4      2022-08-22 [1] CRAN (R 4.2.3)
    ##    codetools              0.2-19     2023-02-01 [1] CRAN (R 4.2.2)
    ##    colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.2.3)
    ##    crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.2.3)
    ##    data.table             1.14.8     2023-02-17 [1] CRAN (R 4.2.3)
    ##    DBI                    1.1.3      2022-06-18 [1] CRAN (R 4.2.3)
    ##    DelayedArray           0.24.0     2022-11-01 [1] Bioconductor
    ##    DESeq2               * 1.38.3     2023-01-19 [1] Bioconductor
    ##    devtools               2.4.5      2022-10-11 [1] CRAN (R 4.2.3)
    ##    digest                 0.6.31     2022-12-11 [1] CRAN (R 4.2.3)
    ##    dplyr                * 1.1.1      2023-03-22 [1] CRAN (R 4.2.3)
    ##    DT                     0.27       2023-01-17 [1] CRAN (R 4.2.3)
    ##    ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.2.3)
    ##    evaluate               0.20       2023-01-17 [1] CRAN (R 4.2.3)
    ##    fansi                  1.0.4      2023-01-22 [1] CRAN (R 4.2.3)
    ##    farver                 2.1.1      2022-07-06 [1] CRAN (R 4.2.3)
    ##    fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.2.3)
    ##    foreach                1.5.2      2022-02-02 [1] CRAN (R 4.2.3)
    ##    foreign                0.8-84     2022-12-06 [1] CRAN (R 4.2.2)
    ##    Formula                1.2-5      2023-02-24 [1] CRAN (R 4.2.2)
    ##    fs                     1.6.1      2023-02-06 [1] CRAN (R 4.2.3)
    ##    geneplotter            1.76.0     2022-11-01 [1] Bioconductor
    ##    generics               0.1.3      2022-07-05 [1] CRAN (R 4.2.3)
    ##    GenomeInfoDb         * 1.34.9     2023-02-02 [1] Bioconductor
    ##    GenomeInfoDbData       1.2.9      2023-04-07 [1] Bioconductor
    ##    GenomicRanges        * 1.50.2     2022-12-27 [1] Bioconductor
    ##    ggplot2              * 3.4.2      2023-04-03 [1] CRAN (R 4.2.3)
    ##    glue                   1.6.2      2022-02-24 [1] CRAN (R 4.2.3)
    ##    gridExtra              2.3        2017-09-09 [1] CRAN (R 4.2.3)
    ##    gtable                 0.3.3      2023-03-21 [1] CRAN (R 4.2.3)
    ##    highr                  0.10       2022-12-22 [1] CRAN (R 4.2.3)
    ##    Hmisc                  5.0-1      2023-03-08 [1] CRAN (R 4.2.3)
    ##    htmlTable              2.4.1      2022-07-07 [1] CRAN (R 4.2.3)
    ##    htmltools              0.5.5      2023-03-23 [1] CRAN (R 4.2.3)
    ##    htmlwidgets            1.6.2      2023-03-17 [1] CRAN (R 4.2.3)
    ##    httpuv                 1.6.9      2023-02-14 [1] CRAN (R 4.2.3)
    ##    httr                   1.4.5      2023-02-24 [1] CRAN (R 4.2.3)
    ##    igraph                 1.4.1      2023-02-24 [1] CRAN (R 4.2.3)
    ##    IRanges              * 2.32.0     2022-11-01 [1] Bioconductor
    ##    iterators              1.0.14     2022-02-05 [1] CRAN (R 4.2.3)
    ##    jsonlite               1.8.4      2022-12-06 [1] CRAN (R 4.2.3)
    ##    kableExtra             1.3.4      2021-02-20 [1] CRAN (R 4.2.3)
    ##    KEGGREST               1.38.0     2022-11-01 [1] Bioconductor
    ##    knitr                  1.42       2023-01-25 [1] CRAN (R 4.2.3)
    ##    labeling               0.4.2      2020-10-20 [1] CRAN (R 4.2.0)
    ##    later                  1.3.0      2021-08-18 [1] CRAN (R 4.2.3)
    ##    lattice              * 0.21-8     2023-04-05 [1] CRAN (R 4.2.3)
    ##    lifecycle              1.0.3      2022-10-07 [1] CRAN (R 4.2.3)
    ##    locfit                 1.5-9.7    2023-01-02 [1] CRAN (R 4.2.3)
    ##    magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.2.3)
    ##    MASS                   7.3-58.3   2023-03-07 [1] CRAN (R 4.2.3)
    ##    Matrix                 1.5-4      2023-04-04 [1] CRAN (R 4.2.3)
    ##    MatrixGenerics       * 1.10.0     2022-11-01 [1] Bioconductor
    ##    matrixStats          * 0.63.0     2022-11-18 [1] CRAN (R 4.2.3)
    ##    memoise                2.0.1      2021-11-26 [1] CRAN (R 4.2.3)
    ##    mgcv                   1.8-42     2023-03-02 [1] CRAN (R 4.2.3)
    ##    mime                   0.12       2021-09-28 [1] CRAN (R 4.2.0)
    ##    miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.2.3)
    ##    multtest               2.54.0     2022-11-01 [1] Bioconductor
    ##    munsell                0.5.0      2018-06-12 [1] CRAN (R 4.2.3)
    ##    NADA                   1.6-1.1    2020-03-22 [1] CRAN (R 4.2.3)
    ##    nlme                 * 3.1-162    2023-01-31 [1] CRAN (R 4.2.3)
    ##    nnet                   7.3-18     2022-09-28 [1] CRAN (R 4.2.3)
    ##    patchwork            * 1.1.2.9000 2023-04-24 [1] Github (thomasp85/patchwork@c14c960)
    ##    permute              * 0.9-7      2022-01-27 [1] CRAN (R 4.2.3)
    ##    phyloseq             * 1.42.0     2022-11-01 [1] Bioconductor
    ##    picante              * 1.8.2      2020-06-10 [1] CRAN (R 4.2.3)
    ##    pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.2.3)
    ##    pkgbuild               1.4.0      2022-11-27 [1] CRAN (R 4.2.3)
    ##    pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.2.3)
    ##    pkgload                1.3.2      2022-11-16 [1] CRAN (R 4.2.3)
    ##    plyr                 * 1.8.8      2022-11-11 [1] CRAN (R 4.2.3)
    ##    png                    0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
    ##    prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.2.3)
    ##    processx               3.8.1      2023-04-18 [1] CRAN (R 4.2.1)
    ##    profvis                0.3.7      2020-11-02 [1] CRAN (R 4.2.3)
    ##    promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.2.3)
    ##    ps                     1.7.4      2023-04-02 [1] CRAN (R 4.2.3)
    ##    purrr                  1.0.1      2023-01-10 [1] CRAN (R 4.2.3)
    ##    qiime2R              * 0.99.6     2023-04-07 [1] Github (jbisanz/qiime2R@2a3cee1)
    ##    R6                     2.5.1      2021-08-19 [1] CRAN (R 4.2.3)
    ##    RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
    ##    Rcpp                   1.0.10     2023-01-22 [1] CRAN (R 4.2.3)
    ##    RCurl                  1.98-1.12  2023-03-27 [1] CRAN (R 4.2.3)
    ##    remotes                2.4.2      2021-11-30 [1] CRAN (R 4.2.3)
    ##    reshape2             * 1.4.4      2020-04-09 [1] CRAN (R 4.2.3)
    ##    rhdf5                  2.42.0     2022-11-01 [1] Bioconductor
    ##  D rhdf5filters           1.10.1     2023-03-24 [1] Bioconductor
    ##    Rhdf5lib               1.20.0     2022-11-01 [1] Bioconductor
    ##    rlang                  1.1.0      2023-03-14 [1] CRAN (R 4.2.3)
    ##    rmarkdown              2.21       2023-03-26 [1] CRAN (R 4.2.3)
    ##    rpart                  4.1.19     2022-10-21 [1] CRAN (R 4.2.3)
    ##    RSQLite                2.3.1      2023-04-03 [1] CRAN (R 4.2.3)
    ##    rstudioapi             0.14       2022-08-22 [1] CRAN (R 4.2.3)
    ##    rvest                  1.0.3      2022-08-19 [1] CRAN (R 4.2.3)
    ##    S4Vectors            * 0.36.2     2023-02-26 [1] Bioconductor
    ##    scales               * 1.2.1      2022-08-20 [1] CRAN (R 4.2.3)
    ##    sessioninfo            1.2.2      2021-12-06 [1] CRAN (R 4.2.3)
    ##    shiny                  1.7.4      2022-12-15 [1] CRAN (R 4.2.3)
    ##    stringi                1.7.12     2023-01-11 [1] CRAN (R 4.2.2)
    ##    stringr                1.5.0      2022-12-02 [1] CRAN (R 4.2.3)
    ##    SummarizedExperiment * 1.28.0     2022-11-01 [1] Bioconductor
    ##    survival               3.5-5      2023-03-12 [1] CRAN (R 4.2.3)
    ##    svglite                2.1.1      2023-01-10 [1] CRAN (R 4.2.3)
    ##    systemfonts            1.0.4      2022-02-11 [1] CRAN (R 4.2.3)
    ##    tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.2.3)
    ##    tidyr                * 1.3.0      2023-01-24 [1] CRAN (R 4.2.3)
    ##    tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.2.3)
    ##    truncnorm              1.0-9      2023-03-20 [1] CRAN (R 4.2.3)
    ##    urlchecker             1.0.1      2021-11-30 [1] CRAN (R 4.2.3)
    ##    usethis                2.1.6      2022-05-25 [1] CRAN (R 4.2.3)
    ##    utf8                   1.2.3      2023-01-31 [1] CRAN (R 4.2.3)
    ##    vctrs                  0.6.1      2023-03-22 [1] CRAN (R 4.2.3)
    ##    vegan                * 2.6-4      2022-10-11 [1] CRAN (R 4.2.3)
    ##    viridis              * 0.6.2      2021-10-13 [1] CRAN (R 4.2.3)
    ##    viridisLite          * 0.4.1      2022-08-22 [1] CRAN (R 4.2.3)
    ##    webshot                0.5.4      2022-09-26 [1] CRAN (R 4.2.3)
    ##    withr                  2.5.0      2022-03-03 [1] CRAN (R 4.2.3)
    ##    xfun                   0.38       2023-03-24 [1] CRAN (R 4.2.3)
    ##    XML                    3.99-0.14  2023-03-19 [1] CRAN (R 4.2.3)
    ##    xml2                   1.3.3      2021-11-30 [1] CRAN (R 4.2.3)
    ##    xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.2.3)
    ##    XVector                0.38.0     2022-11-01 [1] Bioconductor
    ##    yaml                   2.3.7      2023-01-23 [1] CRAN (R 4.2.3)
    ##    zCompositions          1.4.0-1    2022-03-26 [1] CRAN (R 4.2.3)
    ##    zlibbioc               1.44.0     2022-11-01 [1] Bioconductor
    ## 
    ##  [1] C:/Program Files/R/R-4.2.1/library
    ## 
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
