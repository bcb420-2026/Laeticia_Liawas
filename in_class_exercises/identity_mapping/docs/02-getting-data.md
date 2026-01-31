# Getting data from GEO

We download GEO supplementary files to obtain realistic identifiers.


``` r
library(GEOquery)
```

```
## Loading required package: Biobase
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: generics
```

```
## 
## Attaching package: 'generics'
```

```
## The following objects are masked from 'package:base':
## 
##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
##     setequal, union
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
##     unsplit, which.max, which.min
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

``` r
source("./fetch_geo_supp.R")

fetch_geo_supp(gse = "GSE233947")
```

```
## Using locally cached version of supplementary file(s) GSE233947 found here:
## data/GSE233947/GSE233947_FeatureCounts_V31genes_RawCounts_ENSG.tsv.gz
```

```
## Using locally cached version of supplementary file(s) GSE233947 found here:
## data/GSE233947/GSE233947_modulize_3CTG_20CTG_junctions.tsv.gz
```

```
## Using locally cached version of supplementary file(s) GSE233947 found here:
## data/GSE233947/GSE233947_modulize_NT_20CTG_junctions.tsv.gz
```

```
## Using locally cached version of supplementary file(s) GSE233947 found here:
## data/GSE233947/GSE233947_modulize_NT_3CTG_junctions.tsv.gz
```

``` r
fetch_geo_supp(gse = "GSE119732")
```

```
## Using locally cached version of supplementary file(s) GSE119732 found here:
## data/GSE119732/GSE119732_count_table_RNA_seq.txt.gz
```

