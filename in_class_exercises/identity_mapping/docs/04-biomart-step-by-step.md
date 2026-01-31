---
title: "Biomart step by step"
params:
  ensembl_version: 115
---

# biomaRt step-by-step

This follows the following workflow: **choose a mart → choose a dataset → discover filters/attributes → query with getBM() → merge back into your matrix**.


``` r
library(biomaRt)
```


## 1. List marts


``` r
listMarts()
```

```
##                biomart                version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 115
## 2   ENSEMBL_MART_MOUSE      Mouse strains 115
## 3     ENSEMBL_MART_SNP  Ensembl Variation 115
## 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 115
```

## 2. (Optional) pin an Ensembl archive version


``` r
listEnsemblArchives()[1:10,]
```

```
##              name     date                                 url version
## 1  Ensembl GRCh37 Feb 2014          https://grch37.ensembl.org  GRCh37
## 2     Ensembl 115 Sep 2025 https://sep2025.archive.ensembl.org     115
## 3     Ensembl 114 May 2025 https://may2025.archive.ensembl.org     114
## 4     Ensembl 113 Oct 2024 https://oct2024.archive.ensembl.org     113
## 5     Ensembl 112 May 2024 https://may2024.archive.ensembl.org     112
## 6     Ensembl 111 Jan 2024 https://jan2024.archive.ensembl.org     111
## 7     Ensembl 110 Jul 2023 https://jul2023.archive.ensembl.org     110
## 8     Ensembl 109 Feb 2023 https://feb2023.archive.ensembl.org     109
## 9     Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
## 10    Ensembl 107 Jul 2022 https://jul2022.archive.ensembl.org     107
##    current_release
## 1                 
## 2                *
## 3                 
## 4                 
## 5                 
## 6                 
## 7                 
## 8                 
## 9                 
## 10
```

``` r
# Example: pin to Ensembl 114
ensembl <- useEnsembl(biomart = "ensembl",version = params$ensembl_version)
```

If you do not pin versions:


``` r
ensembl <- useMart("ensembl")
```

## 3. List datasets and select human


``` r
datasets <- listDatasets(ensembl)

# filter for human
human <- datasets[grep(datasets$dataset, pattern = "hsapiens"), ]
human
```

```
##                  dataset              description    version
## 80 hsapiens_gene_ensembl Human genes (GRCh38.p14) GRCh38.p14
```

``` r
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
```

## 4. Identify the correct filter

Filters are your **input ID type**.


``` r
all_filters <- listFilters(ensembl)
dim(all_filters)
```

```
## [1] 443   2
```

``` r
# search for Ensembl gene filters
all_filters[grep(all_filters$name, pattern = "ensembl_gene"), ]
```

```
##                       name
## 54         ensembl_gene_id
## 55 ensembl_gene_id_version
##                                                 description
## 54                 Gene stable ID(s) [e.g. ENSG00000000003]
## 55 Gene stable ID(s) with version [e.g. ENSG00000000003.17]
```

For Ensembl gene IDs, use filter: `ensembl_gene_id`.

## 5. Identify the correct attributes

Attributes are the **output columns** you want.


``` r
all_attr <- listAttributes(ensembl)
all_attr[1:10,]
```

```
##                             name                  description         page
## 1                ensembl_gene_id               Gene stable ID feature_page
## 2        ensembl_gene_id_version       Gene stable ID version feature_page
## 3          ensembl_transcript_id         Transcript stable ID feature_page
## 4  ensembl_transcript_id_version Transcript stable ID version feature_page
## 5             ensembl_peptide_id            Protein stable ID feature_page
## 6     ensembl_peptide_id_version    Protein stable ID version feature_page
## 7                ensembl_exon_id               Exon stable ID feature_page
## 8                    description             Gene description feature_page
## 9                chromosome_name     Chromosome/scaffold name feature_page
## 10                start_position              Gene start (bp) feature_page
```

``` r
# search for HGNC
searchAttributes(ensembl, "hgnc")
```

```
##               name        description         page
## 63     hgnc_symbol        HGNC symbol feature_page
## 64         hgnc_id            HGNC ID feature_page
## 95 hgnc_trans_name Transcript name ID feature_page
```

Common attributes:

- `ensembl_gene_id`
- `hgnc_symbol`

## 6. Strip version suffixes (if present)


``` r
strip_ensembl_version <- function(x) sub("\\..*$", "", x)
ids <- c("ENSG00000141510.17", "ENSG00000157764.2")
ids_clean <- strip_ensembl_version(ids)
ids_clean
```

```
## [1] "ENSG00000141510" "ENSG00000157764"
```

## 7. Run the query with `getBM()`


``` r
map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ids_clean,
  mart = ensembl
)
map
```

```
##   ensembl_gene_id hgnc_symbol
## 1 ENSG00000141510        TP53
## 2 ENSG00000157764        BRAF
```

## 8. Cache the mapping (recommended)


``` r
cache_file <- "id_conversion.rds"

if (file.exists(cache_file)) {
  map <- readRDS(cache_file)
} else {
  map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ids_clean,
    mart = ensembl
  )
  saveRDS(map, cache_file)
}
```

## 9. Merge back into your data


``` r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:biomaRt':
## 
##     select
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
df <- tibble(ensembl_gene_id = ids_clean, value = c(1, 2))
df_annot <- left_join(df, map, by = "ensembl_gene_id")
df_annot
```

```
## # A tibble: 2 × 3
##   ensembl_gene_id value hgnc_symbol
##   <chr>           <dbl> <chr>      
## 1 ENSG00000141510     1 TP53       
## 2 ENSG00000157764     2 BRAF
```


