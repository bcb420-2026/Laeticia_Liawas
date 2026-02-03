# Normalization methods

Different assays require different methods.

## Library-size scaling

The simplest approach divides each sample by its total counts (or reads) to account for sequencing depth.

### Counts-per-million (CPM)

There are other versions of this type of normalization (RPKM, FPKM ,TPM,...) - all of them try to normalize for the library size (the number of reads in a given experiment) but some of them incorporate other factors as well.  

CPM - normalizes for just library size.  
TPM - normalzie for gene length, then library size.
FPKM/RPKM - similiar but one uses reads and the other fragments.  normalizes for gene length and library size at the same time.  TPM is preferred.   

## Distribution-based normalization

Methods like **quantile normalization** force sample distributions to match and are common in microarrays.

## RNA-seq composition-aware methods

Two widely used methods:

- **TMM** (edgeR): computes scaling factors to adjust effective library sizes.
- **RLE / size factors** (DESeq2): median-of-ratios approach.

### edgeR: TMM

TMM looks at most genes, ignores the weird ones,
figures out how different two samples really are,
and uses that to scale the samples so they’re comparable.

```r
library(edgeR)

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
```

### DESeq2: size factors

RLE finds what “most genes” are doing compared to a typical sample.
It adjusts each sample so that the median gene looks the same in all samples.

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)
```

## Single-cell note

In scRNA-seq, normalization often starts with per-cell library-size scaling followed by `log1p` and sometimes variance stabilization. The concepts (make samples comparable; avoid composition artifacts) are the same.

