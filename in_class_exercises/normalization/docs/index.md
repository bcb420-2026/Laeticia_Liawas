---
title: "Normalization (BCB420 Lecture)"
author: "Ruth Isserlin"
description: "A short bookdown module introducing normalization for expression data (bulk and single-cell)."
site: bookdown::bookdown_site
documentclass: book
bibliography: [references.bib]
link-citations: true
---

# Welcome {-}

This bookdown module turns the **Lecture – Normalization** slides into a hands-on, reproducible mini-book.

## Learning goals

By the end, learners should be able to:

- Visualize sample distributions with boxplots and density plots.
- Apply appropriate normalization methods for:
  - **bulk RNA-seq** (e.g., TMM, DESeq2 size factors)
  - **single-cell RNA-seq** (library-size scaling + log1p; conceptually)
- Evaluate normalization with diagnostic plots.

## Example datasets

We re-use the same GEO series as practical examples:

- **GSE119732** (bulk RNA-seq)


