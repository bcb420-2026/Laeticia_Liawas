# State base image
FROM risserlin/bcb420-base-image:winter2026

# Install additional packages from Bioconductor
RUN R -e 'BiocManager::install(c("enrichplot","DESeq2","GEOquery","ComplexHeatmap"))'

# Install additional packages from CRAN 
RUN install2.r -d TRUE -r "https://cran.rstudio.com" pheatmap ggplot2
