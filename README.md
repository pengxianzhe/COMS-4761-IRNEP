# COMS-4761-IRNEP

All codes are running in R-3.4.1, co-expression network inference is memory intensive, so it's better to use machines with large memories.

### Code
Data_prep.R: scripts for single-cell data pre-processing and quality control, filtering.

Build_regulons.R: scripts for building the co-expression network and finding the enriched motifs and pruning the regulons

Viper_tsne.R: scripts for VIPER inference and t-SNE analysis, inputs are normalized expression data from Data_prep.R and regulons generated from Build_regulons.R

Seurat_and_tsne.R: scripts for algorithm evaluations and visualizations

### Data
Brain data 1: GSE76381_ESMoleculeCounts.cef.txt.gz from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76381

Brain data 2: GSE76381_EmbryoMoleculeCounts.cef.txt.gz from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76381

Pancreas dataset: a sample subset from cluster 0: cluster0.matrix.rda

Combined CHIP-seq evidence for TF-targets relationship: gene_chea_matrix.txt

Regulons are directly inferred from corresponding expression dataset, if we are using multiple regulons, we combine regulons from two brain datasets, plus the GBM network from https://bioconductor.org/packages/release/data/experiment/html/aracne.networks.html

### Figures
All figures are included in Figures.pptx
