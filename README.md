# METAFlux: Characterizing metabolism from bulk and single-cell RNA-seq data 
<p align="center">
  <img width="200"  src="https://github.com/KChen-lab/METAFlux/blob/main/METAFlux%20logo.jpeg">
</p>

## Capabilities
* Use genome scale metabolic model and FBA (Flux balance analysis) to derive 13,082 metabolic fluxes for bulk and single cell RNA-seq datasets
* Discover different modes of metabolic cooperation and competition in single-cell RNA-seq data
* Enable massive public RNA-seq datasets mining in cancer metabolism


Check out the full paper for more detailed methods and applications: 
**Characterizing metabolism from bulk and single-cell RNA-seq data using METAFlux**

## Installation 
METAFlux R package can be easily installed from Github using devtools:

`devtools::install_github('KChen-lab/METAFlux')`

### Installation of Other Dependencies
* Install the osqp package for optimization using `install.packages('osqp')`, If you encounter any issue during METAFlux installation.
* Install the dplyr package using `install.packages('dplyr')`.
* For single cell data analysis, we provide pipeline to work with Seurat. Please install Seurat package
 by `install.packages('Seurat')`.


## Tutorial
Please check out this link for a full tutorial on using METAFlux:

- [Full tutorial using METFAFlux for bulk and single-cell RNA seq data analysis with detailed explanation](https://htmlpreview.github.io/?https://github.com/KChen-lab/METAFlux/blob/main/Tutorials/pipeline.html)
