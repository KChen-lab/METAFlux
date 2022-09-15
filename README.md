# METAFlux: Characterizing metabolism from bulk and single-cell RNA-seq data 
<p align="center">
  <img width="200"  src="https://github.com/KChen-lab/METAFlux/blob/main/METAFlux%20logo.jpeg">
</p>

## Capabilities
* Use genome scale metabolic model and FBA (Flux balance analysis) with gene expression to derive metabolic fluxes for bulk and single cell RNA-seq datasets
* Discover various modes of metabolic cooperation and competition in single-cell RNA-seq data
* Has the potential to improve the understanding of aberrant metabolism and serve as a preliminary source to investigate specific metabolic targets.

Check out the full paper for more detailed methods and applications: 
**Characterizing metabolism from bulk and single-cell RNA-seq data using METAFlux**

## Installation 
METAFlux R package can be easily installed from Github using devtools:

`devtools::install_github('Rhyf/METAFlux')`

### Installation of Other Dependencies
* Install the osqp package for optimization using `install.packages('osqp')`, If you encounter any issue during METAFlux installation.
* Install the dplyr package using `install.packages('dplyr')`.



## Tutorial
Please check out this link for a full tutorial on using METAFlux:
[METAFlux tutorial](https://rhyf.github.io/METAFlux/)


