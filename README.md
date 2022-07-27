# METAFlux
## Capabilities
* Use genome scale metabolic model and FBA (Flux balance analysis) with gene expression to derive metabolic fluxes for 13,082 reactions for bulk and RNA-seq datasets
* Discover various modes of metabolic cooperation and competition in single-cell data
* Has the potential to improve the understanding of aberrant metabolism and serve as a preliminary source to investigate specific metabolic targets.

Check out the full paper for more detailed methods and applications: 
**Characterizing metabolism from bulk and single-cell RNA-seq data using METAFlux**

## Installation 
METAFlux R package can be easily installed from Github using devtools:

`devtools::install_github('Rhyf/METAFlux',build_vignettes = T)`

### Installation of Other Dependencies
* Install the osqp package for optimization using `install.packages('osqp')`.
* Install the dplyr package using `install.packages('dplyr')`.
* Install the utils package using `install.packages('utils')`.


## Tutorial
Please check out this link for a full tutorial on using METAFlux:
[METAFlux tutorial](https://rhyf.github.io/METAFlux/)


