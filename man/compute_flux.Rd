% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{compute_flux}
\alias{compute_flux}
\title{Final optimization step for flux calculation}
\usage{
compute_flux(mras, medium)
}
\arguments{
\item{mras}{metabolic reaction activity scores}

\item{medium}{input medium file which indicates the nutrients available in the medium.
We provide 2 general mediums if you have no prior knowledge about your medium: cell line medium and human blood medium if prior knowledge of medium is not available.
Please see tutorial for more details.}
}
\value{
Calculated fluxes
}
\description{
Final optimization step for flux calculation
}
