#' Generate bootstrap index data
#'
#' @param celltype colnames of single cell data.Colnames should be labeled as cell type or cluster.
#' @param n  number of bootstrap
#' @import dplyr
#'
#' @return
#'
#'
#' @examples
generate_boots <- function(celltype, n) {
  dt <- data.frame(cluster = celltype, id = 1:length(celltype))
  index <- do.call(cbind, sapply(1:n, function(x) {
    splits <- dt %>%
      group_by(cluster) %>%
      sample_n(dplyr::n(), replace = TRUE) %>%
      ungroup() %>%
      dplyr::select("id")
  }))
  return(index)
}




#
#' Calculate mean expression for one bootstrap
#'
#' @param i index
#' @param raw single cell data
#' @param samples generated bootstrap index data
#'
#' @return
#'
#'
#' @examples
get_ave_exp <- function(i, raw, samples) {
  sample <- as.data.frame(t(raw[, samples[, i, drop = F]]))
  sname <- colnames(raw)[samples[, i, drop = F]]
  k <- split(sample, sname)
  averg <- do.call(cbind, lapply(k, colMeans))
  return(averg)
}



#' Calculate bootstrap mean expression for single cell data
#'
#' @param raw data frame of single cell data with celltype or cluster as colnames
#' @param n_bootstrap number of bootstrap
#' @param seed random seed
#'
#' @return mean expression data
#' @export
#'
#' @examples
calculate_avg_exp <- function(raw, n_bootstrap, seed) {
  set.seed(seed)
  samples <- generate_boots(celltype = colnames(raw), n = n_bootstrap)
  exp <- lapply(1:n_bootstrap, get_ave_exp, raw, samples)
  exp <- do.call(cbind, exp)
  return(exp)
}


