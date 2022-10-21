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
#' @param myseurat single cell Seurat object
#' @param samples generated bootstrap index data
#' @param myident Seurat idents.This will be a character string indicating the grouping of the seurat object
#'
#' @return
#' @import Seurat
#'
#' @examples

get_ave_exp <- function(i, myseurat, samples,myident) {
  meta.data=myseurat@meta.data[samples[,i],]
  sample <-myseurat@assays$RNA@counts[,samples[,i]]
  SeuratObject<-suppressWarnings(
    CreateSeuratObject(count=sample,meta.data = meta.data))
  SeuratObject<-NormalizeData(SeuratObject,verbose = FALSE)
  ave<-AverageExpression(SeuratObject,group.by = myident,return.seurat = T)[["RNA"]]@data
  return(ave)
}


#' Calculate bootstrap mean expression for single cell data
#'
#' @param n_bootstrap number of bootstrap
#' @param seed random seed
#' @param myseurat Seurat object
#' @param myident Seurat idents for grouping.This will be a character string indicating the grouping of the seurat object
#'
#' @return mean expression data
#' @import Seurat
#' @export
#'
#' @examples
calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
  set.seed(seed)
  samples=generate_boots(myseurat@meta.data[,myident],n_bootstrap)
  exp <- lapply(1:n_bootstrap,get_ave_exp,myseurat,samples,myident)
  exp <- do.call(cbind, exp)
  return(exp)
}


