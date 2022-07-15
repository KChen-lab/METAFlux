#' Calculate reaction score for isoenzymes (or relationship only)
#'
#' @param x index
#' @param data input gene expression data
#' @param list isoenzymes reaction list
#' @param gene_num file of number of times one gene has involved in all pathways
#'
#' @return Calculated isoenzyme list scores
#'
#'
#' @examples
calculate_iso_score <- function(x, data, list, gene_num) {
  data_feature <- list[[x]][list[[x]] %in% rownames(data)]
  if (length(data_feature) > 0) {
    vec <- gene_num[data_feature, ]$V1
    expr <- as.matrix(data[data_feature, , drop = F])
    norm <- colSums(sweep(expr, MARGIN = 1, vec, "/"), na.rm = T)
    return(norm)
  } else {
    norm <- matrix(rep(NA, ncol(data)), nrow = 1)
  }
}




#' Calculate reaction score for enzyme complex(and relationship only)
#'
#' @param x index
#' @param data input gene expression data
#' @param list enzyme complex reaction list
#' @param gene_num file of number of times one gene has involved in all pathways
#'
#' @return Calculated enzyme complex list scores
#'
#'
#' @examples
calculate_simple_comeplex_score <- function(x, data, list, gene_num) {
  data_feature <- list[[x]][list[[x]] %in% rownames(data)]
  if (length(data_feature) > 0) {
    vec <- gene_num[data_feature, ]$V1
    expr <- as.matrix(data[data_feature, , drop = F])
    norm <- apply((1/vec) * expr, 2, min, na.rm = T)
    return(norm)
  } else {
    norm <- matrix(rep(NA, ncol(data)), nrow = 1)
  }
}


#' calculate reaction score for complicated complex(with both 'and' and 'or' relationship)
#'
#' @param x index
#' @param data input gene expression data
#' @param gene_num file of number of times one gene has involved in all pathways
#'
#' @return Calculated complicated enzyme complex list scores
#' @importFrom stringi stri_detect_fixed
#' @importFrom stringr str_extract_all
#'
#' @examples
calculate_multi_comp <- function(x, data, gene_num) {
  reaction <- multi_comp[[x]]
  # determine whether it is or and? within the bracket
  com <- stri_detect_fixed(str_extract_all(reaction, "\\([^()]+\\)")[[1]],
                           "or")
  c <- gsub("\\)", "", gsub("\\(", "", str_extract_all(reaction, "\\([^()]+\\)")[[1]]))
  if (unique(com) == T) {
    newiso <- lapply(lapply(c, function(x) {
      unlist(strsplit(x, "or"))
    }), function(x) {
      trimws(x)
    })
    sum_score <- do.call(rbind, lapply(1:length(newiso), calculate_iso_score,
                                       data = data, list = newiso, gene_num = gene_num))
    feature <- trimws(unlist(strsplit(reaction, "and"))[!stri_detect_fixed(unlist(strsplit(reaction,
                                                                                           "and")), "or")])
    data_feature <- feature[feature %in% rownames(data)]
    vec <- gene_num[data_feature, ]$V1
    expr <- as.matrix(data[data_feature, , drop = F])
    whole_score <- rbind((1/vec) * expr, sum_score)
    norm <- apply(whole_score, 2, min, na.rm = T)
  } else if (unique(com) == FALSE) {
    # trim the white space for iso enzyme
    newiso <- lapply(lapply(c, function(x) {
      unlist(strsplit(x, "and"))
    }), function(x) {
      trimws(x)
    })
    sum_score <- do.call(rbind, lapply(1:length(newiso), calculate_simple_comeplex_score,
                                       data = data, list = newiso, gene_num = gene_num))
    feature <- trimws(unlist(strsplit(reaction, "or"))[!stri_detect_fixed(unlist(strsplit(reaction,
                                                                                          "or")), "and")])
    data_feature <- feature[feature %in% rownames(data)]
    vec <- gene_num[data_feature, ]$V1
    expr <- as.matrix(data[data_feature, , drop = F])
    upper_score <- (1/vec) * expr
    whole_score <- rbind(upper_score, sum_score)
    norm <- colSums(whole_score, na.rm = T)
  }
}



#' Normalize scores
#'
#' @param x index
#' @param ...
#'
#' @return
#'
#'
#' @examples
#'
stdize = function(x, ...) {(x ) / (max(x, ...))}


#' Calculate metabolic reaction scores (MRAS) for 13082 reactions
#'
#' @param data gene expression data.1.The gene expression matrix should be gene by sample matrix where row names are human gene names (gene symbols),
#' and column names should be sample names. Please note that METAFlux does not support other gene IDs.
#' 2.The input gene expression matrix should be normalized (e.g., log-transformed, etc.) before using METAFlux.
#' METAflux will not perform any normalization on expression data.
#' 3.Gene expression data cannot have negative values.
#'
#' @return
#' @export
#'
#' @examples  calculate_reaction_score(bulk_test_example)
calculate_reaction_score <- function(data) {
  if (sum(data < 0) > 0)
    stop("Expression data needs to be all positive")
  # make sure features are present
  features <- rownames(data)
  if (sum(features %in% rownames(gene_num)) == 0)
    stop("Requested gene names cannot be found. Rownames of input data should be human gene names.Please check the rownames of input data.")  #change
  message(paste0(round(sum(features %in% rownames(gene_num))/3625 * 100,
                       3), "% metabolic related genes were found......"))
  gene_num <- METAFlux:::gene_num
  Hgem <- METAFlux:::Hgem
  iso <- METAFlux:::iso
  multi_comp <- METAFlux:::multi_comp
  simple_comp <- METAFlux:::simple_comp
  message("Computing metabolic reaction activity scores......")
  core <- do.call(rbind, lapply(1:length(iso), calculate_iso_score, data = data,
                                list = iso, gene_num = gene_num))
  core2 <- do.call(rbind, lapply(1:length(simple_comp), calculate_simple_comeplex_score,
                                 data, list = simple_comp, gene_num = gene_num))
  core3 <- do.call(rbind, lapply(1:length(multi_comp), calculate_multi_comp,
                                 data = data, gene_num = gene_num))
  message("Preparing for score matrix......")
  rownames(core) <- names(iso)
  rownames(core2) <- names(simple_comp)
  rownames(core3) <- names(multi_comp)
  big_score_matrix <- rbind(core, core2, core3)
  big_score_matrix <- apply(big_score_matrix, 2, stdize, na.rm = T)
  # handle NAN
  big_score_matrix[is.na(big_score_matrix)] <- 0
  empty_helper <- as.data.frame(Hgem$Reaction)
  colnames(empty_helper) <- "reaction"
  Final_df <- merge(empty_helper, big_score_matrix, all.x = T, by.x = 1,
                    by.y = 0)
  Final_df[is.na(Final_df)] <- 1
  rownames(Final_df) <- Final_df$reaction
  Final_df$reaction <- NULL
  Final_df <- Final_df[Hgem$Reaction, , drop = F]
  if (all.equal(rownames(Final_df), Hgem$Reaction)) {
    message("Metabolic reaction activity scores successfully calculated \n")
  } else {
    message("Calculation not reliable Check input data format \n")
  }
  Final_df[which(Hgem$LB == 0 & Hgem$UB == 0), ] <- 0  #this is biomass reaction, do not use
  return(Final_df)
}



