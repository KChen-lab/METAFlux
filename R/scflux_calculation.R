#' sc-rna seq flux calculation
#'
#' @param num_cell the number of cell types or clusters
#' @param fraction fraction of each cell types. Fractions need to sum up to 1
#' @param fluxscore calculated metabolic activity score(MRAS)
#' @param medium Medium profile
#' @importFrom Seurat as.sparse
#' @import utils
#'
#' @return Flux score for single cell data
#' @export
#'
#' @examples
compute_sc_flux <- function(num_cell, fraction, fluxscore, medium) {
  if (sum(fraction) != 1)
    stop("Sum of fractions must be eqaul to 1")
  if (length(fraction) != num_cell)
    stop("Number of cell clusters does not match with the length of fraction")
  # perform preparation
  Hgem <- METAFlux:::Hgem
  mat <- Hgem$S
  reaction_name <- Hgem$Reaction
  names(reaction_name) <- NULL
  A_combined <- METAFlux:::A_combined
  D <- as.sparse(matrix(0, 8378, 13082))
  candi_list <- lapply(rapply(list(mat, lapply(1:2, function(x) D)),
                              enquote, how = "unlist"), eval)
  matrix_construct <- diag(1, ncol = num_cell, nrow = num_cell)
  matrix_construct[matrix_construct == 0] <- 2
  celltype_matrix <- NULL
  A_matrix <- NULL
  for (i in 1:num_cell) {
    A_matrix[[i]] <- A_combined
    celltype_matrix[[i]] <- do.call(cbind, candi_list[matrix_construct[,
                                                                       i]])
  }
  message("Preparing for TME S matrix.....")
  final_s <- rbind(do.call(cbind, A_matrix), do.call(rbind, celltype_matrix))
  whole3 <- rbind(as.sparse(diag(-1, 1648, 1648)), as.sparse(matrix(0,
                                                                    num_cell * nrow(mat), 1648)))
  final_s <- cbind(final_s, whole3)
  external <- paste("external_medium", reaction_name[which(Hgem$pathway ==
                                                             "Exchange/demand reactions")])
  reaction_name[which(Hgem$pathway == "Exchange/demand reactions")] <- paste0("internal_medium ",
                                                                              reaction_name[which(Hgem$pathway == "Exchange/demand reactions")])
  cell_reaction <- NULL
  for (i in 1:num_cell) {
    cell_reaction[[i]] <- paste(paste("celltype", i), reaction_name)
  }
  construct_reaction_names <- c(unlist(cell_reaction), external)
  P1 <- as.sparse(diag(1, ncol(final_s), ncol(final_s)))
  message("S matrix completed......")
  flux_vector <- list()
  message("Compute metabolic flux......")
  Seq <- seq(1, ncol(fluxscore), by = num_cell)
  pb <- txtProgressBar(0, length(Seq), style = 3)
  for (t in Seq) {
    setTxtProgressBar(pb, match(t, Seq))
    score <- fluxscore[, c(t:(t + num_cell - 1))]
    P <- as.sparse(diag(c(unlist(Map(rep, fraction, 13082)), rep(1,
                                                                 1648)), ncol(final_s), ncol(final_s)))
    fraction_finals <- final_s %*% P
    q <- rep(0, ncol(final_s))
    q[seq(13015, ncol(final_s), by = 13082)] <- -10000 * fraction
    A <- as.sparse(rbind(fraction_finals, P1))
    ras <- c(as.vector(unlist(score[, 1:num_cell])), rep(1, 1648))
    origlb <- c(rep(Hgem$LB, num_cell), rep(-1, 1648))
    origlb[rep(Hgem$rev == 1, num_cell)] <- (-ras[rep(Hgem$rev == 1,
                                                      num_cell)])
    origlb[rep(Hgem$rev == 0, num_cell)] <- 0
    origub <- ras
    origlb[tail(1:ncol(final_s), 1648)] <- 0
    matches <- sapply(medium$reaction_name, function(x) intersect(which(stri_detect_fixed(construct_reaction_names,
                                                                                               x)), tail(1:ncol(final_s), 1648)))
    origlb[matches] <- -1
    l <- c(rep(0, nrow(final_s)), origlb)
    u <- c(rep(0, nrow(final_s)), origub)
    settings <- osqpSettings(max_iter = 1000000L, eps_abs = 1e-04,
                             eps_rel = 1e-04, verbose = FALSE, adaptive_rho_interval = 50)
    model <- osqp(P, q, A, l, u, settings)
    # Solve problem
    res <- model$Solve()
    flux_vector <- append(flux_vector, list(res$x))
  }
  close(pb)
  flux_vector <- as.data.frame(do.call(cbind, flux_vector))
  rownames(flux_vector) <- construct_reaction_names
  return(flux_vector)

}

