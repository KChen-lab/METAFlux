
#' Final optimization step for flux calculation
#'
#' @param mras metabolic reaction activity scores
#' @param medium input medium file which indicates the nutrients available in the medium.
#' We provide 2 general mediums if you have no prior knowledge about your medium: cell line medium and human blood medium if prior knowledge of medium is not available.
#' Please see tutorial for more details.
#'
#' @return Calculated fluxes
#' @export
#' @import osqp
#'
#' @examples
compute_flux <- function(mras, medium) {
  message("Setting up for optimization.....")
  Hgem <- METAFlux:::Hgem
  P <- diag(1, ncol(Hgem$S), ncol(Hgem$S))
  q <- rep(0, ncol(Hgem$S))
  q[which(Hgem$Obj == 1)] <- -10000
  A <- as.matrix(rbind(Hgem$S, P))
  resu <- list()
  flux_vector <- list()
  message("Computing bulk RNA-seq flux.....")
  Seq <- seq(1, ncol(mras))
  pb <- txtProgressBar(0, length(Seq), style = 3)
  for (i in Seq) {
    setTxtProgressBar(pb, i)
    origlb <- Hgem$LB
    origlb[Hgem$rev == 1] <- -mras[, i][Hgem$rev == 1]
    origlb[Hgem$rev == 0] <- 0
    origlb <- origlb[, 1]
    origub <- mras[, i]
    origlb[Hgem$Reaction %in% Hgem$Reaction[which(Hgem$pathway == "Exchange/demand reactions")]] <- 0
    origlb[Hgem$Reaction %in% medium$reaction_name] <- -1
    l <- c(rep(0, nrow(Hgem$S)), origlb)
    u <- c(rep(0, nrow(Hgem$S)), origub)
    settings <- osqpSettings(max_iter = 1000000L, eps_abs = 1e-04,
                             eps_rel = 1e-04, adaptive_rho_interval = 50, verbose = FALSE)
    model <- osqp(P, q, A, l, u, settings)
    # Solve problem
    res <- model$Solve()
    resu[[i]] <- res
    flux_vector[[i]] <- res$x
  }
  close(pb)
  flux_vector <- do.call(cbind, flux_vector)
  colnames(flux_vector) <- colnames(mras)
  rownames(flux_vector) <- Hgem$Reaction
  return(flux_vector)
}
