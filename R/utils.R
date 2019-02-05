#' Depth normalize RNA-seq gene expression matrix
#'
#' @param x \code{data.frame} containing gene (row) by cell (column) RNA-seq
#' quantification
#' @param sf size factor used for scaling
#'
#' @return depth normalized gene expression matrix
#' @export
normalize_counts <- function(x, sf = 10000){
  log((sweep(x, 2, apply(x, 2, function(x) sum(x)), "/") * sf) + 1)
}


#' Average single cell clusters to one artifical bulk sample aggregrate
#'
#' @param x \code{data.frame} containing gene (row) by cell (column) RNA-seq
#' quantification
#' @param labels a vector of cell labels used for aggregration
#' @param normalize boolean indicating whether to normalize
#'
#' @return \code{data.frame} containing aggregrated single cell gene counts
#' @export
make_scrna_aggregrate <- function(x, labels, normalize = TRUE){
  scrna_aggr <- do.call(cbind.data.frame, lapply(sort(unique(labels)), function(l) {
    cells <- colnames(x[, labels %in% l])
    aggr <- rowSums(x[, colnames(x) %in% cells, drop = FALSE])
  }))
  names(scrna_aggr) <- paste("Cluster", sort(unique(labels)), sep = "_")

  if(normalize == TRUE){
    scrna_aggr <- normalize_counts(scrna_aggr)
  }

  return(scrna_aggr)
}
