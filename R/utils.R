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


#' Read in ChIP-seq signal intensity values
#'
#' @param x file name of ChIP-seq data
#' @param format file format of ChIP-seq data
#' @param genome If not missing, should be “wig” or “bw” (case insensitive).
#' @param chr optional; specify to subset to given chromosome
#'
#' @return a \code{GRanges} object with ChIP-seq signal data
#' @import rtracklayer GenomeInfoDb
#' @importFrom rtracklayer import.bw
#' @importFrom rtracklayer import.wig
#' @importFrom GenomeInfoDb keepSeqlevels
#' @export
read_signal <- function(x, format = "bw", genome = "mm10", chr = NULL) {
  if (format == "bw") {
    f <- rtracklayer::import.bw(x, format = "bw")
  } else if (format == "wig") {
    f <- rtracklayer::import.wig(x, format = "wig", genome = genome)
  } else {
    message("INVALID FORMAT")
  }

  f <- GenomeInfoDb::keepSeqlevels(f, paste0("chr", c(seq(1, 19), "X")), pruning.mode = "coarse")

  if (!is.null(chr)) {
    f <- f[seqnames(f) == chr]
  }

  return(f)
}


#' Generate an in silico convolved ChIP-seq signal track from known ChIP-seq data
#'
#' @param x a list of \code{GRanges} objects with ChIP-seq signal data from different cell types to
#' to be convolved
#' @param log2_scale boolean indicating whether to log2 scale (with pseudocount of 1)
#'
#' @return a \code{GRanges} object with convolved ChIP-seq signal data
#' @export
make_convolved <- function(x, log2_scale = TRUE){
  n_cts <- length(x)
  score <- do.call(cbind.data.frame, lapply(x, function(y){ mcols(y)$score * (1/n_cts)}))
  convolved_score <- rowSums(score)

  if(log2_scale == TRUE){
    convolved_score <- log2(convolved_score + 1)
  }

  convolved <- x[[1]]
  convolved$score <- convolved_score

  return(convolved)
}


#' Cluster genomic loci via accelerated K-Means clustering
#'
#' @param x a \code{GRanges} object with ChIP-seq signal data with rows as
#' genomic bins and columns as cell types
#' @param k the number of k clusters
#' @param iter the number of clustering iterations (about 10 is typically sufficient)
#'
#' @return a \code{GRanges} object containing genomic coordinates tagged with a
#' cluster label
#' @import ClusterR GenomicRanges
#' @importFrom ClusterR KMeans_arma predict_KMeans
#' @importFrom GenomicRanges GRanges
#' @export
cluster_loci <- function(x, k = 500, iter = 30){
  km <- ClusterR::KMeans_arma(as.data.frame(mcols(x)),
                              clusters = k,
                              n_iter = iter,
                              seed_mode = "random_subset",
                              CENTROIDS = NULL, verbose = TRUE, seed = 42)
  pr <- ClusterR::predict_KMeans(as.data.frame(mcols(x)), km)
  gr <- GenomicRanges::GRanges(seqnames = seqnames(x), ranges = ranges(x), CLUSTER = as.vector(pr))

  return(gr)
}


#' Write a BigWig track from either a \code{GRanges} or \code{rtracklayer} object
#'
#' @param x The object to export
#' @param out The connection from which data is saved
#' @param log2_scale boolean indicating whether to log2 scale intensity values (pseudocount of 1)
#' @param unlog2_scale boolean indicating whether to remove log2 scaled intensity values
#' (assuming pseudocount of 1)
#'
#' @return nothing
#' @importFrom rtracklayer export.bw
#' @export
write_bw <- function(x, out, log2_scale = FALSE, unlog2_scale = FALSE){
  if(log2_scale == TRUE){
    mcols(x)$score <- log2(mcols(x)$score + 1)
  }

  if(unlog2_scale == TRUE){
    mcols(x)$score <- 2^mcols(x)$score - 1
  }

  rtracklayer::export.bw(object=x, con=out)
}
