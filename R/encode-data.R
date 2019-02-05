#' Mouse ENCODE H3K27ac raw signal values
#'
#' Data from 87 mouse cell types source from ENCODE contains scaled (log2) H3K27ac
#' signal indensity values with a genomic resolution of 200 bps covering all
#' of chr16 (491038 non-overlapping 200bp bins).
#'
#' @docType data
#'
#' @usage data(encode_signal)
#'
#' @format An object of class \code{"GRanges"}; see \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}.
#'
#' @keywords datasets
#'
#' @references Yue et al. (2014) Nature 515:355-364
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/25409824}{PubMed})
#'
#' @source \href{http://www.mouseencode.org}{ENCODE Mouse Project}
#'
#' @examples
#' data(encode_signal)
"encode_signal"



#' Mouse ENCODE RNA-seq gene expression values
#'
#' Data from 87 mouse cell types source from ENCODE contains RNA-seq depth-normalized
#' TPM values with replicates aggregrated for 14539 genes.
#'
#' @docType data
#'
#' @usage data(encode_expr)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords datasets
#'
#' @references Yue et al. (2014) Nature 515:355-364
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/25409824}{PubMed})
#'
#' @source \href{http://www.mouseencode.org}{ENCODE Mouse Project}
#'
#' @examples
#' data(encode_expr)
"encode_expr"
