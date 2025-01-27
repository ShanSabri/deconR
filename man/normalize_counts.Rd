% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{normalize_counts}
\alias{normalize_counts}
\title{Depth normalize RNA-seq gene expression matrix}
\usage{
normalize_counts(x, sf = 10000)
}
\arguments{
\item{x}{\code{data.frame} containing gene (row) by cell (column) RNA-seq
quantification}

\item{sf}{size factor used for scaling}
}
\value{
depth normalized gene expression matrix
}
\description{
Depth normalize RNA-seq gene expression matrix
}
