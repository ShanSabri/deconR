% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_scrna_aggregrate}
\alias{make_scrna_aggregrate}
\title{Average single cell clusters to one artifical bulk sample aggregrate}
\usage{
make_scrna_aggregrate(x, labels, normalize = TRUE)
}
\arguments{
\item{x}{\code{data.frame} containing gene (row) by cell (column) RNA-seq
quantification}

\item{labels}{a vector of cell labels used for aggregration}

\item{normalize}{boolean indicating whether to normalize}
}
\value{
\code{data.frame} containing aggregrated single cell gene counts
}
\description{
Average single cell clusters to one artifical bulk sample aggregrate
}
