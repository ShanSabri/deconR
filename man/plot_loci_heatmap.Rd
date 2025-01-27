% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_loci_heatmap}
\alias{plot_loci_heatmap}
\title{Plot a heatmap of signal intensity values for a given cluster of loci}
\usage{
plot_loci_heatmap(x, bins, data_split, min = -3, max = 3)
}
\arguments{
\item{x}{a \code{GRanges} object containing ChIP-seq signal intensity values}

\item{bins}{a \code{GRanges} object containing the genomic bins within a
given cluster of loci. This object is used for subsetting the original signal
\code{GRanges} object}

\item{data_split}{a list of cell types to be used for TRAINING and VALIDATION}

\item{min}{floor of min Z-score}

\item{max}{ceiling of max Z-score}
}
\value{
a \code{pheatmap} object containing heatmap data
}
\description{
Plot a heatmap of signal intensity values for a given cluster of loci
}
