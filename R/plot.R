#' Plot a scatterplot of model predicted values against true values
#'
#' @param x a \code{data.frame} contains model predictions
#' @param bins a \code{GRanges} object containing the genomic bins within a
#' given cluster of loci. This object is used for cluster and n_loci lookup.
#' @param title plot title/identifier
#'
#' @return a \code{ggplot} scatterplot object
#' @import ggplot2 ggrepel
#' @importFrom stats cor lm predict sd
#' @export
plot_model_results <- function(x, bins, title = ""){
  mse <- mean((x$True - x$Predicted)^2)
  fit <- lm(True ~ Predicted, data = x)
  interval <- predict(fit, interval="confidence")
  x$inside <- ifelse(x$True < interval[,"upr"] & x$True > interval[,"lwr"], "", as.character(x$Cell_type))
  P <- ggplot(x, aes(Predicted, True)) +
    geom_point() +
    geom_point(data = subset(x, Cell_type %in% c("embryonicfibroblast_13.5week")), aes(Predicted, True), colour = 'blue') +
    geom_point(data = subset(x, Cell_type %in% c("ES.Bruce4_unknown")), aes(Predicted, True), colour = 'red') +
    geom_smooth(method = "lm", colour = "forestgreen", se=TRUE) + geom_text_repel(aes(label = inside), size=1.5, segment.size = 0.2) +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank()) + geom_abline(slope=1, intercept=0, colour="black", linetype="dotted", size=0.5) +
    labs(title = title,
         subtitle = paste0("Cluster ", unique(bins$CLUSTER), " with ", length(bins), " loci",
                           "\nPearson Corr: ", round(cor(x$Predicted, x$True, method = "pearson"), 5),
                           "\nMSE: ", round(mse, 5)))

  return(P)
}


#' Plot a heatmap of signal intensity values for a given cluster of loci
#'
#' @param x a \code{GRanges} object containing ChIP-seq signal intensity values
#' @param bins a \code{GRanges} object containing the genomic bins within a
#' given cluster of loci. This object is used for subsetting the original signal
#' \code{GRanges} object
#' @param data_split a list of cell types to be used for TRAINING and VALIDATION
#' @param min floor of min Z-score
#' @param max ceiling of max Z-score
#'
#' @return a \code{pheatmap} object containing heatmap data
#' @import pheatmap RColorBrewer
#' @importFrom grDevices colorRampPalette
#' @export
plot_loci_heatmap <- function(x, bins, data_split, min = -3, max = 3){
  bs <- seq(min, max, 0.2)
  ph_data <- as.data.frame(mcols(x[start(x) %in% start(bins)]))
  ph_data <- ph_data[, c(data_split$TRAINING, data_split$VALIDATION)]
  annot <- data.frame(row.names = c(data_split$TRAINING, data_split$VALIDATION),
                      SET = c(rep("TRAIN", length(data_split$TRAINING)), rep("VALIDATION", length(data_split$VALIDATION))),
                      AVG_SIGNAL = colMeans(ph_data))
  ph_data <- ph_data[apply(ph_data, 1, FUN = function(x) sd(x) != 0), ]
  P <- pheatmap::pheatmap(ph_data,
                          show_rownames = FALSE,
                          scale = "row",
                          annotation_col = annot,
                          annotation_names_col = FALSE,
                          cluster_rows = ifelse(nrow(ph_data) > 1, TRUE, FALSE),
                          color = colorRampPalette((brewer.pal(9,"BrBG")))(length(bs)),
                          breaks = bs,
                          na_col = 'black',
                          annotation_names_row = FALSE)

  return(P)
}
