#' Predicts ChIP-seq profiles from gene expression by learned relationships between RNA and
#' ChIP within ENCODE data
#'
#' @param clusters A \code{data.frame} containing genomic bins associated with a given
#' cluster label
#' @param encode_signal An object of class \code{"GRanges"} containing normalized
#' ChIP-seq signal indensity values sourced from Mouse ENCODE
#' @param encode_expr An object of class \code{"data.frame"} containing normalized
#' RNA-seq gene expression information source from Mouse ENCODE
#' @param data_split a list of cell types to be used for TRAINING and VALIDATION
#' @param output_dir optional; if specified, then function will output model/loci-specific
#' scatterplots with modeling statistics to the specified directory
#' @param scrna_aggr \code{data.frame} containing aggregrated single cell gene counts
#' @param genes a character vector of genes to be used as features to the model
#' @param seed set seed for model reproducibility
#'
#' @return a list of lists contaning predicted and true, cluster-averaged predictions for each cell
#' type. Also included are average prediction from using scRNA aggregrates.
#' @import ggplotify FNN
#' @importFrom ggplotify as.ggplot
#' @importFrom FNN knn.reg
#' @importFrom cowplot plot_grid save_plot
#' @export
model_across_celltypes <- function(clusters, scrna_aggr, genes = NULL,
                                   encode_signal, encode_expr, data_split,
                                   output_dir = NULL, seed = 42){
  start_time <- Sys.time()
  n_clusters <- max(mcols(clusters)$CLUSTER)

  if(is.null(genes)){
    genes <- row.names(encode_expr)
  }

  output <- lapply(seq_along(n_clusters), function(y) {
    set.seed(seed)
    print(paste(paste("Cluster", y), paste0(round(difftime(Sys.time(), start_time, units='mins'), 4), " mins elapsed"), sep=": "))
    bins <- subset(clusters, CLUSTER == y)

    # train
    Y_train <- encode_signal[start(encode_signal) %in% start(bins), data_split$TRAINING]
    X_train <- encode_expr[row.names(encode_expr) %in% genes, data_split$TRAINING]
    train_model <- FNN::knn.reg(train = t(X_train), test = t(X_train), y = colMeans(as.data.frame(mcols(Y_train))), k = 2)$pred
    train_results <- cbind.data.frame(Cell_type = names(X_train),
                                      Predicted = train_model,
                                      True = colMeans(as.data.frame(mcols(Y_train))))

    # validation
    Y_val <- encode_signal[start(encode_signal) %in% start(bins), data_split$VALIDATION]
    X_val <- encode_expr[row.names(encode_expr) %in% genes, data_split$VALIDATION]
    validation_model <- FNN::knn.reg(train = t(X_train), test = t(X_val), y = colMeans(as.data.frame(mcols(Y_train))), k = 2)$pred
    validation_results <- cbind.data.frame(Cell_type = names(X_val),
                                           Predicted = validation_model,
                                           True = colMeans(as.data.frame(mcols(Y_val))))

    # plot
    if(!is.null(output_dir)){
      train_plot <- plot_model_results(train_results, bins, title = "Predicting TRAIN:TRAIN")
      validation_plot <- plot_model_results(validation_results, bins, title = "Predicting TRAIN:VALIDATION")
      heatmap <- plot_loci_heatmap(encode_signal, bins, data_split, min = -3, max = 3)
      top_row <- cowplot::plot_grid(train_plot, validation_plot, labels = "AUTO")
      bottom_row <- cowplot::plot_grid(ggplotify::as.ggplot(heatmap), labels = c('C'), align = 'h')
      cowplot::plot_grid(top_row, bottom_row, ncol=1, rel_heights = c(1, 1.5)) +
        cowplot::save_plot(file.path(output_dir, paste0("ACROSS-MODEL-CLUSTER-", y, ".png")), height = 12, width = 14)
    }

    # test
    test_model <- FNN::knn.reg(train = t(X_train), test = t(scrna_aggr[row.names(scrna_aggr) %in% genes, ]),
                               y = colMeans(as.data.frame(mcols(Y_train))), k = 2)$pred
    names(test_model) <- names(scrna_aggr)
    test_model <- c(test_model, CLUSTER = as.integer(y))

    results <- list(train_results = train_results,
                    validation_results = validation_results,
                    cluster_prediction_mu = test_model)

    return(results)
  })
  names(output) <- paste("Cluster", seq(n_clusters), sep = "_")

  return(output)
}
