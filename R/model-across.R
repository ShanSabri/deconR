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
  clusters_list <- split(clusters, f = mcols(clusters)$CLUSTER)
  n_clusters <- length(clusters_list)

  if(is.null(genes)){
    genes <- row.names(encode_expr)
  }

  output <- lapply(seq_along(clusters_list), function(y) {
    set.seed(seed)
    print(paste(paste("Cluster", y), paste0(round(difftime(Sys.time(), start_time, units='mins'), 4), " mins elapsed"), sep=": "))

    bins <- clusters_list[[y]]

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
      plot <- cowplot::plot_grid(top_row, bottom_row, ncol=1, rel_heights = c(1, 1.5))
      cowplot::save_plot(file.path(output_dir, paste0("ACROSS-MODEL-CLUSTER-", y, ".png")), plot, base_height = 12, base_width = 14)
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


#' Compute cell type ratios to split convolved signal based on across-cell type
#' predictions
#'
#' @param x a list containing cluster-specific, average loci intensity predictions,
#' where each element in the list corresponds to a cluster.
#'
#' @return a list containing predictions that are convered to split proportions
#' @export
compute_split_proportions <- function(x){
  cluster_split_proportions <- lapply(seq_along(x), function(y){
    cluster <- names(x)[y]
    predictions <- as.data.frame(t(x[[y]]$cluster_prediction_mu))
    predictions$CLUSTER <- NULL
    total_prediction <- sum(predictions)
    return(predictions / total_prediction)
  })
  names(cluster_split_proportions) <- names(x)

  return(cluster_split_proportions)
}



#' Compute deconvolution from across-cell type modeling based on the
#' cell type ratios derived from \code{compute_split_proportions}
#'
#' @param x A \code{GRanges} object returned from the function \code{compute_split_proportions}
#' @param cluster_labels A \code{data.frame} or \code{GRanges} object containing
#' genomic bins associated with a given cluster label
#' @param convolved a \code{GRanges} object with convolved ChIP-seq signal data
#' @param verbose boolean indicating whether to flush output to console
#'
#' @return a \code{GRanges} object containing metadata of across-cell type
#' deconvolution
#' @importFrom GenomicRanges GRanges
#' @export
compute_across_model_decon <- function(x, cluster_labels, convolved, verbose = TRUE){
  across_decon <- do.call(c, lapply(seq_along(x), function(y){
    clust_id <- names(x)[y]
    if(verbose == TRUE) message(clust_id)

    clust_bins <- cluster_labels[cluster_labels$CLUSTER == y, ]
    convolved_bins <- convolved[start(convolved) %in% start(clust_bins)]
    split_proportion <- as.list(x[[y]])

    cluster_decon <- sapply(split_proportion, function(z){
      mcols(convolved_bins)$score * z},
      USE.NAMES = TRUE)

    if(!is.null(dim(cluster_decon))){
      gr <- GenomicRanges::GRanges(seqnames = seqnames(convolved_bins),
                                   ranges = ranges(convolved_bins),
                                   strand = NULL, seqinfo = seqinfo(convolved_bins),
                                   cluster_decon)
    } else {
      gr <- GenomicRanges::GRanges(seqnames = seqnames(convolved_bins),
                                   ranges = ranges(convolved_bins),
                                   strand = NULL, seqinfo = seqinfo(convolved_bins),
                                   t(as.matrix(cluster_decon)))
    }
    gr
  }))

  return(sort(across_decon))
}
