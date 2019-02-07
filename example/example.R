# Clean environment
rm(list = ls(all = TRUE))
options(warn = -1)



# Install
# if(!require(devtools)) {install.packages(devtools)}
# devtools::install_github("ShanSabri/deconR")
library(deconR)



# Load ENCODE data
data(encode_signal)
data(encode_expr)
data_split <- readRDS(file.path("~/Dropbox/ErnstLab/decon_mouse", "data", "DATA-SPLIT.rds"))



# Load scRNA data
scrna <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/TC1-RAW-FILTERED-DGE.rds")
scrna_pheno <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/TC1-PHENODATA.rds")
identical(row.names(scrna_pheno), names(scrna))



# Create artifical scRNA bulk samples
clusters_of_interest <- c(1, 23)
scrna_aggr <- make_scrna_aggregrate(scrna, labels = as.vector(scrna_pheno$Cluster), normalize = TRUE)
scrna_aggr <- scrna_aggr[row.names(encode_expr), paste("Cluster", clusters_of_interest, sep = "_")]



# Subset to differential genes (OPTIONAL)
diff_exp <- read.table(gzfile("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/Clustering/Drop7/2016-12-21-final-data/diffexp/Drop7.negbinom.cluster.v.all.sig.genes.tsv.gz"), header=TRUE)
diff_exp_genes <- as.vector(subset(diff_exp, cluster %in% c("C1", "C23"))$gene)
encode_expr <- encode_expr[row.names(encode_expr) %in% diff_exp_genes, ]
scrna_aggr <- scrna_aggr[row.names(scrna_aggr) %in% diff_exp_genes, ]
identical(row.names(scrna_aggr), row.names(encode_expr))



# # How well do these data correlate (ENCODE vs aggregrated scRNA)?
# corrs <- list(p_cor = cor(cbind.data.frame(scrna_aggr, encode_expr), method = "pearson")[,1:2],
#               s_cor = cor(cbind.data.frame(scrna_aggr, encode_expr), method = "spearman")[,1:2])
# heatmaps <- lapply(seq_along(corrs), function(x){
#   ph <- pheatmap::pheatmap(corrs[[x]],
#                      display_numbers = TRUE,
#                      cluster_cols = FALSE,
#                      main = toupper(names(corrs)[x]),
#                      filename = file.path("~/Dropbox/ErnstLab/deconR/example", paste0(toupper(names(corrs)[x]), "-COR-BULK-W-SCRNA.png")),
#                      height = 14,
#                      width = 5)
# })



# Create in silico convolved signal track
MEF <- read_signal("~/Dropbox/ErnstLab/deconvolution_old/CHIP/signal/full/200bp-windows/liftOver-mm10/MEF_K27AC.sorted.wig", format = "wig", chr = "chr16")
ESC <- read_signal("~/Dropbox/ErnstLab/deconvolution_old/CHIP/signal/full/200bp-windows/liftOver-mm10/ES_K27AC.sorted.wig", format = "wig", chr = "chr16")
cts <- list(ct1 = MEF, ct2 = ESC)
convolved <- make_convolved(cts, log2_scale = TRUE)
write_bw(convolved, out = file.path(getwd(), "example", "tracks", "CONVOLVED-LOG2.bw")) # already log2 scale
write_bw(convolved, out = file.path(getwd(), "example", "tracks", "CONVOLVED.bw"), unlog2_scale = TRUE)
write_bw(MEF, out = file.path(getwd(), "example", "tracks", "MEF.bw"))
write_bw(MEF, out = file.path(getwd(), "example", "tracks", "MEF-LOG2.bw"), log2_scale = TRUE)
write_bw(ESC, out = file.path(getwd(), "example", "tracks", "ESC.bw"))
write_bw(ESC, out = file.path(getwd(), "example", "tracks", "ESC-LOG2.bw"), log2_scale = TRUE)



# Create BASELINE model (50/50 split of cell types above)
baseline_model <- convolved
baseline_model$score <- (2^baseline_model$score - 1) / length(cts)
write_bw(baseline_model, out = file.path(getwd(), "example", "tracks", "BASELINE-MODEL.bw"))
write_bw(baseline_model, out = file.path(getwd(), "example", "tracks", "BASELINE-MODEL-LOG2.bw"), log2_scale = TRUE)



# Cluster loci to reduce dimensionality *only cluster training data*
n_clusters <- 500
# clusters <- cluster_loci(encode_signal[, data_split$TRAINING], k = n_clusters, iter = 30)
# saveRDS(clusters, file.path(getwd(), "example", "CLUSTERS.rds"), compress = TRUE)
clusters <- readRDS(file.path(getwd(), "example", "CLUSTERS.rds"))

# Train across-cell type models *one model per cluster*
# across_model_output <- model_across_celltypes(clusters, scrna_aggr, encode_signal, encode_expr,
#                                               data_split, genes = diff_exp_genes, seed = 42,
#                                               output_dir = "~/Dropbox/ErnstLab/deconR/example/ACROSS-MODEL-PLOTS")
# saveRDS(across_model_output, file.path(getwd(), "example", "ACROSS-MODEL-OUTPUT.rds"), compress = TRUE)
across_model_output <- readRDS(file.path(getwd(), "example", "ACROSS-MODEL-OUTPUT.rds"))
cluster_split_proportions <- compute_split_proportions(across_model_output)
across_model_decon <- compute_across_model_decon(cluster_split_proportions, cluster_labels = clusters,
                                                 convolved = convolved, verbose = TRUE)




############## ANALYSIS ############################ ANALYSIS ############################ ANALYSIS ##############
############## ANALYSIS ############################ ANALYSIS ############################ ANALYSIS ##############
############## ANALYSIS ############################ ANALYSIS ############################ ANALYSIS ##############

## PLOT CLUSTER AVERAGES FOR EACH CLUSTER AND COMPARE WITH BASELINE AND TRUE
scRNA_loci_clust_avgs <- do.call(rbind.data.frame, lapply(seq_along(across_model_output), function(y){
  predictions <- as.data.frame(t(across_model_output[[y]]$cluster_prediction_mu))
  predictions$CLUSTER <- as.numeric(predictions$CLUSTER)
  predictions
}))

clusters_list <- split(clusters, ~CLUSTER)
bulk_loci_clust_avgs <- do.call(rbind.data.frame, lapply(seq_along(clusters_list), function(x){
  message(x)
  strt <- start(clusters_list[[x]])
  true_clust_avgs <- data.frame(MEF = colMeans(data.frame(mcols(MEF[start(MEF) %in% strt]))),
                                ESC = colMeans(data.frame(mcols(ESC[start(ESC) %in% strt]))),
                                CLUSTER = as.numeric(x))
}))

clusters_list <- split(clusters, ~CLUSTER)
null_loci_clust_avgs <- do.call(rbind.data.frame, lapply(seq_along(clusters_list), function(x){
  message(x)
  strt <- start(clusters_list[[x]])
  null_clust_avgs <- data.frame(BASELINE = colMeans(data.frame(mcols(baseline_model[start(baseline_model) %in% strt]))),
                                CLUSTER = as.numeric(x))
}))

to_plot = Reduce(function(...) merge(..., all=T), list(scRNA_loci_clust_avgs, bulk_loci_clust_avgs, null_loci_clust_avgs))

library(ggplot2)
P1 <- ggplot(to_plot, aes(log2(MEF + 1), log2(BASELINE + 1))) +
  geom_text(alpha = 0.7, aes(label = CLUSTER), size = 2) +
  geom_smooth(method = "lm", colour = 'blue') +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank()) +
  labs(title = "Chr16: TRUE vs NULL MODEL", subtitle = paste0("Pearson Corr: ", round(cor(log2(to_plot$BASELINE + 1), log2(to_plot$MEF + 1), method = "pearson"), 5)))
P2 <- ggplot(to_plot, aes(log2(MEF + 1), Cluster_1)) +
  geom_text(alpha = 0.7, aes(label = CLUSTER), size = 2) +
  geom_smooth(method = "lm", colour = 'blue') +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank()) +
  labs(title = "Chr16: TRUE vs DECON MODEL", subtitle = paste0("Pearson Corr: ", round(cor(to_plot$Cluster_1, log2(to_plot$MEF + 1), method = "pearson"), 5)))
P <- cowplot::plot_grid(P1, P2, labels = "AUTO") +
  ggsave(file.path(getwd(), "example", "MEF-VS-NULL-VS-DECON.png"), height = 5, width = 9)

P3 <- ggplot(to_plot, aes(log2(ESC + 1), log2(BASELINE + 1))) +
  geom_text(alpha = 0.7, aes(label = CLUSTER), size = 2) +
  geom_smooth(method = "lm", colour = 'red') +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank()) +
  labs(title = "Chr16: TRUE vs NULL MODEL", subtitle = paste0("Pearson Corr: ", round(cor(log2(to_plot$BASELINE + 1), log2(to_plot$ESC + 1), method = "pearson"), 5)))
P4 <- ggplot(to_plot, aes(log2(ESC + 1), Cluster_23)) +
  geom_text(alpha = 0.7, aes(label = CLUSTER), size = 2) +
  geom_smooth(method = "lm", colour = 'red') +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank()) +
  labs(title = "Chr16: TRUE vs DECON MODEL", subtitle = paste0("Pearson Corr: ", round(cor(to_plot$Cluster_23, log2(to_plot$ESC + 1), method = "pearson"), 5)))
P <- cowplot::plot_grid(P3, P4, labels = "AUTO") +
  ggsave(file.path(getwd(), "example", "ESC-VS-NULL-VS-DECON.png"), height = 5, width = 9)


############## TO DO ############################ TO DO ############################ TO DO ############################ TO DO ##############
############## TO DO ############################ TO DO ############################ TO DO ############################ TO DO ##############
############## TO DO ############################ TO DO ############################ TO DO ############################ TO DO ##############
