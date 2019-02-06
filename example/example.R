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



# How well do these data correlate (ENCODE vs aggregrated scRNA)?
corrs <- list(p_cor = cor(cbind.data.frame(scrna_aggr, encode_expr), method = "pearson")[,1:2],
              s_cor = cor(cbind.data.frame(scrna_aggr, encode_expr), method = "spearman")[,1:2])
heatmaps <- lapply(seq_along(corrs), function(x){
  ph <- pheatmap::pheatmap(corrs[[x]],
                     display_numbers = TRUE,
                     cluster_cols = FALSE,
                     main = toupper(names(corrs)[x]),
                     filename = file.path("~/Dropbox/ErnstLab/deconR/example", paste0(toupper(names(corrs)[x]), "-COR-BULK-W-SCRNA.png")),
                     height = 14,
                     width = 5)
})



# Create in silico convolved signal track
MEF <- read_signal("~/Dropbox/ErnstLab/deconvolution_old/CHIP/signal/full/200bp-windows/liftOver-mm10/MEF_K27AC.sorted.wig", format = "wig", chr = "chr16")
ESC <- read_signal("~/Dropbox/ErnstLab/deconvolution_old/CHIP/signal/full/200bp-windows/liftOver-mm10/ES_K27AC.sorted.wig", format = "wig", chr = "chr16")
convolved <- make_convolved(list(ct1 = MEF, ct2 = ESC), log2_scale = TRUE)



# Cluster loci to reduce dimensionality *only cluster training data*
n_clusters <- 500
clusters <- cluster_loci(encode_signal[, data_split$TRAINING], k = n_clusters)



# Train across-cell type models *one model per cluster*
across_model_output <- model_across_celltypes(clusters, scrna_aggr, encode_signal, encode_expr,
                                              data_split, genes = diff_exp_genes, seed = 42,
                                              output_dir = "~/Dropbox/ErnstLab/deconR/example/ACROSS-MODEL-PLOTS")










# NOTES
# save(encode_signal, file = file.path("~/Dropbox/ErnstLab/deconR/data", "encode_signal.RData"), compress='xz')
# encode_signal_metadata <- cbind.data.frame(chr = seqnames(encode_signal), ranges(encode_signal))
