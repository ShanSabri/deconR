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



# Load scRNA data
scrna <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/TC1-RAW-FILTERED-DGE.rds")
scrna_pheno <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/TC1-PHENODATA.rds")
identical(row.names(scrna_pheno), names(scrna))



# Create artifical scRNA bulk samples
clusters_of_interest <- c(1, 23)
scrna_aggr <- make_scrna_aggregrate(scrna, labels = as.vector(scrna_pheno$Cluster), normalize = TRUE)
scrna_aggr <- scrna_aggr[row.names(encode_expr), paste("Cluster", clusters_of_interest, sep = "_")]



# How well do these data correlate (ENCODE vs aggregrated scRNA)?
corrs <- list(p_cor = cor(cbind.data.frame(scrna_aggr, encode_expr), method = "pearson")[,1:2],
              s_cor = cor(cbind.data.frame(scrna_aggr, encode_expr), method = "spearman")[,1:2])
lapply(seq_along(corrs), function(x){
  ph <- pheatmap::pheatmap(corrs[[x]],
                     display_numbers = TRUE,
                     cluster_cols = FALSE,
                     main = toupper(names(corrs)[x]),
                     filename = file.path("~/Dropbox/ErnstLab/deconR/example", paste0(toupper(names(corrs)[x]), "-COR-BULK-W-SCRNA.png")),
                     height = 14,
                     width = 5)
})








# NOTES
# save(encode_signal, file = file.path("~/Dropbox/ErnstLab/deconR/data", "encode_signal.RData"), compress='xz')
