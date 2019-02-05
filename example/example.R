# CLEAN ENV
rm(list = ls(all = TRUE))
options(warn = -1)


# INSTALL AND LOAD LIB
# if(!require(devtools)) {install.packages(devtools)}
# devtools::install_github("ShanSabri/deconR")
library(deconR)

# LOAD ENCODE DATA
data(encode_signal)
data(encode_expr)

head(encode_expr)




# encode_expr <- readRDS(paste("~/Dropbox/ErnstLab/decon_mouse", "data", "ENCODE-RNA-COUNTS-GENES-MEAN-REPLICATES-WITH-SCRNA-DEPTHNORM.rds", sep = "/"))
#
# a <- readRDS(paste("~/Dropbox/ErnstLab/decon_mouse/", "data", "ENCODE-CHIP-H3K27AC-CHR16-LOG2.rds", sep="/"))
# md <- a[,1:5]
# a <- a[,seq(6, ncol(a))]
# a <- 2^a - 1
# encode_signal <- a
#
# common_ct <- sort(intersect(names(encode_expr), names(mcols(encode_signal))))
# encode_expr <- encode_expr[, common_ct]
# encode_signal <- encode_signal[, common_ct]
# identical(names(encode_expr), names(encode_signal))
#
#
# library(GenomicRanges)
# encode_signal <- makeGRangesFromDataFrame(cbind.data.frame(md, encode_signal),
#                                           keep.extra.columns=TRUE,
#                                           ignore.strand=FALSE,
#                                           seqinfo=NULL,
#                                           seqnames.field=c("seqnames", "seqname",
#                                                            "chromosome", "chrom",
#                                                            "chr", "chromosome_name",
#                                                            "seqid"),
#                                           start.field="start",
#                                           end.field=c("end", "stop"),
#                                           strand.field="strand",
#                                           starts.in.df.are.0based=FALSE)
# class(encode_signal)
# save(encode_signal, file = file.path("~/Dropbox/ErnstLab/deconR/data", "encode_signal.RData"), compress='xz')
# save(encode_expr, file = file.path("~/Dropbox/ErnstLab/deconR/data", "encode_expr.RData"), compress='xz')
#
#
# names(encode_expr); length(names(encode_expr))
# names(mcols(encode_signal)); length(names(mcols(encode_signal)))
# identical(names(encode_expr), names(mcols(encode_signal)))
