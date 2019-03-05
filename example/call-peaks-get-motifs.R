# Clean environment
rm(list = ls(all = TRUE))
options(warn = -1)

library(GenomicRanges)
library(tidyverse)


# CALL PEAKS
input_files <- list.files(path = file.path(getwd(), "example", "tracks"), pattern = "\\.bw", ignore.case = TRUE, full.names = TRUE)
atlas <- do.call("c", lapply(seq_along(input_files), function(x){
  bedgraph <- file.path(dirname(input_files[[x]]), "peaks", "tmp.bedGraph")
  output_file <- gsub(".bw", ".bed", basename(input_files[[x]]))
  rtracklayer::export(rtracklayer::import(input_files[[x]]), format="bedGraph", con=bedgraph)
  cmd <- paste("macs2 bdgbroadcall",
               "-i", bedgraph,
               "--outdir", dirname(bedgraph),
               "-o", output_file
  )
  message(cmd)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE); system(paste("rm", bedgraph))
  peaks <- rtracklayer::import(file.path(dirname(bedgraph), output_file))
  mcols(peaks) <- data.frame(ID = gsub(".bw", "", basename(input_files[[x]])))
  return(peaks)
}))
saveRDS(atlas, compress = TRUE, file.path(getwd(), "example", "tracks", "peaks", "PEAK-ATLAS-ACROSS-MODEL.rds"))




# FILTER PEAKS
atlas <- readRDS(file.path(getwd(), "example", "tracks", "peaks", "PEAK-ATLAS-ACROSS-MODEL.rds"))
atlas_split <- split(atlas, f = mcols(atlas))


# # Method 1: UNIQ
# f <- 'MEF-ACROSS-MODEL-LOG2-UNIQ.bed'
# a <- atlas_split$`MEF-LOG2`
# b <- atlas_split$`ESC-LOG2`
# a <- data.frame(chr = seqnames(a), start = start(a), end = end(a)); a$chr <- as.character(a$chr)
# b <- data.frame(chr = seqnames(b), start = start(b), end = end(b)); b$chr <- as.character(b$chr)
# only_a <- bedr::bedr.subtract.region(a, b)
# write.table(only_a, file=file.path(getwd(), "example", "tracks", "peaks", f), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Method 2: by max/average signal intensity within a block
atlas_split_filtered <- lapply(seq_along(atlas_split), function(x){
  id <- names(atlas_split)[x]
  peaks <- atlas_split[[x]]
  signal <- rtracklayer::import(file.path(getwd(), "example", "tracks", paste0(id, ".bw")))
  message(id)

  hits <- findOverlaps(signal, peaks, ignore.strand=FALSE)

  signal_in_peaks <- tibble::as.tibble(data.frame(
    peak_id    = hits@to,
    peak_start = start(peaks)[hits@to],
    peak_end   = end(peaks)[hits@to],
    signal     = score(signal)[hits@from]
  )) %>%
    group_by(peak_id) %>%
    summarise(avg_signal = mean(signal),
              max_signal = max(signal),
              peak_start = unique(peak_start),
              peak_end = unique(peak_end))

  signal_in_peaks %>%
    reshape2::melt(., id.vars = c("peak_id", "peak_start", "peak_end")) %>%
    ggplot(aes(value, fill = variable)) +
      facet_wrap(~variable, nrow = 1, scales = "free") +
      geom_density(aes(y = ..scaled..), alpha = 0.5) +
      labs(x = "Signal Intensity", y = "Density Distribution", title = id) +
      theme_bw(base_size = 14) + theme(panel.grid = element_blank(), legend.position = "none", axis.text = element_text(size = 10)) +
      ggsave(file.path(getwd(), "example", "tracks", "peaks", "signal-in-peak-dist", paste0(id, "-SIGNAL-IN-PEAK-DIST.png")), height = 5, width = 9)

  threshold <- quantile(signal_in_peaks$avg_signal, 0.50)
  filtered_peaks <- subset(signal_in_peaks, avg_signal >= threshold)

  filtered_peaks_gr <- peaks[start(peaks) %in% filtered_peaks$peak_start]
  return(filtered_peaks_gr)
})
names(atlas_split_filtered) <- names(atlas_split)
saveRDS(atlas_split_filtered, compress = TRUE, file.path(getwd(), "example", "tracks", "peaks", "filtered-peaks", "PEAK-ATLAS-ACROSS-MODEL.rds"))

tibble::as.tibble(data.frame(
  id = names(atlas_split),
  before = sapply(atlas_split, length),
  after = sapply(atlas_split_filtered, length),
  diff = sapply(atlas_split, length) - sapply(atlas_split_filtered, length)))

out <- lapply(seq_along(atlas_split_filtered), function(x){
  id <- names(atlas_split_filtered)[x]; message(id)
  rtracklayer::export(atlas_split_filtered[[x]],
                      file.path(getwd(), "example", "tracks", "peaks", "filtered-peaks",  paste0(id, ".bed")), format = "bed")
})



# MEF & ESC SPECIFIC PEAK - REMOVE OVERLAP
f <- 'ESC-ACROSS-MODEL-LOG2-UNIQ.bed'
a <- atlas_split_filtered$`ESC-LOG2`
b <- atlas_split_filtered$`MEF-LOG2`
a <- data.frame(chr = seqnames(a), start = start(a), end = end(a)); a$chr <- as.character(a$chr)
b <- data.frame(chr = seqnames(b), start = start(b), end = end(b)); b$chr <- as.character(b$chr)
only_a <- bedr::bedr.subtract.region(a, b)
write.table(only_a, file=file.path(getwd(), "example", "tracks", "peaks", "filtered-peaks", f), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



# FIND MOTIFS
run_homer <- function(file, fa, out, dry = TRUE) {
  Sys.setenv(PATH=paste0("/usr/local/opt/libxml2/bin:/Users/shansabri/miniconda3/bin:/Library/Frameworks/Python.framework/Versions/3.4/bin:",
                         "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/ncbi/blast/bin:/Library/TeX/texbin:/usr/local/MacGPG2/bin:",
                         "/opt/X11/bin:/Users/shansabri/bin:/Users/shansabri/edirect:/Users/shansabri/homer/bin/:/Users/shansabri/bin:/Users/shansabri/edirect"))
  cmd <- sprintf("findMotifsGenome.pl %s %s %s -p 7 -size given", file, fa, out)
  if(dry == TRUE) {
    return(cmd)
  } else {
    try(system(cmd))
    return(TRUE)
  }
}
lapply(c("MEF-ACROSS-MODEL-LOG2-UNIQ.bed", "ESC-ACROSS-MODEL-LOG2-UNIQ.bed",
         "CLUSTER_1-ACROSS-MODEL-LOG2-UNIQ.bed", "CLUSTER_23-ACROSS-MODEL-LOG2-UNIQ.bed",
         "BASELINE-MODEL-LOG2.bed", "MEF-LOG2.bed", "ESC-LOG2.bed"), function(f){
  run_homer(file = file.path(getwd(), "example", "tracks", "peaks", "filtered-peaks", f),
            fa = "~/Dropbox/Binhouse/mm9.fa",
            out = file.path(getwd(),  "example", "motifs", "from-filtered-peaks", tools::file_path_sans_ext(f)),
            dry = FALSE)
})




# PLOT MOTIF ENRICHMENTS
motif_dir <- file.path(getwd(),  "example", "motifs")
res <- list.dirs(motif_dir, full.names = TRUE, recursive = FALSE)
all_motif_results <- do.call(rbind.data.frame, lapply(seq_along(res), function(x) {
  # x = 1
  message(paste(x, basename(res[[x]])))
  # denovo <- read_denovo_results(path = x, homer_dir = TRUE)
  known <- marge::read_known_results(path = res[[x]], homer_dir = TRUE)
  known$ID <- basename(res)[x]
  known <- subset(known, database == "Homer")
  return(known)
  # }
}))
saveRDS(all_motif_results, compress = TRUE, file.path(motif_dir, "RESULTS-AGGREGRATED.rds"))




# WHAT ARE THE TOP 10 MOTIFS FOR EACH GENE SET?
all_motif_results %>%
  group_by(ID) %>%
  top_n(log_p_value, n = 10) %>%
  dplyr::select(., motif_name, motif_family, experiment, accession, log_p_value) %>%
  ungroup() %>%
  write_tsv(file.path(motif_dir, "TOP10-MOTIFS-BY-ID.txt"),
            na = "NA",  quote_escape = "double")


motif_df <- reshape2::dcast(all_motif_results, motif_name + motif_family + accession + consensus ~ ID,
                            value.var = "log_p_value", fun.aggregate = max
)
row.names(motif_df) <- paste(motif_df$motif_name, motif_df$motif_family, motif_df$accession, motif_df$consensus, sep = "/")
motif_df$motif_name <- NULL
motif_df$motif_family <- NULL
motif_df$accession <- NULL
motif_df$consensus <- NULL
write.table(motif_df, file.path(motif_dir, "MOTIF-BY-PROG.txt"), col.names = NA, quote = FALSE, sep = "\t")

motif_df$max_p_val <- apply(motif_df, 1, max)
ggplot(motif_df, aes(max_p_val)) + geom_density(fill = 'yellow'); summary(motif_df$max_p_val)
maxs <- seq(1, 10, by = 1)
pdf(paste(motif_dir, paste0("ENRICHED_MOTIFS_MAX_THRESHOLD.pdf"), sep = "/"), height = 15, width = 10)
lapply(maxs, function(x) {
  message(x)
  tmp <- subset(motif_df, max_p_val >= x)
  tmp$max_p_val <- NULL
  gplots::heatmap.2(as.matrix(tmp),
                    Rowv = TRUE,
                    Colv = TRUE,
                    trace = "none",
                    na.color = "grey",
                    margins = c(35, 20),
                    col = cm.colors(255),
                    main = paste("Max log_p_value >=", x, "\n", nrow(tmp), "motifs")
  )
})
dev.off(); graphics.off()

pdf(paste(motif_dir, paste0("ENRICHED_MOTIFS_MAX_THRESHOLD-LOG10.pdf"), sep = "/"), height = 15, width = 20)
lapply(maxs, function(x) {
  message(x)
  tmp <- subset(motif_df, max_p_val >= x)
  tmp$max_p_val <- NULL
  h <- as.matrix(log10(tmp + 1))
  h[is.infinite(h)] <- 0
  gplots::heatmap.2(h,
                    Rowv = TRUE,
                    Colv = TRUE,
                    trace = "none",
                    na.color = "grey",
                    margins = c(35, 20),
                    col = cm.colors(255),
                    main = paste("Max log_p_value >=", x, "\n", nrow(tmp), "motifs")
  )
})
dev.off(); graphics.off()
