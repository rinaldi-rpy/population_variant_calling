args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5){
  stop("Usage:
  Rscript roh_island_cli.R roh.hom pel.txt mad.txt gtf out_prefix")
}

roh_file <- args[1]
pel_file <- args[2]
mad_file <- args[3]
gtf_file <- args[4]
out_prefix <- args[5]

library(detectRUNS)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(rtracklayer)

# =========================
# ID
# =========================
id_pelaihari <- trimws(as.character(read.table(pel_file)$V1))
id_madura    <- trimws(as.character(read.table(mad_file)$V1))

runs <- readExternalRuns(roh_file, program="plink")

runs$population <- ifelse(runs$id %in% id_madura, "Madura", "Pelaihari")

# =========================
# GRanges
# =========================
roh_gr <- GRanges(
  seqnames = runs$chr,
  ranges   = IRanges(start = runs$from, end = runs$to),
  population = runs$population
)

chr_lengths <- tapply(runs$to, runs$chr, max)

windows <- tileGenome(
  seqlengths = chr_lengths,
  tilewidth = 40000,
  cut.last.tile.in.chrom = TRUE
)

roh_island <- lapply(unique(runs$population), function(pop){

  roh_pop <- roh_gr[roh_gr$population == pop]
  n_ind   <- length(unique(roh_pop$id))

  hits <- findOverlaps(windows, roh_pop)

  df <- data.frame(
    window = queryHits(hits),
    id     = roh_pop$id[subjectHits(hits)]
  )

  df %>%
    distinct(window, id) %>%
    count(window) %>%
    mutate(freq = n / n_ind, population = pop)

}) %>% bind_rows()

roh_island <- roh_island %>% filter(freq >= 0.30)

# =========================
# SAVE
# =========================
write.table(
  roh_island,
  paste0(out_prefix, "_roh_island.txt"),
  sep="\t", row.names=FALSE, quote=FALSE
)