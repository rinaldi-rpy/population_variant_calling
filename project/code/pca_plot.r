#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--meta", type="character"),
  make_option("--eigenvec", type="character"),
  make_option("--eigenval", type="character"),
  make_option("--out", type="character", default="PCA.png")
)

opt <- parse_args(OptionParser(option_list=option_list))

meta_path     <- opt$meta
eigenvec_path <- opt$eigenvec
eigenval_path <- opt$eigenval
out_path      <- opt$out

meta <- read.csv(meta_path, stringsAsFactors = FALSE)
stopifnot(all(c("ID","Breed","Bangsa") %in% names(meta)))

meta$ID     <- trimws(meta$ID)
meta$Breed  <- trimws(meta$Breed)
meta$Bangsa <- trimws(meta$Bangsa)

eig <- read.table(eigenvec_path, header = FALSE, stringsAsFactors = FALSE)
n_pc <- ncol(eig) - 2
stopifnot(n_pc >= 2)

colnames(eig) <- c("FID","IID",paste0("PC", seq_len(n_pc)))

eigval <- scan(eigenval_path, quiet = TRUE)
pct <- eigval / sum(eigval) * 100

xlab <- sprintf("PC1 (%.2f%%)", pct[1])
ylab <- sprintf("PC2 (%.2f%%)", pct[2])

dat <- merge(eig, meta, by.x="IID", by.y="ID", all.x=TRUE)

cat("\n=== Missing metadata ===\n")
print(dat[is.na(dat$Breed) | is.na(dat$Bangsa), c("IID","PC1","PC2")])

cat("\n=== Missing PC ===\n")
print(dat[is.na(dat$PC1) | is.na(dat$PC2), c("IID","Breed","Bangsa")])

dat$Bangsa <- factor(dat$Bangsa)
dat$Breed  <- factor(dat$Breed)

########################################
# AUTO SHAPE & COLOR (GENERALIZED)
########################################

n_group <- length(levels(dat$Bangsa))

shape_values <- rep(c(16,17,15,3,0,8,1,2,7,9), length.out = n_group)

########################################
# PLOT
########################################

p <- ggplot(dat, aes(PC1, PC2)) +
  geom_point(aes(color = Bangsa, shape = Bangsa),
             size = 3, alpha = 0.9, na.rm = TRUE) +
  scale_shape_manual(values = shape_values) +
  labs(x = xlab, y = ylab, color = "Bangsa", shape = "Bangsa") +
  theme_bw()

print(p)

ggsave(out_path, plot = p, width = 8, height = 5, dpi = 600)
