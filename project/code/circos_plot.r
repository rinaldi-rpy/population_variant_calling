# =========================================================
# MASTER VISUALIZATION (FINAL CLI VERSION)
# VENN → CIRCOS → VENN PERCENT
# TANPA PERUBAHAN LOGIKA
# =========================================================

library(data.table)
library(circlize)
library(ggplot2)
library(grid)
library(readr)
library(ggforce)

# =========================================================
# 0. ARGUMENT INPUT
# =========================================================

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 6){
  stop("Usage:
  Rscript script.R madura_snps pelaihari_snps madura_count pelaihari_count shared_count output_prefix")
}

madura_file      <- args[1]
pelaihari_file   <- args[2]
madura_count_file    <- args[3]
pelaihari_count_file <- args[4]
shared_count_file    <- args[5]
out_prefix       <- args[6]

# =========================================================
# 1. LOAD SNP COUNT (SOURCE OF TRUTH)
# =========================================================

madura_count    <- as.numeric(fread(madura_count_file)$V1)
pelaihari_count <- as.numeric(fread(pelaihari_count_file)$V1)
shared_count    <- as.numeric(fread(shared_count_file)$V1)

# =========================================================
# 2. BUILD VENN (COUNT BASED)
# =========================================================

circle_df <- function(center, r, n=200){
  theta <- seq(0, 2*pi, length.out = n)
  data.frame(
    x = center[1] + r*cos(theta),
    y = center[2] + r*sin(theta)
  )
}

r1 <- sqrt((madura_count + shared_count)/pi)
r2 <- sqrt((pelaihari_count + shared_count)/pi)

c1 <- circle_df(c(0,0), r1);        c1$group <- "Madura"
c2 <- circle_df(c(r1*0.8,0), r2);   c2$group <- "Pelaihari"

circles <- rbind(c1, c2)

labels <- data.frame(
  x = c(-r1*0.4, r1*0.4, r1*1.1),
  y = c(0,0,0),
  label = c(madura_count, shared_count, pelaihari_count)
)

ven_plot <- ggplot() +
  geom_polygon(data=circles,
               aes(x=x, y=y, group=group, fill=group),
               alpha=0.5, color="black") +
  geom_text(data=labels,
            aes(x=x, y=y, label=label),
            size=6, fontface="bold") +
  scale_fill_manual(values=c("Madura"="red","Pelaihari"="blue")) +
  theme_void() +
  coord_equal()

# =========================================================
# 3. LOAD AF DATA (CIRCOS)
# =========================================================

madura <- fread(madura_file)
pelaihari <- fread(pelaihari_file)

setnames(madura,   c("SNP","AF_pel","AF_mad"))
setnames(pelaihari,c("SNP","AF_pel","AF_mad"))

madura[,   c("chr","pos") := tstrsplit(SNP,":")]
pelaihari[,c("chr","pos") := tstrsplit(SNP,":")]

madura[,   pos := as.numeric(pos)]
pelaihari[,pos := as.numeric(pos)]

madura_gen <- madura[,.(chr=as.numeric(chr), start=pos, end=pos, AF=AF_mad)]
pelaihari_gen <- pelaihari[,.(chr=as.numeric(chr), start=pos, end=pos, AF=AF_pel)]

chr_levels <- as.character(1:29)

madura_gen[, chr := factor(chr, levels = chr_levels)]
pelaihari_gen[, chr := factor(chr, levels = chr_levels)]

all_snps <- rbind(
  madura_gen[,.(chr,start)],
  pelaihari_gen[,.(chr,start)]
)

chr_size <- all_snps[,.(max_pos = max(start)), by=chr]
chr_size[, chr := factor(chr, levels = chr_levels)]

# =========================================================
# 4. DRAW CIRCOS
# =========================================================

png(paste0(out_prefix, "_circos.png"), width=3000, height=3000, res=300)

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree   = 2,
  track.margin = c(0.01,0.01)
)

circos.initialize(
  factors = chr_size$chr,
  xlim    = cbind(rep(0,nrow(chr_size)), chr_size$max_pos)
)

# OUTER CHR
circos.trackPlotRegion(
  ylim = c(0,1),
  track.height = 0.05,
  bg.border = "black",
  panel.fun = function(x,y){
    chr <- CELL_META$sector.index
    circos.text(
      mean(CELL_META$xlim),
      0.5,
      paste0("Chr", chr),
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.6
    )
  }
)

# MADURA
circos.genomicTrackPlotRegion(
  madura_gen,
  ylim = c(0,1),
  track.height = 0.12,
  panel.fun = function(region, value, ...){
    circos.genomicPoints(region, value, col="red", pch=16, cex=0.4)
  }
)

# PELAIHARI
circos.genomicTrackPlotRegion(
  pelaihari_gen,
  ylim = c(0,1),
  track.height = 0.12,
  panel.fun = function(region, value, ...){
    circos.genomicPoints(region, value, col="blue", pch=16, cex=0.4)
  }
)

# INSERT VENN
grid.draw(ggplotGrob(ven_plot))

dev.off()

# =========================================================
# 5. VENN PERCENT (SEPARATE FIGURE)
# =========================================================

madura_only    <- as.numeric(readLines(madura_count_file))
pelaihari_only <- as.numeric(readLines(pelaihari_count_file))
shared         <- as.numeric(readLines(shared_count_file))

total <- madura_only + pelaihari_only + shared

pct_madura    <- round(madura_only / total * 100, 2)
pct_pelaihari <- round(pelaihari_only / total * 100, 2)
pct_shared    <- round(shared / total * 100, 2)

circle_df2 <- data.frame(
  x0 = c(-1, 1),
  y0 = c(0, 0),
  r  = c(2, 2),
  group = c("Madura", "Pelaihari")
)

vendiagram <- ggplot() +
  geom_circle(data = circle_df2,
              aes(x0 = x0, y0 = y0, r = r, fill = group),
              alpha = 0.4,
              color = "black") +
  annotate("text", x = -1.5, y = 0,
           label = paste0(pct_madura, "%"), size = 8) +
  annotate("text", x = 0, y = 0,
           label = paste0(pct_shared, "%"), size = 8) +
  annotate("text", x = 1.5, y = 0,
           label = paste0(pct_pelaihari, "%"), size = 8) +
  coord_fixed() +
  scale_fill_manual(values = c(
    "Madura" = "lightgreen",
    "Pelaihari" = "lightpink"
  )) +
  theme_void()

ggsave(
  paste0(out_prefix, "_venn_percentage.png"),
  vendiagram,
  width = 6,
  height = 6,
  dpi = 600
)