#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

########################################
# ARGUMENT
########################################

option_list <- list(
  make_option("--fam", type="character"),
  make_option("--meta", type="character"),
  make_option("--cv", type="character"),
  make_option("--prefix", type="character"),
  make_option("--k", type="character"),
  make_option("--out", type="character", default="admixture.png")
)

opt <- parse_args(OptionParser(option_list=option_list))

########################################
# LOAD DATA
########################################

fam <- fread(opt$fam, header=FALSE)
colnames(fam) <- c("FID","ID","PID","MID","SEX","PHENO")

meta <- fread(opt$meta)

meta2 <- meta[match(fam$ID, meta$ID)]

data <- cbind(
  fam,
  meta2[, setdiff(names(meta2),"ID"), with=FALSE]
)

########################################
# POPULATION STRUCTURE
########################################

pop_info <- rle(data$Bangsa)

pop_counts <- pop_info$lengths
pop_names  <- pop_info$values

pop_breaks <- cumsum(pop_counts)
pop_centers <- pop_breaks - pop_counts/2

########################################
# CV ERROR
########################################

cv <- readLines(opt$cv)

cv <- gsub("CV error \\(K=", "", cv)
cv <- gsub("\\): ", " ", cv)

cv <- fread(text=cv)
colnames(cv) <- c("K","CV")

########################################
# K LIST (GENERAL)
########################################

Klist <- as.numeric(unlist(strsplit(opt$k, ",")))

########################################
# COLOR AUTO
########################################

get_colors <- function(K){
  base <- c("cyan","red","gold","blue","green","purple","orange","brown")
  rep(base, length.out=K)
}

########################################
# PLOT
########################################

png(opt$out, width=10, height=6, units="in", res=600)

par(
  mfrow=c(length(Klist),1),
  mar=c(1,4,1,6),
  oma=c(3,0,0,0)
)

for(K in Klist){

  Q <- fread(paste0(opt$prefix,".",K,".Q"))
  Q <- as.matrix(Q)

  cols <- get_colors(K)

  barplot(
    t(Q),
    border=NA,
    space=0,
    col=cols,
    axes=FALSE
  )

  axis(2, las=1)

  abline(v = pop_breaks, lwd=3)

  cv_val <- cv$CV[cv$K==K]

  mtext(
    paste0("K = ",K,"\nCV = ",round(cv_val,4)),
    side=4,
    line=2
  )
}

axis(
  1,
  at = pop_centers,
  labels = pop_names,
  tick = FALSE,
  line = 1
)

dev.off()