#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

########################################
# ARGUMENT
########################################

option_list <- list(
  make_option("--madura", type="character"),
  make_option("--pelaihari", type="character"),
  make_option("--fst", type="character"),
  make_option("--roh", type="character"),
  make_option("--gtf", type="character"),
  make_option("--outdir", type="character", default="result")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, showWarnings=FALSE)

########################################
# FUNCTION
########################################

norm_chr <- function(x){
  x <- gsub("^chr","",x,ignore.case=TRUE)
  toupper(x)
}

allowed_chr <- as.character(1:29)

########################################
# LOAD AF
########################################

madura <- fread(opt$madura)
pelaihari <- fread(opt$pelaihari)

setnames(madura,c("SNP","AF_pel","AF_mad"))
setnames(pelaihari,c("SNP","AF_pel","AF_mad"))

madura[,c("CHROM","POS") := tstrsplit(SNP,":")]
pelaihari[,c("CHROM","POS") := tstrsplit(SNP,":")]

madura$POS <- as.numeric(madura$POS)
pelaihari$POS <- as.numeric(pelaihari$POS)

madura$CHROM <- norm_chr(madura$CHROM)
pelaihari$CHROM <- norm_chr(pelaihari$CHROM)

madura <- madura[CHROM %in% allowed_chr]
pelaihari <- pelaihari[CHROM %in% allowed_chr]

########################################
# DELTA AF
########################################

madura[, delta_AF := abs(AF_mad - AF_pel)]
pelaihari[, delta_AF := abs(AF_mad - AF_pel)]

########################################
# LOAD FST
########################################

fst <- fread(opt$fst)
fst <- fst[,.(CHROM,BIN_START,BIN_END,WEIGHTED_FST)]
setnames(fst,c("CHROM","START","END","FST"))

fst$CHROM <- norm_chr(fst$CHROM)
fst <- fst[CHROM %in% allowed_chr]

setkey(fst,CHROM,START,END)

########################################
# LOAD ROH
########################################

roh <- fread(opt$roh)

setnames(roh,
c("population","chr","start","end","freq"),
c("POP","CHROM","START","END","ROH_FREQ"))

roh$CHROM <- norm_chr(roh$CHROM)
roh <- roh[CHROM %in% allowed_chr]

roh_mad <- roh[POP=="Madura"]
roh_pel <- roh[POP=="Pelaihari"]

setkey(roh_mad,CHROM,START,END)
setkey(roh_pel,CHROM,START,END)

########################################
# SNP INTERVAL
########################################

madura[,`:=`(SNP_START=POS,SNP_END=POS)]
pelaihari[,`:=`(SNP_START=POS,SNP_END=POS)]

setkey(madura,CHROM,SNP_START,SNP_END)
setkey(pelaihari,CHROM,SNP_START,SNP_END)

########################################
# OVERLAP
########################################

mad <- foverlaps(madura,fst,type="within",nomatch=NA)
mad <- foverlaps(mad,roh_mad,type="within",nomatch=NA)

pel <- foverlaps(pelaihari,fst,type="within",nomatch=NA)
pel <- foverlaps(pel,roh_pel,type="within",nomatch=NA)

########################################
# THRESHOLD
########################################

fst_thr <- quantile(fst$FST,0.99,na.rm=TRUE)
mad_thr <- quantile(mad$delta_AF,0.99,na.rm=TRUE)
pel_thr <- quantile(pel$delta_AF,0.99,na.rm=TRUE)

mad <- mad[FST>=fst_thr & delta_AF>=mad_thr]
pel <- pel[FST>=fst_thr & delta_AF>=pel_thr]

########################################
# LOAD GTF
########################################

gtf <- fread(opt$gtf,sep="\t",header=FALSE,skip="#")

setnames(gtf,
c("CHROM","SOURCE","FEATURE","START","END","SCORE","STRAND","FRAME","ATTRIBUTE"))

gtf$CHROM <- norm_chr(gtf$CHROM)

genes <- gtf[FEATURE=="gene"]

genes[, gene_name := sub('.*gene_name "([^"]+)".*','\\1',ATTRIBUTE)]
genes <- genes[,.(CHROM,START,END,gene_name)]

setkey(genes,CHROM,START,END)

########################################
# ANNOTATION
########################################

mad <- foverlaps(mad,genes,type="within",nomatch=NA)
pel <- foverlaps(pel,genes,type="within",nomatch=NA)

########################################
# CLEAN
########################################

mad <- mad[!is.na(gene_name) & !grepl("ENSBTAG",gene_name)]
pel <- pel[!is.na(gene_name) & !grepl("ENSBTAG",gene_name)]

########################################
# BEST SNP PER GENE
########################################

mad <- mad[, .SD[which.max(FST)], by=gene_name]
pel <- pel[, .SD[which.max(FST)], by=gene_name]

mad <- mad[order(-FST,-delta_AF)]
pel <- pel[order(-FST,-delta_AF)]

########################################
# SAVE
########################################

fwrite(mad[,.(CHROM,POS,FST,delta_AF,gene_name)],
       file.path(opt$outdir,"Madura_candidate_genes.csv"))

fwrite(pel[,.(CHROM,POS,FST,delta_AF,gene_name)],
       file.path(opt$outdir,"Pelaihari_candidate_genes.csv"))