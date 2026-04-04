args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5){
  stop("Usage:
  Rscript roh_froh_cli.R roh.hom pel.txt mad.txt bim out_prefix")
}

roh_file <- args[1]
pel_file <- args[2]
mad_file <- args[3]
bim_file <- args[4]
out_prefix <- args[5]

library(detectRUNS)
library(dplyr)
library(openxlsx)

# =========================
# ID
# =========================
id_pelaihari <- trimws(read.table(pel_file)$V1)
id_madura    <- trimws(read.table(mad_file)$V1)

meta <- data.frame(
  id = c(id_pelaihari, id_madura),
  population = c(
    rep("Pelaihari", length(id_pelaihari)),
    rep("Madura", length(id_madura))
  )
)

# =========================
# MAP FILE
# =========================
bim <- read.table(bim_file)
map_file <- paste0(out_prefix, "_map.txt")

write.table(
  bim[bim$V1 %in% 1:29, 1:4],
  map_file,
  sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE
)

# =========================
# ROH
# =========================
runs <- readExternalRuns(roh_file, program="plink")

froh <- Froh_inbreedingClass(runs, map_file, Class = 1)

final_froh <- meta %>%
  left_join(froh, by="id")

# =========================
# STATS
# =========================
runs$population <- ifelse(runs$id %in% id_madura, "Madura", "Pelaihari")
runs$length_Mb <- runs$lengthBps / 1e6

indiv_stats <- runs %>%
  group_by(id, population) %>%
  summarise(
    nROH = n(),
    Total_Mb = sum(length_Mb),
    Avg_Mb = mean(length_Mb),
    .groups="drop"
  )

journal_table <- indiv_stats %>%
  group_by(population) %>%
  summarise(
    n = n(),
    nROH = paste0(round(mean(nROH),2)," ± ",round(sd(nROH),2)),
    Total_Mb = paste0(round(mean(Total_Mb),2)," ± ",round(sd(Total_Mb),2)),
    Avg_Mb = paste0(round(mean(Avg_Mb),2)," ± ",round(sd(Avg_Mb),2)),
    .groups="drop"
  )

# =========================
# SAVE
# =========================
wb <- createWorkbook()

addWorksheet(wb, "FROH")
writeData(wb, "FROH", final_froh)

addWorksheet(wb, "ROH_stats")
writeData(wb, "ROH_stats", indiv_stats)

addWorksheet(wb, "Journal")
writeData(wb, "Journal", journal_table)

saveWorkbook(wb, paste0(out_prefix, "_ROH.xlsx"), overwrite = TRUE)