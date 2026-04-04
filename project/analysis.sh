#!/bin/bash
set -euo pipefail

############################################
# INPUT
############################################

VCF_MADURA="${1:?}"
VCF_BOSJAVA="${2:?}"
VCF_INDIA="${3:?}"
TXT_PEL="${4:?}"
TXT_MAD="${5:?}"
META="${6:?}"
GTF="${7:?}"
MAF="${8:?}"
MISS="${9:?}"

############################################
# DIR
############################################

BASE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

DATA="$BASE_DIR/data"
INPUT="$BASE_DIR/input"
WORK="$BASE_DIR/work"
RES="$BASE_DIR/result"
FIG="$RES/figure"
TABLE="$RES/tables"
CODE="$BASE_DIR/code"

mkdir -p "$WORK"/{01_qc,02_filter,03_merge,04_merge_filter,05_pca,05_pca_madura,06_subset,07_isec,07_isec_stats,09_fst,10_roh,12_admixture,14_af_divergence}
mkdir -p "$FIG" "$TABLE"

############################################
# 1. QC
############################################

for POP in madura bosjava india; do
  VCF_VAR="VCF_${POP^^}"
  vcffilter -f 'QUAL > 20 & DP > 4' "${!VCF_VAR}" \
    | bgzip -c > "$WORK/01_qc/${POP}.qc.vcf.gz"
  tabix -p vcf "$WORK/01_qc/${POP}.qc.vcf.gz"
done

############################################
# 2. FILTER PER POPULASI
############################################

for POP in madura bosjava india; do
  vcftools --gzvcf "$WORK/01_qc/${POP}.qc.vcf.gz" \
    --maf "$MAF" \
    --max-missing "$MISS" \
    --recode --recode-INFO-all \
    --out "$WORK/02_filter/${POP}"

  bgzip -f "$WORK/02_filter/${POP}.recode.vcf"
  tabix -p vcf "$WORK/02_filter/${POP}.recode.vcf.gz"
done

############################################
# 3. MERGE VCF
############################################

bcftools merge \
  "$WORK/02_filter/bosjava.recode.vcf.gz" \
  "$WORK/02_filter/india.recode.vcf.gz" \
  "$WORK/02_filter/madura.recode.vcf.gz" \
  -Oz -o "$WORK/03_merge/merged.vcf.gz"

tabix -p vcf "$WORK/03_merge/merged.vcf.gz"

############################################
# 4. FILTER MERGED (UNTUK PCA GLOBAL)
############################################

vcftools --gzvcf "$WORK/03_merge/merged.vcf.gz" \
  --maf "$MAF" \
  --max-missing "$MISS" \
  --recode --recode-INFO-all \
  --out "$WORK/04_merge_filter/merged.filtered"

bgzip -f "$WORK/04_merge_filter/merged.filtered.recode.vcf"
tabix -p vcf "$WORK/04_merge_filter/merged.filtered.recode.vcf.gz"

############################################
# 5. PCA GLOBAL
############################################

plink --vcf "$WORK/04_merge_filter/merged.filtered.recode.vcf.gz" \
  --allow-extra-chr --chr-set 29 \
  --snps-only --double-id \
  --make-bed --out "$WORK/05_pca/all"

plink --bfile "$WORK/05_pca/all" \
  --indep-pairwise 50 10 0.2 \
  --out "$WORK/05_pca/all"

plink --bfile "$WORK/05_pca/all" \
  --extract "$WORK/05_pca/all.prune.in" \
  --make-bed \
  --out "$WORK/05_pca/all.pruned"

plink --bfile "$WORK/05_pca/all.pruned" \
  --pca 4 \
  --out "$WORK/05_pca/all.pruned"

############################################
# 6. ADMIXTURE
############################################

plink --bfile "$WORK/05_pca/all.pruned" \
  --geno 0.95 \
  --make-bed \
  --chr-set 29 \
  --out "$WORK/12_admixture/all.clean"

for K in 2 3 4; do
  admixture --cv "$WORK/12_admixture/all.clean.bed" $K \
    | tee "$WORK/12_admixture/log${K}.out"
done

grep -h CV "$WORK"/12_admixture/log*.out \
  > "$WORK/12_admixture/cv_error.txt"

############################################
# 7. PCA MADURA
############################################

plink --vcf "$WORK/02_filter/madura.recode.vcf.gz" \
  --allow-extra-chr --chr-set 29 \
  --snps-only --double-id \
  --make-bed \
  --out "$WORK/05_pca_madura/madura"

plink --bfile "$WORK/05_pca_madura/madura" \
  --indep-pairwise 50 10 0.2 \
  --out "$WORK/05_pca_madura/madura"

plink --bfile "$WORK/05_pca_madura/madura" \
  --extract "$WORK/05_pca_madura/madura.prune.in" \
  --make-bed \
  --out "$WORK/05_pca_madura/madura.pruned"

plink --bfile "$WORK/05_pca_madura/madura.pruned" \
  --pca 4 \
  --out "$WORK/05_pca_madura/madura.pruned"

############################################
# 8. SUBSET MADURA
############################################

bcftools view -S "$TXT_PEL" "$WORK/02_filter/madura.recode.vcf.gz" \
  -Oz -o "$WORK/06_subset/pel.vcf.gz"
tabix -p vcf "$WORK/06_subset/pel.vcf.gz"

bcftools view -S "$TXT_MAD" "$WORK/02_filter/madura.recode.vcf.gz" \
  -Oz -o "$WORK/06_subset/mad.vcf.gz"
tabix -p vcf "$WORK/06_subset/mad.vcf.gz"

############################################
# 9. ISEC (VENN)
############################################

bcftools isec -p "$WORK/07_isec" \
  "$WORK/06_subset/pel.vcf.gz" \
  "$WORK/06_subset/mad.vcf.gz"

for i in 0000 0001 0002; do
  bcftools view -H "$WORK/07_isec/$i.vcf" | wc -l \
    > "$WORK/07_isec_stats/${i}_snp_count.txt"
done

############################################
# 10. FST
############################################

vcftools \
  --gzvcf "$WORK/02_filter/madura.recode.vcf.gz" \
  --weir-fst-pop "$TXT_PEL" \
  --weir-fst-pop "$TXT_MAD" \
  --fst-window-size 40000 \
  --fst-window-step 20000 \
  --out "$WORK/09_fst/fst"

############################################
# 11. ROH (PLINK)
############################################

plink --bfile "$WORK/05_pca_madura/madura.pruned" \
  --chr-set 29 \
  --homozyg \
  --out "$WORK/10_roh/roh_madura"

############################################
# 12. AF DIVERGENCE
############################################

vcftools --gzvcf "$WORK/06_subset/pel.vcf.gz" --freq --out "$WORK/14_af_divergence/pel"
vcftools --gzvcf "$WORK/06_subset/mad.vcf.gz" --freq --out "$WORK/14_af_divergence/mad"

awk 'NR>1{match($5,/:[0-9.]+/);af=substr($5,RSTART+1,RLENGTH-1)+0;print $1":"$2,af}' \
"$WORK/14_af_divergence/pel.frq" | sort > "$WORK/14_af_divergence/pel.txt"

awk 'NR>1{match($5,/:[0-9.]+/);af=substr($5,RSTART+1,RLENGTH-1)+0;print $1":"$2,af}' \
"$WORK/14_af_divergence/mad.frq" | sort > "$WORK/14_af_divergence/mad.txt"

join "$WORK/14_af_divergence/pel.txt" "$WORK/14_af_divergence/mad.txt" \
> "$WORK/14_af_divergence/merged.txt"

awk '$2>0.6 && $3<0.3' "$WORK/14_af_divergence/merged.txt" \
> "$WORK/14_af_divergence/pel_specific.txt"

awk '$3>0.6 && $2<0.3' "$WORK/14_af_divergence/merged.txt" \
> "$WORK/14_af_divergence/mad_specific.txt"

############################################
# 13. R VISUALIZATION
############################################

Rscript "$CODE/pca_plot.r" \
  --meta "$META" \
  --eigenvec "$WORK/05_pca/all.pruned.eigenvec" \
  --eigenval "$WORK/05_pca/all.pruned.eigenval" \
  --out "$FIG/PCA_all.png"

Rscript "$CODE/pca_plot.r" \
  --meta "$META" \
  --eigenvec "$WORK/05_pca_madura/madura.pruned.eigenvec" \
  --eigenval "$WORK/05_pca_madura/madura.pruned.eigenval" \
  --out "$FIG/PCA_madura.png"

Rscript "$CODE/admixture_plot.r" \
  --fam "$WORK/12_admixture/all.clean.fam" \
  --meta "$META" \
  --cv "$WORK/12_admixture/cv_error.txt" \
  --prefix "$WORK/12_admixture/all.clean" \
  --k 2,3,4 \
  --out "$FIG/admixture.png"

Rscript "$CODE/circos_plot.r" \
  --madura "$WORK/14_af_divergence/mad_specific.txt" \
  --pelaihari "$WORK/14_af_divergence/pel_specific.txt" \
  --fst "$WORK/09_fst/fst.windowed.weir.fst" \
  --out "$FIG/circos.png"

Rscript "$CODE/roh_froh_cli.R" \
  "$WORK/10_roh/roh_madura.hom" \
  "$TXT_PEL" \
  "$TXT_MAD" \
  "$WORK/05_pca_madura/madura.pruned.bim" \
  "$TABLE/roh"

Rscript "$CODE/roh_island_cli.R" \
  "$WORK/10_roh/roh_madura.hom" \
  "$TXT_PEL" \
  "$TXT_MAD" \
  "$GTF" \
  "$FIG/roh"

Rscript "$CODE/candidate_gene.r" \
  --madura "$WORK/14_af_divergence/mad_specific.txt" \
  --pelaihari "$WORK/14_af_divergence/pel_specific.txt" \
  --fst "$WORK/09_fst/fst.windowed.weir.fst" \
  --roh "$WORK/10_roh/roh_madura.hom" \
  --gtf "$GTF" \
  --outdir "$TABLE"

############################################
echo "PIPELINE SELESAI"