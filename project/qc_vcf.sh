#!/bin/bash
set -euo pipefail

############################################
# FUNCTION: QC VCF
############################################

qc_vcf() {
  local VCF="$1"
  local OUT="$2"

  mkdir -p "$OUT"

  echo "[INFO] QC dimulai"
  echo "[INFO] VCF: $VCF"
  echo "[INFO] Output: $OUT"

  ##########################################
  # 1. JUMLAH SNP
  ##########################################
  echo "[STEP 1] Jumlah SNP"
  bcftools view -H "$VCF" | wc -l > "$OUT/snp_count.txt"

  ##########################################
  # 2. DEPTH
  ##########################################
  echo "[STEP 2] Depth per site"
  vcftools --gzvcf "$VCF" \
    --site-mean-depth \
    --out "$OUT/depth"

  awk 'NR>1 {sum+=$3; n++} END {if(n>0) print sum/n; else print "NA"}' \
    "$OUT/depth.ldepth.mean" > "$OUT/mean_depth.txt"

  awk 'NR>1 {print $3}' "$OUT/depth.ldepth.mean" | sort -n | \
  awk '{a[NR]=$1} END {if(NR>0) print a[int(NR/2)]; else print "NA"}' \
    > "$OUT/median_depth.txt"

  ##########################################
  # 3. MISSINGNESS
  ##########################################
  echo "[STEP 3] Missingness"
  vcftools --gzvcf "$VCF" \
    --missing-site \
    --out "$OUT/missing"

  awk 'NR>1 {sum+=$6; n++} END {if(n>0) print sum/n; else print "NA"}' \
    "$OUT/missing.lmiss" > "$OUT/mean_missing_site.txt"

  ##########################################
  # 4. QUAL
  ##########################################
  echo "[STEP 4] QUAL"
  bcftools query -f '%QUAL\n' "$VCF" > "$OUT/qual_values.txt"

  awk '{sum+=$1; n++} END {if(n>0) print sum/n; else print "NA"}' \
    "$OUT/qual_values.txt" > "$OUT/mean_qual.txt"

  sort -n "$OUT/qual_values.txt" | \
  awk '{a[NR]=$1} END {if(NR>0) print a[int(NR/2)]; else print "NA"}' \
    > "$OUT/median_qual.txt"

  ##########################################
  # 5. MAF
  ##########################################
  echo "[STEP 5] MAF"
  vcftools --gzvcf "$VCF" \
    --freq \
    --out "$OUT/maf"

  awk 'NR>1 {
    split($5,a,":"); maf=a[2];
    if(maf<0.01) a1++;
    else if(maf<0.05) a2++;
    else if(maf<0.1) a3++;
    else if(maf<0.2) a4++;
    else a5++;
  }
  END{
    print "MAF<0.01:",a1;
    print "0.01-0.05:",a2;
    print "0.05-0.1:",a3;
    print "0.1-0.2:",a4;
    print ">0.2:",a5
  }' "$OUT/maf.frq" > "$OUT/maf_distribution.txt"

  ##########################################
  # 6. FILTER MISSING
  ##########################################
  echo "[STEP 6] Missing filter"
  awk 'NR>1 && $6<=0.1 {a++}
  NR>1 && $6<=0.2 {b++}
  END{
    print "max-missing 0.9:",a;
    print "max-missing 0.8:",b
  }' "$OUT/missing.lmiss" > "$OUT/missing_filter.txt"

  ##########################################
  # 7. FILTER MAF
  ##########################################
  echo "[STEP 7] MAF filter"
  awk 'NR>1 {
    split($5,a,":"); maf=a[2];
    if(maf>=0.05) b++;
    if(maf>=0.1) c++;
  }
  END{
    print "MAF >=0.05:",b;
    print "MAF >=0.1:",c
  }' "$OUT/maf.frq" > "$OUT/maf_filter.txt"

  echo "[INFO] QC selesai"
}

############################################
# OPTIONAL: RUN DIRECTLY
############################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  qc_vcf "${1:?VCF required}" "${2:?OUT required}"
fi