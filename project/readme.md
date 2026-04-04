## Usage (Downstream Analysis)

Run QC first, then run the main pipeline.

---

### 1. QC Simulation (to determine thresholds)

Run QC on each VCF:

```bash
bash qc_vcf.sh data/madura.vcf.gz qc/madura/
bash qc_vcf.sh data/bosjava.vcf.gz qc/bosjava/
bash qc_vcf.sh data/india.vcf.gz qc/india/
```

Use the outputs (e.g. `maf_distribution.txt`, `missing_filter.txt`) to decide:

* MAF threshold (e.g. 0.01 or 0.05)
* Missing threshold (e.g. 0.9 or 0.95)

---

### 2. Run Main Pipeline

```bash
bash analysis.sh \
  data/madura.vcf.gz \
  data/bosjava.vcf.gz \
  data/india.vcf.gz \
  input/pelaihari.txt \
  input/madura.txt \
  input/meta_data.csv \
  input/Bos_taurus.gtf \
  0.05 0.9
```

Where:

* `0.05` = MAF threshold
* `0.9` = max-missing threshold

---

QC is used only to guide threshold selection; the pipeline applies those thresholds for analysis.
