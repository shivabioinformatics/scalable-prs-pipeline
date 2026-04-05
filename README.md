# Scalable Cloud Architecture for Polygenic Risk Scoring (PRS)

A complete, end-to-end Nextflow DSL2 pipeline that takes a patient VCF, performs quality control, normalizes to a common genome build, imputes missing genotypes, computes Polygenic Risk Scores using GWAS summary statistics, stores results in a database, and generates publication-quality visualizations.

---

## Problem Motivation

**Why Polygenic Risk Scores matter:**

Most common diseases — coronary artery disease, type 2 diabetes, breast cancer, Alzheimer's — aren't caused by a single gene. They're driven by the combined effects of thousands of common genetic variants, each contributing a tiny amount of risk. A Polygenic Risk Score (PRS) aggregates all of these small effects into a single number that estimates how much genetic risk a person carries.

**The clinical impact is real:**
- A person in the top 8% of PRS for coronary artery disease has the same risk as someone with a monogenic mutation in *LDLR* (familial hypercholesterolemia) — roughly 3x the population average
- Unlike family history, PRS is quantitative and actionable from birth
- The UK Biobank, Genomics England, and the NIH All of Us program are all integrating PRS into their research platforms

**Why a pipeline is needed:**
- Genotype data is massive: millions of variants × thousands of patients
- Different datasets use different genome builds (hg19 vs hg38) — coordinates must be harmonized
- Missing genotypes must be imputed using reference panels
- The scoring itself requires matching, aligning, and weighting thousands of variants
- Results need to be stored in a queryable database and visualized for clinical interpretation

Without a reproducible pipeline, every step is a manual, error-prone process. This project automates the entire workflow.

---

## Pipeline Architecture

```
                    ┌─────────────────────────────────────────────────────┐
                    │           PRS PIPELINE (Nextflow DSL2)             │
                    └─────────────────────────────────────────────────────┘
                                          │
          ┌───────────────────────────────┼───────────────────────────────┐
          ▼                               ▼                               ▼
    ┌───────────┐                  ┌────────────┐                  ┌───────────┐
    │  Input    │                  │   GWAS     │                  │ Reference │
    │  VCF      │                  │   Summary  │                  │   Panel   │
    │ (Patient) │                  │   Stats    │                  │  (1000G)  │
    └─────┬─────┘                  └─────┬──────┘                  └─────┬─────┘
          │                               │                               │
          ▼                               │                               │
    ┌───────────┐                         │                               │
    │  Step 1   │                         │                               │
    │    QC     │──► call rate, MAF,      │                               │
    │ BCFtools  │    depth filtering      │                               │
    └─────┬─────┘                         │                               │
          │                               │                               │
          ▼                               │                               │
    ┌───────────┐                         │                               │
    │  Step 2   │                         │                               │
    │ Normalize │──► liftover hg19→hg38   │                               │
    │ CrossMap  │    + allele alignment   │                               │
    └─────┬─────┘                         │                               │
          │                               │                               │
          ▼                               │                               │
    ┌───────────┐                         │                               │
    │  Step 3   │◄────────────────────────┼───────────────────────────────┘
    │  Impute   │──► fill missing calls
    │  Beagle   │    using LD patterns
    └─────┬─────┘
          │                               │
          ▼                               │
    ┌───────────┐                         │
    │  Step 4   │◄────────────────────────┘
    │   Score   │──► PRS = Σ(dosage × beta)
    │  Python   │    z-score normalization
    └─────┬─────┘
          │
          ├──────────────────────┐
          ▼                      ▼
    ┌───────────┐         ┌───────────┐
    │  Step 5   │         │  Step 6   │
    │ Database  │         │ Visualize │
    │  SQLite   │         │  Python   │
    └───────────┘         └───────────┘
          │                      │
          ▼                      ▼
    ┌───────────┐         ┌───────────┐
    │  3 Tables │         │  5 Plots  │
    │ • samples │         │ • QC      │
    │ • variants│         │ • Distrib │
    │ • qc_meta │         │ • Risk    │
    └───────────┘         │ • Manhtn  │
                          │ • Ranked  │
                          └───────────┘
```

---

## AWS Cloud Architecture

```
    ┌─────────────────────────────────────────────────────────────────────┐
    │                         AWS Cloud                                  │
    │                                                                     │
    │   ┌──────────┐     ┌──────────┐     ┌───────────────────────────┐  │
    │   │  S3      │────►│  Lambda  │────►│       AWS Batch           │  │
    │   │  Input   │     │ Trigger  │     │                           │  │
    │   │  Bucket  │     │          │     │  ┌──────────────────────┐ │  │
    │   │          │     └──────────┘     │  │   Nextflow Pipeline  │ │  │
    │   │ patient  │                      │  │                      │ │  │
    │   │  .vcf    │                      │  │  QC → Norm → Impute │ │  │
    │   └──────────┘                      │  │     → Score → DB    │ │  │
    │                                     │  └──────────────────────┘ │  │
    │                                     │                           │  │
    │                                     │  EC2 On-Demand  EC2 Spot  │  │
    │                                     │  (scoring)     (impute)   │  │
    │                                     └──────────┬────────────────┘  │
    │                                                │                   │
    │        ┌───────────────────┬────────────────────┤                   │
    │        ▼                   ▼                    ▼                   │
    │   ┌──────────┐     ┌────────────┐      ┌────────────┐             │
    │   │  S3      │     │    RDS     │      │ CloudWatch │             │
    │   │  Output  │     │ PostgreSQL │      │ Dashboard  │             │
    │   │  Bucket  │     │ (or SQLite │      │ + Alarms   │             │
    │   │          │     │  /DynamoDB)│      │            │             │
    │   │ results/ │     └────────────┘      └────────────┘             │
    │   │ figures/ │                                                     │
    │   └──────────┘          ┌────────────┐                            │
    │        │                │ AWS Budget │                            │
    │        │ Lifecycle      │ $100/month │                            │
    │        │ Rules          │ + alerts   │                            │
    │        ▼                └────────────┘                            │
    │   ┌──────────┐                                                    │
    │   │ S3       │                                                    │
    │   │ Glacier  │  ◄── auto-archive after 90 days                   │
    │   │ Archive  │      HIPAA: retain ≥ 6 years                      │
    │   └──────────┘                                                    │
    │                                                                     │
    └─────────────────────────────────────────────────────────────────────┘
```

**Key Design Decisions:**

| Component | Choice | Why |
|---|---|---|
| Compute | AWS Batch | Auto-scales EC2 instances, shuts down when done — no idle costs |
| Spot Instances | Imputation step only | 60-90% cheaper; Nextflow retries on interruption |
| Storage | S3 + Lifecycle Rules | Standard → IA (30d) → Glacier (90d) → Deep Archive (365d) |
| Database | RDS PostgreSQL / DynamoDB | Variant weights in DynamoDB for O(1) lookups; scores in RDS for SQL queries |
| Trigger | S3 → Lambda | Event-driven: drop a VCF, pipeline runs automatically |
| Monitoring | CloudWatch + Budget | Dashboards for memory/CPU, $100/month budget with email alerts |

---

## Variant Selection for the PRS Model

Not all GWAS variants should go into a PRS. The selection strategy matters:

**P-value Thresholding (C+T Method):**
- Include only variants below a significance threshold (e.g., p < 5×10⁻⁸ for genome-wide significant, or p < 0.05 for more liberal)
- Simple but effective. This pipeline uses p-value filtering by default.

**LD Clumping:**
- Nearby variants are often correlated (linkage disequilibrium)
- Including both would double-count the same genetic signal
- PLINK `--clump` selects the most significant variant in each LD block
- In production, I'd add a PLINK clumping step before scoring

**Advanced Methods (not implemented, but worth knowing):**
- **LDpred2**: Bayesian method that adjusts betas for LD structure
- **PRS-CS**: Uses continuous shrinkage priors on effect sizes
- **These give better prediction than simple C+T but are computationally heavier**

---

## Normalization — What It Means in PRS Context

Normalization happens at two levels:

### 1. Coordinate Normalization (Liftover)
Different datasets are built on different reference genomes (hg19 vs hg38). If my patient data is on hg19 and the GWAS summary stats are on hg38, the positions won't match. CrossMap converts coordinates using a chain file that maps old positions to new positions.

### 2. Batch / Population Normalization
A raw PRS of 0.47 means nothing by itself. It only has meaning relative to a population:
- A European with PRS = 0.47 might be at the 60th percentile in a European reference
- An African American with PRS = 0.47 might be at the 90th percentile in an African reference

**Z-score normalization** within ancestry groups fixes this:
```
z_score = (raw_prs - ancestry_mean) / ancestry_std
```

This pipeline computes z-scores across all input samples. In production with known ancestry labels, you'd normalize within each group separately.

---

## S3 Lifecycle Rules & Data Archiving

This addresses the feedback about discussing archiving:

| Time After Pipeline Run | S3 Storage Class | Cost (per TB/month) | Access Time |
|---|---|---|---|
| 0–30 days | STANDARD | $23 | Instant |
| 30–90 days | STANDARD_IA | $12.50 | Instant |
| 90–365 days | GLACIER | $4 | 3-5 hours |
| 365+ days | DEEP_ARCHIVE | $1 | 12 hours |

**Why this matters:**
- Genomic data is large — a single whole-genome VCF is ~100GB
- HIPAA requires retaining patient data for ≥6 years
- Active results need fast access, but year-old runs can be archived
- Lifecycle rules automate this — no manual intervention needed

The lifecycle policy is defined in `nextflow.config` and `conf/aws.config`.

---

## Project Structure

```
prs_pipeline/
├── main.nf                      # Pipeline entry point (Nextflow DSL2)
├── nextflow.config              # Parameters, profiles, S3 lifecycle config
├── README.md                    # This file
│
├── modules/                     # One Nextflow process per file
│   ├── qc.nf                   # BCFtools-based QC filtering
│   ├── normalize.nf            # CrossMap liftover + batch normalization
│   ├── impute.nf               # Beagle genotype imputation
│   ├── score.nf                # PRS weighted sum calculation
│   ├── database.nf             # SQLite database storage
│   └── visualize.nf            # Matplotlib visualization
│
├── subworkflows/
│   └── prs_workflow.nf         # Chains all modules: QC → Norm → Impute → Score → DB → Viz
│
├── bin/                         # Python scripts (auto-added to PATH in processes)
│   ├── generate_test_data.py   # Creates simulated VCF + GWAS stats for demo
│   ├── qc_filter.py            # QC: call rate, MAF, depth filtering
│   ├── calculate_prs.py        # Core PRS scoring engine
│   ├── store_results.py        # SQLite database operations
│   └── visualize_prs.py        # 5 publication-quality plots
│
├── conf/
│   ├── resources.config        # Per-process CPU/memory allocation
│   └── aws.config              # AWS Batch, S3, CloudWatch configuration
│
├── data/                        # Test data (generated by generate_test_data.py)
│   ├── sample_input.vcf        # 500 variants × 20 patients (~5% missing)
│   ├── gwas_summary_stats.tsv  # Effect sizes + p-values for 500 variants
│   ├── reference_panel.vcf     # 50-sample reference for imputation
│   └── chain_file.chain        # Placeholder liftover chain
│
└── results/                     # Pipeline outputs (generated by running the pipeline)
    ├── qc/                     # Filtered VCF + QC stats
    ├── normalized/             # Build-normalized VCF
    ├── imputed/                # Imputed VCF
    ├── scores/                 # PRS scores per patient
    ├── database/               # SQLite database
    └── figures/                # 5 visualization PNGs
```

---

## How to Run

### Prerequisites
- Python 3.8+ with `matplotlib` and `numpy`
- Nextflow 22.10+ (for running the full pipeline)

### Quick Start (Python scripts only — no Nextflow needed)

```bash
# 1. generate test data
cd prs_pipeline
python3 bin/generate_test_data.py

# 2. run QC
python3 bin/qc_filter.py \
  --vcf data/sample_input.vcf \
  --out-vcf results/qc/filtered.vcf \
  --out-stats results/qc/qc_stats.tsv

# 3. calculate PRS (skipping normalize/impute for simplicity)
python3 bin/calculate_prs.py \
  --vcf results/qc/filtered.vcf \
  --gwas data/gwas_summary_stats.tsv \
  --output results/scores/prs_scores.tsv

# 4. store in database
python3 bin/store_results.py \
  --scores results/scores/prs_scores.tsv \
  --gwas data/gwas_summary_stats.tsv \
  --qc-stats results/qc/qc_stats.tsv \
  --db results/database/prs_results.db

# 5. generate visualizations
python3 bin/visualize_prs.py \
  --scores results/scores/prs_scores.tsv \
  --qc-stats results/qc/qc_stats.tsv \
  --gwas data/gwas_summary_stats.tsv \
  --outdir results/figures
```

### Full Pipeline (with Nextflow)

```bash
# local execution
nextflow run main.nf -profile local

# docker execution (more reproducible)
nextflow run main.nf -profile docker

# on AWS
nextflow run main.nf -profile aws_batch \
  --input_vcf s3://my-bucket/patient.vcf \
  --outdir s3://my-bucket/results
```

---

## Outputs

| Output | Location | Description |
|---|---|---|
| Filtered VCF | `results/qc/filtered.vcf` | Variants passing QC thresholds |
| QC Report | `results/qc/qc_stats.tsv` | Per-variant and per-sample QC metrics |
| Normalized VCF | `results/normalized/normalized.vcf` | Genome build–harmonized variants |
| Imputed VCF | `results/imputed/imputed.vcf` | Missing genotypes filled in |
| PRS Scores | `results/scores/prs_scores.tsv` | Per-patient raw score, z-score, percentile, risk category |
| Database | `results/database/prs_results.db` | SQLite with samples, variants, qc_metrics tables |
| QC Summary Plot | `results/figures/01_qc_summary.png` | Call rates and depth per sample |
| PRS Distribution | `results/figures/02_prs_distribution.png` | Histogram with risk thresholds |
| Risk Stratification | `results/figures/03_risk_stratification.png` | Patients per risk category |
| Variant Effects | `results/figures/04_variant_effects.png` | Manhattan-style effect size plot |
| Score Comparison | `results/figures/05_score_comparison.png` | Ranked lollipop chart |

---

## Database Schema

```sql
-- patient PRS scores (what clinicians query)
SELECT sample_id, prs_zscore, percentile, risk_category
FROM samples
WHERE risk_category = 'HIGH_RISK';

-- most significant variants in the model
SELECT rsid, chromosome, beta, p_value
FROM variants
ORDER BY p_value ASC
LIMIT 10;

-- samples with questionable QC
SELECT s.sample_id, s.prs_zscore, q.call_rate
FROM samples s
JOIN qc_metrics q ON s.sample_id = q.sample_id
WHERE q.call_rate < 0.95;
```

---

## Cost Management

| Strategy | Tool | Details |
|---|---|---|
| Budget alerts | AWS Budgets | $100/month cap with email notifications |
| Spot instances | AWS Batch | 60-90% savings on imputation compute |
| Auto-shutdown | AWS Batch | EC2 instances terminate when pipeline finishes |
| Data archiving | S3 Lifecycle | Auto-transition to Glacier after 90 days |
| Monitoring | CloudWatch | CPU/memory dashboards + duration alarms |

**Estimated cost per run** (1M variants × 100 samples): **~$0.36** with spot instances.

---

## Tools & References

| Tool | Purpose | Reference |
|---|---|---|
| Nextflow | Workflow orchestration | [nextflow.io](https://nextflow.io) |
| BCFtools | VCF quality control | [samtools.github.io/bcftools](https://samtools.github.io/bcftools) |
| CrossMap | Genome build liftover | [crossmap.sourceforge.net](http://crossmap.sourceforge.net) |
| Beagle | Genotype imputation | [faculty.washington.edu/browning/beagle](https://faculty.washington.edu/browning/beagle) |
| PGS Catalog | Published PRS models | [pgscatalog.org](https://www.pgscatalog.org) |
| pgsc_calc | Reference PRS pipeline | [github.com/PGScatalog/pgsc_calc](https://github.com/PGScatalog/pgsc_calc) |
| matplotlib | Visualization | [matplotlib.org](https://matplotlib.org) |
| SQLite | Local database | [sqlite.org](https://sqlite.org) |

---

## License

This project was built for BINFX 410 — Final Project.
