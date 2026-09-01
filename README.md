# Ultra-Lean Differential Expression Pipeline

[![Python Version](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Core Dependencies](https://img.shields.io/badge/Dependencies-pandas%20%7C%20numpy%20%7C%20matplotlib-brightgreen.svg)]()
[![Functional Enrichment](https://img.shields.io/badge/Enrichment-Enrichr%20REST%20API-purple.svg)](https://maayanlab.cloud/Enrichr/)

An ultra-lightweight, standalone differential expression (DE) and functional pathway analysis engine designed for bulk RNA-seq and single-cell pseudo-bulk transcriptomics.

---

## Abstract and Overview

Standard differential expression workflows frequently rely on complex runtime stacks with heavy third-party dependencies, compilation requirements, and platform-specific binaries. This pipeline provides a mathematically rigorous, self-contained implementation of differential expression analysis, quality control, dimensional reduction, and functional enrichment using only three core Python packages: `pandas`, `numpy`, and `matplotlib`.

All statistical distributions, matrix decompositions, multiple hypothesis adjustments, and web API queries are implemented natively or via Python standard library modules without requiring external specialized frameworks.

---

## Architectural Highlights

* **Minimal Dependency Stack**: Requires exclusively `pandas`, `numpy`, and `matplotlib`.
* **Zero Heavy Frameworks**:
  * **Principal Component Analysis (PCA)**: Computed natively via Singular Value Decomposition (`numpy.linalg.svd`), producing identical projections to `scikit-learn` without overhead.
  * **Exact $p$-Value Calculation**: Two-tailed $p$-values for Welch's $t$-distribution are computed directly from the regularized incomplete beta function using modified Lentz continued fraction expansions to machine precision ($< 10^{-15}$), removing `scipy` dependency.
  * **Multiple Testing Correction**: Vectorized step-up Benjamini–Hochberg False Discovery Rate (FDR) procedure in pure NumPy.
  * **Functional Pathway Queries**: Native REST client utilizing Python standard library `urllib.request` and `json` to query Enrichr databases (KEGG, GO Biological Process, Reactome, WikiPathways), eliminating `requests`.
  * **Statistical Engine**: Log2-Counts Per Million ($\log_2\text{CPM}$) variance-stabilized Welch's $t$-test providing robust sample-level and gene-level contrast testing without compilation bottlenecks.
* **Comprehensive Reporting**: Automatically renders publication-grade vector/raster figures, structured CSV tables, machine-readable JSON metadata, and an offline, self-contained HTML report with embedded Base64 graphics.
* **Audit and Provenance**: Exports SHA-256 checksums, conda environment definitions, and exact command invocations for complete analytical reproducibility.

---

## Mathematical Methodology

### 1. Library Size Scaling and Normalization
Raw sequencing counts $C_{g,s}$ for gene $g$ in sample $s$ are normalized for differences in sequencing depth using Counts Per Million (CPM) with a pseudo-count offset:

$$\text{CPM}_{g,s} = \frac{C_{g,s}}{\sum_{j} C_{j,s}} \times 10^6$$

$$Y_{g,s} = \log_2(\text{CPM}_{g,s} + 1)$$

### 2. Dimensional Reduction via SVD
Sample-level variance structure is evaluated on centered normalized expression matrix $\mathbf{Y}_{\text{centered}} = \mathbf{Y} - \boldsymbol{\mu}$:

$$\mathbf{Y}_{\text{centered}} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T$$

Principal component coordinates $\mathbf{Z} = \mathbf{U} \mathbf{\Sigma}$ and percentage of explained variance $\lambda_i = \sigma_i^2 / \sum \sigma_j^2$ are computed across orthogonal axes.

### 3. Welch's $t$-Statistic and Degrees of Freedom
For groups $A$ (numerator) and $B$ (denominator) with sample sizes $n_A, n_B$, sample means $\bar{Y}_A, \bar{Y}_B$, and sample variances $s_A^2, s_B^2$:

$$t = \frac{\bar{Y}_A - \bar{Y}_B}{\sqrt{\frac{s_A^2}{n_A} + \frac{s_B^2}{n_B}}}$$

Effective degrees of freedom $\nu$ are determined via the Welch–Satterthwaite approximation:

$$\nu = \frac{\left(\frac{s_A^2}{n_A} + \frac{s_B^2}{n_B}\right)^2}{\frac{(s_A^2/n_A)^2}{n_A - 1} + \frac{(s_B^2/n_B)^2}{n_B - 1}}$$

### 4. Continuous $p$-Value Evaluation via Incomplete Beta Function
The two-tailed survival probability for the calculated $t$-statistic is evaluated through the regularized incomplete beta function $I_x(a, b)$:

$$x = \frac{\nu}{\nu + t^2}, \quad p = I_x\left(\frac{\nu}{2}, \frac{1}{2}\right) = \frac{\text{B}\left(x; \frac{\nu}{2}, \frac{1}{2}\right)}{\text{B}\left(\frac{\nu}{2}, \frac{1}{2}\right)}$$

The incomplete beta integral is computed using the continued fraction expansion:

$$I_x(a, b) = \frac{x^a (1-x)^b}{a \, \text{B}(a, b)} \cdot \cfrac{1}{1 + \cfrac{d_1}{1 + \cfrac{d_2}{1 + \dots}}}$$

### 5. False Discovery Rate (FDR) Adjustment
Ranked raw $p$-values $p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$ for $m$ tested genes are adjusted using the Benjamini–Hochberg procedure:

$$\text{padj}_{(i)} = \min \left( 1, \, \min_{j \ge i} \left( \frac{m \cdot p_{(j)}}{j} \right) \right)$$

---

## Installation

### Prerequisites
* Python $\ge 3.10$

### Setup via pip
```bash
pip install pandas numpy matplotlib
```

### Setup via Conda
```bash
conda create -n depipeline -c conda-forge python=3.11 pandas numpy matplotlib -y
conda activate depipeline
```

---

## Quick Start

### 1. Verification with Built-in Demo Dataset
Run the integrated benchmark suite to verify installation integrity:

```bash
python3 de_pipeline.py --demo --output results/demo_run
```

### 2. Standard Analysis Invocation
Execute differential expression analysis on custom count matrices and sample sheets:

```bash
python3 de_pipeline.py \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --formula "~ batch + condition" \
  --contrast "condition,treated,control" \
  --output results/my_experiment
```

---

## Command-Line Interface (CLI) Specification

```
usage: de_pipeline.py [-h] [--counts COUNTS] [--metadata METADATA]
                      [--formula FORMULA] [--contrast CONTRAST]
                      [--min-count MIN_COUNT] [--min-samples MIN_SAMPLES]
                      --output OUTPUT [--demo] [--no-enrichment]
```

### Argument Reference

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--counts` | File path | `None` | Path to raw count matrix (CSV or TSV format; genes in rows, samples in columns). |
| `--metadata` | File path | `None` | Path to sample metadata table (CSV or TSV format; must include `sample_id` column). |
| `--formula` | String | `~ condition` | Design formula describing experimental factors (e.g., `~ condition` or `~ batch + condition`). |
| `--contrast` | String | `condition,treated,control` | Contrast definition formatted as `factor,numerator,denominator`. |
| `--min-count` | Integer | `10` | Minimum raw count threshold required for gene retention. |
| `--min-samples` | Integer | `2` | Minimum number of samples that must meet the `--min-count` threshold. |
| `--output` | Directory path | *(Required)* | Target directory for all figures, tables, reports, and checksums. |
| `--demo` | Flag | `False` | Executes analysis on the bundled 10-gene synthetic verification dataset. |
| `--no-enrichment` | Flag | `False` | Disables automated Enrichr REST API functional pathway queries. |

---

## Input File Format Requirements

### Count Matrix (`--counts`)
A comma- or tab-separated table where the first column contains unique gene symbols or Ensembl identifiers, and subsequent columns contain raw integer read counts.

```csv
gene,ctrl_1,ctrl_2,trt_1,trt_2
DUSP1,120,115,850,890
SPARCL1,45,50,720,690
GAPDH,8500,8620,8490,8530
```

### Metadata Table (`--metadata`)
A comma- or tab-separated table containing one row per sample. Must contain a `sample_id` column matching column names in the count matrix.

```csv
sample_id,condition,batch
ctrl_1,control,B1
ctrl_2,control,B2
trt_1,treated,B1
trt_2,treated,B2
```

---

## Output Structure and Artifacts

```
results/experiment_run/
├── report.html                  # Standalone interactive HTML report with embedded Base64 figures
├── report.md                    # Structured Markdown summary report
├── result.json                  # Machine-readable run execution metadata
├── figures/
│   ├── pca.png                  # Principal component analysis sample separation plot
│   ├── volcano.png              # Volcano plot (padj < 0.05, |log2FC| >= 1.0)
│   ├── ma_plot.png              # MA plot (log2FC vs log10 baseMean)
│   ├── heatmap.png              # Z-score normalized expression heatmap of top DE genes
│   └── enrichment_bubble.png    # Over-representation bubble plot for functional pathways
├── tables/
│   ├── qc_summary.csv           # Sample-level library sizes, detection rates, and medians
│   ├── normalized_counts.csv    # Variance-stabilized log2(CPM + 1) expression matrix
│   ├── de_results.csv           # Complete statistical results (baseMean, log2FC, stat, pvalue, padj)
│   └── enrichment_results.csv   # Over-representation statistics from Enrichr libraries
└── reproducibility/
    ├── commands.sh              # Exact shell invocation string
    ├── environment.yml          # Minimal Conda environment specification
    └── checksums.sha256         # SHA-256 cryptographic hashes of all input and output files
```

---

## Empirical Benchmark: NCBI GEO GSE52778

The pipeline was benchmarked against the NCBI Gene Expression Omnibus dataset **[GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)** (*Himes et al., 2014*), comprising human airway smooth muscle cell lines treated with dexamethasone versus untreated vehicle controls across four cell lines.

### Benchmark Results
* **Input Dimensions**: 8 samples $\times$ 35,268 genes (15,202 retained after filtering).
* **Execution Time**: Under 10 seconds for end-to-end QC, SVD PCA, differential expression, figure generation, and Enrichr query.
* **Statistical Output**: 131 significant differentially expressed genes ($\text{FDR} < 0.05, |\log_2\text{FC}| \ge 1.0$; 71 upregulated, 60 downregulated).
* **Validated Biological Targets**:
  * `SPARCL1`: $\log_2\text{FC} = +4.003, \, \text{padj} = 0.011$
  * `DUSP1`: $\log_2\text{FC} = +2.845, \, \text{padj} = 0.011$ (Dual specificity phosphatase 1; primary glucocorticoid anti-inflammatory mediator)
  * `PER1`: $\log_2\text{FC} = +2.643, \, \text{padj} = 0.011$ (Period circadian regulator 1)
  * `ARHGEF2`: $\log_2\text{FC} = -1.108, \, \text{padj} = 0.008$
* **Top Enriched Biological Pathways**:
  * Negative Regulation of p38MAPK Cascade (`DUSP1`, `DUSP10`)
  * Metallothionein Metal Binding and Zinc Ion Response (`MT2A`, `MT1X`, `MT1E`)
  * Negative Regulation of Leukocyte Chemotaxis (`DUSP1`, `CCN3`)

---

## License and Disclaimer

* **License**: This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
* **Disclaimer**: This software is intended exclusively for scientific research and educational purposes. It is not approved for clinical diagnostic procedures.
