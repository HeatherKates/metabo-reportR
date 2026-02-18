# SECIM Metabolomics Analysis Pipeline

**Southeast Center for Integrated Metabolomics (SECIM)**  
University of Florida

A modular, reproducible pipeline for untargeted metabolomics data analysis with interactive reporting.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Directory Structure](#directory-structure)
- [Installation](#installation)
- [Input File Requirements](#input-file-requirements)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Overview

This pipeline provides end-to-end analysis of untargeted LC-MS metabolomics data, from raw peak tables (MZmine output) to publication-ready interactive reports. The workflow is designed for reproducibility and portability across computing environments.

**Two-Stage Architecture:**

1. **Stage 1 (Computation):** Data processing, normalization, statistical analysis, and metabolite annotation
2. **Stage 2 (Reporting):** Interactive HTML report generation with Plotly visualizations and searchable tables

---

## Features

- ✅ **Supports multiple statistical tests:** t-test, ANOVA, linear models, linear mixed models
- ✅ **Automated metabolite annotation:** SECIM internal library + MS1 database matching (metid)
- ✅ **Interactive visualizations:** PCA, volcano plots, heatmaps (Plotly)
- ✅ **Batch correction:** ComBat or limma methods
- ✅ **Portable configuration:** YAML-based study parameters
- ✅ **Fast serialization:** Uses `qs` package for efficient data storage
- ✅ **Self-contained reports:** Single HTML file with embedded resources

---

## Directory Structure

```
SECIM_Metabolomics_Pipeline/
├── Run_Pipeline.R              # Master execution script
├── config/
│   ├── study_params.yaml       # Template configuration file
│   └── [your_study].yaml       # Study-specific configs
├── R/
│   ├── Generate_Reporting_Inputs.R   # Main processing function
│   ├── SECIM_Metabolomics.R          # Statistical analysis engine
│   ├── Norm_Plots.R                  # Normalization QC plots
│   ├── Generate_network.R            # Pathway network visualization
│   ├── general_data_utils.R          # MetaboAnalystR utilities
│   ├── metid_annotation.R            # Custom metid functions
│   └── SanityCheck.HRK.R             # Data validation
├── Templates/
│   └── Report_Template.qmd           # Quarto report template
├── InputFiles/
│   └── [your_data].xlsx              # MZmine output (user-provided)
├── reference/
│   ├── library_positive1_MK_QE2.csv  # SECIM library (Pos, Evospere-PFP)
│   ├── Positive_ACE-PFPs.csv         # SECIM library (Pos, ACE-PFP)
│   ├── library_negative1_MK02.csv    # SECIM library (Neg, Evospere-PFP)
│   ├── Negative_ACE-PFP.csv          # SECIM library (Neg, ACE-PFP)
│   ├── kegg_ms1_database0.0.3.rda    # Metid database (KEGG)
│   ├── hmdb_database0.0.3.rda        # Metid database (HMDB)
│   └── bloodexposome_database1.0.rda # Metid database (Blood Exposome)
└── Results/
    └── [client_name]/
        ├── [client]_ReportData.qs    # Serialized analysis results
        └── [client]_Report.html      # Interactive HTML report
```

---

## Installation

### System Requirements

- **R version:** ≥ 4.2.0
- **Operating System:** Linux, macOS, or Windows
- **RAM:** ≥ 16 GB recommended
- **Disk Space:** ≥ 10 GB for reference databases

### R Package Installation

```r
# Install required CRAN packages
install.packages(c(
  "yaml", "qs", "dplyr", "tidyr", "stringr", "readxl",
  "ggplot2", "plotly", "DT", "heatmaply", "pheatmap",
  "RColorBrewer", "gridExtra", "car", "foreach", "parallel",
  "broom", "emmeans", "lme4", "lmerTest", "nlme",
  "impute", "htmltools", "downloadthis", "quarto"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "MetaboAnalystR", "KEGGREST", "pathview", "org.Hs.eg.db",
  "SummarizedExperiment", "S4Vectors"
))

# Install GitHub packages
remotes::install_github("tidymass/metid")
remotes::install_github("andreasmock/MetaboDiff")
remotes::install_github("cran/omu")
```

### Quarto Installation

Download and install Quarto from [https://quarto.org/docs/get-started/](https://quarto.org/docs/get-started/)

---

## Input File Requirements

### Excel File Structure

Your input Excel file **must** contain exactly **three sheets**:

#### 1. `Sample.data` (Metadata)

| Sample.Name | Class | Subject | Batch |
|-------------|-------|---------|-------|
| Sample_001  | Control | P001 | 1 |
| Sample_002  | Treatment | P002 | 1 |
| Sample_003  | Control | P003 | 2 |

- **Required columns:**
  - `Sample.Name`: Unique sample identifier (must match peak table column names)
  - `Class`: Grouping variable for statistical comparisons
  
- **Optional columns:**
  - `Subject`: For paired analyses or random effects models
  - `Batch`: For batch correction

#### 2. `Peaktable.neg` (Negative Ion Mode)

MZmine output with the following structure:

| row ID | row m/z | row retention time | compound | Sample_001 | Sample_002 | ... |
|--------|---------|--------------------|-----------:|------------|------------|-----|
| 1 | 180.0634 | 1.23 | Glucose | 123456 | 234567 | ... |
| 2 | 146.0459 | 2.45 | | 45678 | 56789 | ... |

- **MZmine 2 format:** `row ID`, `row m/z`, `row retention time`, `compound`, [samples]
- **MZmine 3 format:** `id`, `rt`, `mz`, `internal_ID`, `spectral_ID`, [samples]

#### 3. `Peaktable.pos` (Positive Ion Mode)

Same structure as `Peaktable.neg`

### Important Notes

- **No special characters** in sample names except `_` and `.`
- **Sample names** must match between metadata and peak tables
- **Column names** must not be duplicated
- **Numeric data only** in peak intensity columns

---

## Configuration

### Creating a Configuration File

Copy the template and modify for your study:

```bash
cp config/study_params.yaml config/MyStudy_2024.yaml
```

### Example Configuration

```yaml
# Project Information
project:
  client: "Smith_TBI_2024"
  PI_name: "Dr. Jane Smith"

# Input Files
input:
  file: "InputFiles/Smith_TBI_MZmine3.xlsx"
  mzmine_version: 3

# Sample Filtering
filtering:
  samples_to_drop_pre_norm:
    - "Sample_045"  # Failed QC
    - "Sample_023"  # Outlier
  samples_to_drop_post_norm: null

# Statistical Configuration
statistics:
  test_type: "t.test"  # Options: t.test, anova, lm, lme, nostats
  contrast_var: "Condition"
  reference_level: "Control"
  paired: true
  num_meta: 1
  subset: null
  anova_formula: null
  lm_model: null

# Metabolite Annotation
annotation:
  SECIM_column: "ACE-PFP"  # Options: ACE-PFP, Evospere-PFP
  metid_DB_file: "kegg_ms1_database0.0.3.rda"

# Data Processing
processing:
  batch_correct: "ComBat"  # Options: null, ComBat, limma
  rowNorm: "SumNorm"
  transNorm: "LogNorm"
  scaleNorm: "ParetoNorm"
```

### Configuration Parameters

#### Statistical Test Types

- **`t.test`**: Pairwise comparison (paired or unpaired)
- **`anova`**: Analysis of variance with post-hoc pairwise comparisons
- **`lm`**: Linear model
- **`lme`**: Linear mixed effects model
- **`nostats`**: Fold-change analysis only

#### Advanced Model Formulas

For `anova`:
```yaml
anova_formula: "id ~ Class + Error(Subject)"
```

For `lm` or `lme`:
```yaml
lm_model: "Metabolite ~ Class + (1|Subject)"
```

#### Normalization Methods

- **rowNorm:** `SumNorm`, `MedianNorm`, `QuantileNorm`, `CompNorm`, `SpecNorm`
- **transNorm:** `LogNorm`, `CrNorm` (Cubic Root)
- **scaleNorm:** `ParetoNorm`, `AutoNorm`, `MeanCenter`, `RangeNorm`

---

## Usage

### Basic Usage

```bash
# Using default config location (config/study_params.yaml)
Rscript Run_Pipeline.R

# Using a specific config file
Rscript Run_Pipeline.R config/MyStudy_2024.yaml
```

### From RStudio

```r
# Set working directory to project root
setwd("/path/to/SECIM_Metabolomics_Pipeline")

# Run with default config
source("Run_Pipeline.R")

# Or specify config
args <- "config/MyStudy_2024.yaml"
source("Run_Pipeline.R")
```

### Re-rendering Report Only

If you've already run Stage 1 and just want to update the report:

```r
library(quarto)

quarto_render(
  input = "Templates/Report_Template.qmd",
  output_file = "Smith_TBI_2024_Report.html",
  execute_params = list(
    client = "Smith_TBI_2024",
    PI_name = "Dr. Jane Smith",
    test_type = "t.test",
    contrast_var = "Condition",
    data_file = "Results/Smith_TBI_2024/Smith_TBI_2024_ReportData.qs"
  )
)
```

---

## Output Files

### Results Directory Structure

```
Results/
└── [client_name]/
    ├── [client]_ReportData.qs        # Serialized data (Stage 1 output)
    ├── [client]_Report.html          # Interactive HTML report
    └── temp/                         # Temporary files (auto-cleaned)
```

### Report Contents

The HTML report includes:

1. **Project Summary:** Sample counts, PI information, processing time
2. **Methods:** Instrumental details, normalization, statistical approach
3. **PCA Plot:** Interactive principal component analysis
4. **Results Table:** Searchable/sortable table with Excel export
5. **Volcano Plots:** Interactive volcano plots per contrast
6. **Heatmaps:** Top changed metabolites
7. **Normalization QC:** Before/after normalization plots
8. **Session Info:** R package versions for reproducibility

### Excel Download

The "Download Results" button in the report generates an Excel file with:

- **README:** Column definitions
- **report_results:** Combined statistical results
- **Pos.emmeans.results.metab / Neg.emmeans.results.metab:** Mode-specific results
- **Pos.processed.data / Neg.processed.data:** QC-filtered peak intensities
- **Pos.normalized.data / Neg.normalized.data:** Normalized peak intensities
- **metadata:** Sample information

---

## Troubleshooting

### Common Issues

#### 1. "Data file not found" error

**Problem:** The `.qs` file path is incorrect

**Solution:**
```r
# Check if file exists
file.exists("Results/MyStudy/MyStudy_ReportData.qs")

# Verify path in config matches actual location
```

#### 2. "Multiple or no matches found for sample name"

**Problem:** Sample names in metadata don't match peak table columns

**Solution:**

- Ensure exact match (case-sensitive)
- Check for extra spaces or special characters
- Verify MZmine prefix format matches your data

#### 3. Memory errors during processing

**Problem:** Large dataset exceeds available RAM

**Solution:**

```r
# Increase memory limit (Windows)
memory.limit(size = 32000)

# Use fewer samples or filter features more aggressively
```

#### 4. Quarto rendering fails

**Problem:** Quarto not installed or not in PATH

**Solution:**

```bash
# Install Quarto
# Visit: https://quarto.org/docs/get-started/

# Verify installation
quarto --version

# In R, check path
Sys.which("quarto")
```

#### 5. Reference database not found

**Problem:** Metabolite database file missing

**Solution:**

```bash
# Check reference directory
ls reference/

# Download missing databases
# Contact SECIM or visit https://github.com/tidymass/metid
```

---

## Best Practices

### 1. Version Control Your Config Files

```bash
git add config/MyStudy_2024.yaml
git commit -m "Add config for Smith TBI study"
```

### 2. Document Sample Exclusions

Always document why samples were excluded:

```yaml
filtering:
  samples_to_drop_pre_norm:
    - "Sample_045"  # Failed QC: RSD > 30%
    - "Sample_023"  # Technical replicate outlier
```

### 3. Archive Raw Data

Keep a backup of:
- Original MZmine output
- Raw `.mzML` files
- Configuration file used

### 4. Test with Subset First

For large studies, test with a subset:

```yaml
filtering:
  samples_to_drop_pre_norm: [all but 10 samples]
```

---

## Pipeline Validation

### Minimal Working Example

Test the pipeline with example data:

```bash
# 1. Create test config
cat > config/test.yaml << EOF
project:
  client: "Test_2024"
  PI_name: "Test User"
input:
  file: "InputFiles/test_data.xlsx"
  mzmine_version: 3
statistics:
  test_type: "t.test"
  contrast_var: "Class"
annotation:
  SECIM_column: "ACE-PFP"
  metid_DB_file: "kegg_ms1_database0.0.3.rda"
processing:
  batch_correct: null
EOF

# 2. Run pipeline
Rscript Run_Pipeline.R config/test.yaml

# 3. Check output
ls Results/Test_2024/
```

Expected output:
```
Test_2024_ReportData.qs
Test_2024_Report.html
```

---

## Performance Benchmarks

| Dataset Size | Samples | Features | Processing Time | Peak RAM |
|--------------|---------|----------|-----------------|----------|
| Small        | 20      | 500      | ~2 min          | 2 GB     |
| Medium       | 50      | 1500     | ~8 min          | 6 GB     |
| Large        | 100     | 3000     | ~20 min         | 12 GB    |
| Very Large   | 200     | 5000     | ~45 min         | 24 GB    |

*Benchmarks on Intel Xeon 2.4 GHz, 32 GB RAM*

---

## Contributing

### Reporting Issues

Please report bugs or request features via:

- Email: secim@metabolomics.ufl.edu
- GitHub Issues: [repository URL]

### Development Workflow

```bash
# 1. Create feature branch
git checkout -b feature/new-stat-method

# 2. Make changes
# Edit R/SECIM_Metabolomics.R

# 3. Test with example data
Rscript Run_Pipeline.R config/test.yaml

# 4. Commit and push
git add R/SECIM_Metabolomics.R
git commit -m "Add support for Kruskal-Wallis test"
git push origin feature/new-stat-method
```

---

## Citation

If you use this pipeline in your research, please cite:

```
Southeast Center for Integrated Metabolomics (SECIM)
University of Florida
Metabolomics Analysis Pipeline v2.0
https://secim.ufl.edu
```

### Key Dependencies to Cite

- **MetaboAnalystR:** Pang et al. (2021) *Nucleic Acids Research*
- **metid:** Shen et al. (2022) *Nature Communications*
- **emmeans:** Lenth (2023) R package version 1.8.4
- **Quarto:** Allaire et al. (2023) https://quarto.org

---

## License

This pipeline is provided for academic and research use.

For commercial licensing inquiries, contact: info@SECIM.ufl.edu

---

## Contact

**Southeast Center for Integrated Metabolomics (SECIM)**  
University of Florida  
Email:  info@SECIM.ufl.edu
Web: https://secim.ufl.edu

**Technical Support:**  
Heather Kates, Bioinformatics Analyst
Email: hkates@ufl.edu

---

## Changelog

### Version 2.0 (2026)

- ✅ Migrated to YAML configuration files
- ✅ Implemented two-stage architecture (compute + report)
- ✅ Added Quarto-based interactive reporting
- ✅ Removed all hardcoded paths for portability
- ✅ Improved metabolite annotation pipeline
- ✅ Added support for linear mixed models

### Version 1.0 (2024)

- Initial RMarkdown-based pipeline
