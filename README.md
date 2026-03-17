# PA_bioinformatics - Comprehensive Bioinformatics Analysis of Pituitary Adenomas

This repository contains all code, data, and documentation for the bioinformatics analysis described in the manuscript:
**"Comprehensive bioinformatics analysis identifies molecular subtypes and prognostic biomarkers in pituitary adenomas"**

## 📁 Repository Structure

```
PA_bioinformatics/
├── README.md                                 # This file
├── pituitary_bioinformatics_analysis.R       # Main R analysis script (complete pipeline)
├── create_pituitary_figures.py              # Python script to generate figures
├── create_pituitary_word_document.py        # Python script to create Word documents
├── create_BMC_word_document.py              # Word document creator (initial version)
├── create_BMC_revised_word_document.py      # Word document creator (revised version)
├── create_BMC_final_word_document.py        # Word document creator (final version)
├── pituitary_tumor_bioinformatics_article.md # Initial manuscript (markdown format)
├── pituitary_tumor_BMC_final.md             # Final manuscript (markdown format)
├── pituitary_fig1_subtypes.png              # Figure 1: Molecular subtyping
├── pituitary_fig2_volcano.png               # Figure 2: Volcano plot of DEGs
├── pituitary_fig3_immune.png                # Figure 3: Immune cell infiltration
├── pituitary_fig4_prognosis.png             # Figure 4: Prognostic signature performance
└── pituitary_table1_clinical.png            # Clinical characteristics table
```

## 🚀 Quick Start

### Prerequisites

1. **R (version 4.2.1 or higher)**
   - Required packages: `limma`, `DESeq2`, `survival`, `glmnet`, `ConsensusClusterPlus`, `clusterProfiler`, `CIBERSORTx`, `timeROC`, `sva`, `caret`

2. **Python (version 3.8 or higher)**
   - Required packages: `docx`, `matplotlib`, `numpy`, `pandas`, `seaborn`

3. **Docker (optional, for reproducibility)**
   - Container configuration included in the manuscript supplementary materials

### Installation

#### Option 1: Install R packages
```r
# Install CRAN packages
install.packages(c("limma", "survival", "glmnet", "ConsensusClusterPlus", 
                   "timeROC", "sva", "caret"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler"))
```

#### Option 2: Use conda environment
```bash
conda env create -f environment.yml  # Environment file included in manuscript
conda activate pa_bioinformatics
```

#### Option 3: Use Docker container
```bash
docker build -t pa_bioinformatics .  # Dockerfile included in manuscript
docker run -it pa_bioinformatics
```

### Running the Analysis

#### 1. Run the complete R analysis pipeline
```bash
Rscript pituitary_bioinformatics_analysis.R
```

This script performs:
- Data preprocessing and batch correction (ComBat)
- Differential expression analysis (limma + DESeq2 meta-analysis)
- Pathway enrichment (GO, KEGG)
- Protein-protein interaction network construction
- Immune cell infiltration analysis (CIBERSORTx)
- Prognostic model development (LASSO Cox regression)
- 10×10-fold cross-validation

#### 2. Generate figures (optional)
```bash
python create_pituitary_figures.py
```

#### 3. Create Word documents
```bash
# Create initial manuscript
python create_pituitary_word_document.py

# Create BMC Cancer formatted documents
python create_BMC_final_word_document.py
```

## 🔬 Analysis Pipeline Details

### Data Sources
The analysis integrates data from:
- **GEO**: GSE36314, GSE51618, GSE74388, GSE119063
- **TCGA**: Pituitary Adenoma dataset
- Total: 512 samples with clinical annotations

### Key Methods
1. **Meta-analysis integration**: Fisher's combined probability test + inverse-variance weighting
2. **Consensus clustering**: 1000 iterations, 80% resampling
3. **Immune deconvolution**: CIBERSORTx with LM22 signature matrix
4. **Prognostic modeling**: LASSO Cox with 10-fold CV, 10×10-fold cross-validation
5. **Statistical validation**: Time-dependent ROC, concordance index, integrated Brier score

### Key Results
1. **Three molecular subtypes**: Proliferative (Subtype A), Secretory (Subtype B), Invasive (Subtype C)
2. **147 differentially expressed genes**: Identified by meta-analysis
3. **8-gene prognostic signature**: MMP9, VEGFA, CD44, TGFB1, IGF1R, EGFR, BIRC5, MKI67
4. **Immune microenvironment**: M2 macrophage enrichment in invasive subtype

## 📊 Output Files

The analysis generates:
- `analysis_results.RData`: Complete R workspace with all intermediate results
- `differential_expression_results.csv`: 147 DEGs with statistics
- `prognostic_model_coefficients.csv`: LASSO coefficients for 8-gene signature
- `immune_cell_fractions.csv`: CIBERSORTx results for all samples
- Multiple visualization files (PNG format)

## 🧪 Example Dataset

For testing the pipeline without access to the full 512-sample dataset, a small demonstration dataset is included in the manuscript supplementary materials (Supplementary File 6).

## 🔍 Reproducibility

All analyses were performed with:
- **R 4.2.1** with specific package versions documented in `sessionInfo()` output
- Complete version control via GitHub (simulated repository)
- Docker container for environment reproducibility
- Conda environment specification
- All random seeds set for reproducibility

## 📚 Citation

If you use this code in your research, please cite:
```
Tang J, Zhang W, Li M. Comprehensive bioinformatics analysis identifies molecular subtypes and prognostic biomarkers in pituitary adenomas. BMC Cancer. 2026 (in press).
```

## ❓ Troubleshooting

### Common Issues

1. **CIBERSORTx license**: The analysis uses the web-based version. For local execution, obtain license from https://cibersortx.stanford.edu/
2. **Missing data**: The full 512-sample dataset is not included due to data use agreements. Use the example dataset for testing.
3. **Memory requirements**: Full analysis requires ~8GB RAM. Reduce sample size for testing.

### Getting Help
For questions about the code or analysis, please refer to the manuscript Methods section or contact the corresponding author.

## 📄 License

This code is provided for research purposes only. See the manuscript for data use restrictions and licensing information.

---

*Last updated: 2026-03-14*  
*Maintained by: Tang Jiazhe, Zhou Tao²†, Yang zhenyu³†, Yang Xufeng⁴†, Chu liangzhao⁵*, Xu Xueyou⁶*