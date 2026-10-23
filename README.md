# Code Repository

## Manuscript

Computational comparative genomic analysis and module reconstruction of Lactiplantibacillus plantarum reveals its neurometabolic potential at strain level

---

## Description

This repository contains the custom Python scripts used for data processing,
statistical analysis, and figure generation for the above manuscript.

---

## Software

Python 3.12.13

---

## Required Python Packages

- pandas
- numpy
- matplotlib
- scipy
- openpyxl

---

## Folder Structure

```
Data/
│
├── README.md
├── 01_pangenome_curve.py
├── 02_heatmap1.py
├── 03_heatmap2.py
├── 04_bacteriocin_curve.py
├── 05_bacteriocin_graph.py
├── 06_Module_completeness.py
├── 07_chi_square.py
├── 08_COG_category.py
├── Supplementary data/
└──	Table1.xlsx
        Table2.xlsx
	Table3.xlsx
	Table4.xlsx
	Table5.xlsx
        Table6.xlsx
	Table7.xlsx
	Table8.xlsx
	Table legends.docx
```

---

## Script Description

### 01_pangenome_curve.py

Generates the pangenome accumulation curve, pie chart of gene classification and gene frequency histogram.

### 02_heatmap1.py

Generates neuroactive gene presence/absence heatmaps.

### 03_heatmap2.py

Generate neuroactive gene presence/absence heatmap on isolation sources.

### 04_Bacteriocin_curve.py

Generate bacteriocin graph for bacteriocin burden per genome

### 05_Bacteriocin_graph.py

Generate bacteriocin heatmap on bacteriocin classes, bacteriocin richness, class 11 bacteriocin richness.

### 06_GBM_module_completeness.py

Calculates GBM module completeness analysis

### 07_Statistical_Analysis.py

Provide module completeness of all 318 genomes and statistical analysis to analyze prevalence of modules.

### 08_COG_category.py

generate graph of COG categories percentage.

### 09_Module_completeness_graph.py

Generate module completeness graph

### Prokka.py

run prokka tool.

---

## Input Files

The scripts require:

- Panaroo gene presence/absence matrix
- eggNOG annotation results
- presence/absence matrix for all modules
- KEGG module information
- Bacteriocin results
- neuroactive genes presence/absence matrix
- 318 genomes isolation source, accession number and strain name file

---

## Output

The scripts generate:

- Heatmaps
- Pangenome curves
- Module completeness analysis
- Statistical results
- Publication-quality figures

---

## Citation

If you use these scripts, please cite:

Asim A. et al.
Computational Comparative Genomic Analysis and Module Reconstruction of Lactiplantibacillus plantarum Reveals Its Neurometabolic Potential at the Strain Level.