# Master-s-Thesis
# HGSOC Spatial Transcriptomics Pipeline

Reproducible spatial transcriptomics analysis pipeline applied to high-grade serous ovarian carcinoma (HGSOC) samples, developed as part of a Master's thesis (TFM) in Bioinformatics at CIMA / Universidad de Navarra, within the PITAGORAS project.

## Overview

This repository contains the full analysis pipeline used to process and analyze four 10x Genomics Visium spatial transcriptomics samples from HGSOC patients: **GIN63, GIN65, GIN67, GIN71**. The pipeline combines Python (Jupyter notebooks, scanpy/squidpy) and R (Quarto `.qmd` scripts) and covers the full workflow from raw Space Ranger output to cell-cell communication analysis.

A key design principle of this pipeline is **reproducibility**: quality-control filtering thresholds are applied uniformly and identically across all four samples rather than being tuned per sample.

## Samples

| Sample | Visium capture area | Slide |
|--------|---------------------|-------|
| GIN63  | A1                  | V14M06-322 |
| GIN65  | B1                  | V14M06-322 |
| GIN67  | C1                  | V14M06-322 |
| GIN71  | D1                  | V14M06-322 |

Reference genome: **GRCh38-2024-A** (Space Ranger).

Pathologist annotations were derived from an adjacent H&E section (not the exact Visium section) and are used as a visual reference and spatial safeguard rather than as classification ground truth.

## Pipeline structure

| Script | Description | Output |
|--------|-------------|--------|
| `00` | Data wrangling of pathologist annotations | `GIN{XX}_patologist_annotation.h5ad` |
| `01` | QC filtering and normalization | `GIN{XX}_filtered_normalized.h5ad` |
| `02a` | Clustering (Leiden / Louvain / GraphST) | `GIN{XX}_louvain_leiden_graphst.h5ad` |
| `02b` | Clustering (BayesSpace) | `GIN{XX}_bayesspace.h5ad` |
| `03` | Tumor purity estimation (ESTIMATE) | `GIN{XX}_estimate.h5ad` |
| `04` | Copy number inference (InferCNV) | `GIN{XX}_infercnv.h5ad` |
| `05` | Cell type deconvolution (cell2location) | `GIN{XX}_cell2location.h5ad` |
| `06` | Cell-cell communication (CellChat) | Interaction networks (no h5ad output) |

Each script reads the h5ad produced by the previous step and writes a new, distinctly named h5ad file — no file is overwritten, so intermediate results remain traceable at every stage.

## Repository structure

```
.
├── Scripts/
│   ├── 00_Patologist_annotation.ipynb
│   ├── 01_Filtrado_normalizado.ipynb
│   ├── 02a_Clustering_Leiden_Louvain_GraphST.ipynb
│   ├── 02b_Clustering_BayesSpace.qmd
│   ├── 03_ESTIMATE.qmd
│   ├── 04_InferCNV.qmd
│   ├── 05_Cell2location.ipynb
│   └── 06_CellChat.qmd
├── h5ad_outputs/          # intermediate .h5ad files (not tracked, see below)
├── Figuras/
│   └── figuras_{SAMPLE_ID}/
│       ├── Filtrado_normalizado/
│       ├── Clustering/
│       └── ...
├── environment.yml        # conda environment specification
└── README.md
```

## Requirements

- Python and R, managed through a single conda environment (`spatial_analysis`)
- Key Python packages: `scanpy`, `squidpy`, `GraphST`, `cell2location` (scvi-tools)
- Key R packages: `BayesSpace`, `ESTIMATE`, `infercnv`, `CellChat`, `zellkonverter`
- GPU recommended for `cell2location` (script 05)
- Developed and tested on a local workstation (ThinkStation P5) and the UNAV HPC cluster (Slurm)

```bash
conda env create -f environment.yml
conda activate spatial_analysis
```

## Usage

Scripts are numbered and intended to be run in order, from `00` through `06`. Each notebook/`.qmd` file expects the h5ad produced by the previous step as input; paths are configured at the top of each script.

```bash
# Example: run script 01 for a given sample
jupyter nbconvert --to notebook --execute Scripts/01_Filtrado_normalizado.ipynb
```

R scripts (`.qmd`) can be rendered with Quarto:

```bash
quarto render Scripts/03_ESTIMATE.qmd
```

## Known limitations

- **Genome build check pending**: the InferCNV gene order file is labeled hg19 while the rest of the pipeline uses hg38; this needs to be resolved before any karyotype/cytoband analysis.
- **BayesSpace** runs require 50,000 MCMC iterations for stable results (increased from the package default).
- **CellChat**: `cellchat@data.signaling` currently returns `NULL` after `subsetData()` in CellChat v2.2.0.9001; under investigation.
- `TumorPurity` from ESTIMATE is a proxy value derived from a normalized `ESTIMATEScore`, not the original Affymetrix cosine-based metric (which is not valid for spatial RNA-seq data).

## Data availability

Raw sequencing data, Space Ranger outputs, and pathologist annotations are not included in this repository due to patient data privacy. Processed/intermediate `.h5ad` files are also excluded (see `.gitignore`). [Adjust this section depending on whether/where raw data will be deposited, e.g. GEO/EGA accession.]

## Authors and supervision

- **Irati Martínez** — TFM student, CIMA / Universidad de Navarra
- **Dra. Sandra Hervás Stubbs** — Primary thesis supervisor
- **Dr. Ángel Mario Martínez Montes** — Co-supervisor

## Acknowledgments

This work was developed within the PITAGORAS project at CIMA, Universidad de Navarra.

## License

[Choose a license, e.g. MIT — see step-by-step guide]
