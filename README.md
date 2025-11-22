# Axolotl Limb Regeneration: Single-Cell Transcriptomic Analysis

A comprehensive bioinformatics pipeline for analyzing single-cell RNA-seq data from axolotl limb regeneration, reproducing and extending the findings of Li et al. (2021).

## Project Overview

This repository contains a reproducible computational workflow for analyzing single-cell transcriptomic data from the axolotl (*Ambystoma mexicanum*) limb regeneration model. The project aims to:

- **Reproduce and extend** Li et al. (2021) using updated genome assemblies and analytical methods
- **Reconstruct neural cell trajectories** during regeneration
- **Map ligand-receptor interactions** between neural cells and immune cells

### Dataset

- **Source**: Li et al. (2021) - NCBI SRA PRJNA589484 / CNSA CNP0000706
- **Scale**: ~80,000 single cells across 8 timepoints
- **Technology**: 10X Genomics Chromium (v2 chemistry)
- **Timepoints**: 0h (control), 3h, 1d, 3d, 7d, 14d, 22d, 33d post-amputation

## Table of Contents

- [Quick Start](#quick-start)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [References](#references)

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/ehda-dolatabadi/axolotl-regeneration-scrna.git
cd axolotl-regeneration-scrna

# Run Cell Ranger alignment (Step 1)
sbatch scripts/cellranger/slurm/01_mkref.sbatch
sbatch scripts/cellranger/slurm/02_cellranger.sbatch

# Run Seurat analysis (Step 2)
bash scripts/seurat/run_pipeline.sh

# Additional analyses coming soon:
# - Monocle3 trajectory inference
# - CellChat cell-cell communication
```

---

## Pipeline Architecture

### Pipeline Components

| Component | Purpose | Key Tools | Documentation |
|-----------|---------|-----------|---------------|
| **Cell Ranger** | FASTQ → Count Matrix | `mkref`, `count`, `aggr` | [Cell Ranger README](scripts/cellranger/README.md) |
| **Seurat Analysis** | QC → Clustering → Markers | Seurat, SCTransform, Leiden | [Seurat README](scripts/seurat/README.md) |
| **Monocle3** | Trajectory Inference & Pseudotime | Monocle3, UMAP | [Monocle3 README](scripts/monocle3/README.md) *(Coming soon)* |
| **CellChat** | Cell-Cell Communication | CellChat | [CellChat README](scripts/cellchat/README.md) *(Coming soon)* |

---

## Prerequisites

**Software Requirements:**
- Cell Ranger
- R with packages: Seurat, dplyr, ggplot2, future
- Monocle3 (for trajectory analysis - coming soon)
- CellChat (for communication analysis - coming soon)
- SLURM job scheduler

---

## Usage

### Full Pipeline Execution

**Step 1: Build Reference Genomes**
```bash
sbatch scripts/cellranger/slurm/01_mkref.sbatch
```

**Step 2: Process FASTQ Files**
```bash
sbatch scripts/cellranger/slurm/02_cellranger.sbatch
```

**Step 3: Run Seurat Analysis**
```bash
bash scripts/seurat/run_pipeline.sh
```

**Step 4: Trajectory Inference with Monocle3** *(Coming soon)*
```bash
bash scripts/monocle3/run_monocle3.sh
```

**Step 5: Cell-Cell Communication with CellChat** *(Coming soon)*
```bash
bash scripts/cellchat/run_cellchat.sh
```

### Pipeline Stages

#### Stage 1: Cell Ranger (Alignment & Quantification)
- **Input**: Raw FASTQ files from 10X Genomics sequencing
- **Output**: Gene-barcode count matrices (`filtered_feature_bc_matrix/`)
- **Time**: ~15-20 hours per sample

#### Stage 2: Seurat (scRNA-seq Analysis)
- **Input**: Count matrices from Cell Ranger
- **Output**: Clustered cells with annotations, marker genes, QC reports
- **Time**: ~2-3 hours for full dataset

#### Stage 3: Monocle3 (Trajectory Inference) *(Coming soon)*
- **Input**: Annotated Seurat objects
- **Output**: Pseudotime trajectories, developmental paths, gene dynamics
- **Time**: TBD

#### Stage 4: CellChat (Cell-Cell Communication) *(Coming soon)*
- **Input**: Annotated Seurat objects
- **Output**: Ligand-receptor networks, signaling pathways, communication patterns
- **Time**: TBD

---

## Project Structure

```
.
├── scripts/
│   ├── cellranger/           # Cell Ranger alignment pipeline
│   │   ├── slurm/            # SLURM batch scripts
│   │   │   ├── 01_mkref.sbatch
│   │   │   ├── 02_cellranger.sbatch
│   │   │   ├── 03_aggr.sbatch
│   │   │   └── 03_aggr_normalize.sbatch
│   │   └── README.md         # Cell Ranger documentation
│   │
│   ├── seurat/               # Seurat scRNA-seq analysis pipeline
│   │   ├── 00_map_features.R      # Gene ID mapping
│   │   ├── 01_qc.R                # Quality control
│   │   ├── 02_filter.R            # Cell filtering
│   │   ├── 03_normalize.R         # SCTransform normalization
│   │   ├── 04_integrate.R         # Batch correction
│   │   ├── 05_pca.R               # PCA
│   │   ├── 06_cluster.R           # Clustering
│   │   ├── 07_umap.R              # UMAP
│   │   ├── 08_find_markers.R      # Marker identification
│   │   ├── 09_score_markers.R     # Cell type scoring
│   │   ├── 10_subset_clusters.R   # Cluster subsetting
│   │   ├── 11_subset_cells.R      # Cell type filtering
│   │   ├── run_pipeline.sh        # Pipeline orchestration
│   │   ├── slurm/                 # SLURM batch scripts
│   │   └── README.md              # Seurat documentation
│   │
│   ├── monocle3/                  # Trajectory inference (Coming soon)
│   │   └── README.md              # Monocle3 documentation
│   │
│   └── cellchat/             	   # Cell-cell communication (Coming soon)
│       └── README.md              # CellChat documentation
│
├── data/                     # Input data (not tracked)
│   ├── Li_dataset/           # Raw sequencing data
│   └── ref_genomes/	      # Reference genomes and annotations
│
├── results/                  # Analysis outputs (not tracked)
│   ├── cellranger/           # Alignment results
│   ├── seurat/               # Analysis results
│   ├── monocle3/             # Trajectory analysis results
│   └── cellchat/             # Cell communication results
│
└── README.md                 # This file
```

---

## Reference Genomes

This project uses two axolotl reference genomes:

| Genome | Source | Version | Notes |
|--------|--------|---------|-------|
| **UKY_AmexF1_1** | NCBI | GCF_040938575.1 | Primary analysis genome (2024) |
| **AmexG_v6.0-DD** | Axolotl-omics | v6.0 / v47 transcriptome | Comparison genome |

### Download Links

**UKY_AmexF1_1 (NCBI)**
- Genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/938/575/GCF_040938575.1_UKY_AmexF1_1/GCF_040938575.1_UKY_AmexF1_1_genomic.fna.gz
- Annotation: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/938/575/GCF_040938575.1_UKY_AmexF1_1/GCF_040938575.1_UKY_AmexF1_1_genomic.gtf.gz

**AmexG_v6.0-DD (Axolotl-omics)**
- Genome: https://www.axolotl-omics.org/dl/AmexG_v6.0-DD.fa.gz
- Annotation: https://www.axolotl-omics.org/dl/AmexT_v47-AmexG_v6.0-DD.gtf.gz

---

## References

### Primary Publication
- **Li, H., Wei, X., Zhou, L., Zhang, W., Wang, C., Guo, Y., Li, D., Chen, J., Liu, T., Zhang, Y., Ma, S., Wang, C., Tan, F., Xu, J., Liu, Y., Yuan, Y., Chen, L., Wang, Q., Qu, J., Shen, Y., et al. (2021).** Dynamic cell transition and immune response landscapes of axolotl limb regeneration revealed by single-cell analysis. *Protein & Cell, 12*(1), 57-66. https://doi.org/10.1007/s13238-020-00763-1

### Related Axolotl Regeneration Studies
- **Gerber, T., Murawala, P., Knapp, D., Masselink, W., Schuez, M., Hermann, S., Gac-Santel, M., Nowoshilow, S., Kageyama, J., Khattak, S., Currie, J. D., Camp, J. G., Tanaka, E. M., & Treutlein, B. (2018).** Single-cell analysis uncovers convergence of cell identities during axolotl limb regeneration. *Science, 362*(6413), eaaq0681. https://doi.org/10.1126/science.aaq0681

- **Leigh, N. D., Dunlap, G. S., Johnson, K., Mariano, R., Oshiro, R., Wong, A. Y., Bryant, D. M., Miller, B. M., Ratner, A., Chen, A., Ye, W. W., Haas, B. J., & Whited, J. L. (2018).** Transcriptomic landscape of the blastema niche in regenerating adult axolotl limbs at single-cell resolution. *Nature Communications, 9*(1), 5153. https://doi.org/10.1038/s41467-018-07604-0

- **Rodgers, A. K., Smith, J. J., & Voss, S. R. (2020).** Identification of immune and non-immune cells in regenerating axolotl limbs by single-cell sequencing. *Experimental Cell Research, 394*(2), 112149. https://doi.org/10.1016/j.yexcr.2020.112149

### Regeneration Biology Background
- **McCusker, C., Bryant, S. V., & Gardiner, D. M. (2015).** The axolotl limb blastema: cellular and molecular mechanisms driving blastema formation and limb regeneration in tetrapods. *Regeneration, 2*(2), 54-71. https://doi.org/10.1002/reg2.32

- **Farkas, J. E., & Monaghan, J. R. (2017).** A brief history of the study of nerve dependent regeneration. *Neurogenesis, 4*(1), e1302216. https://doi.org/10.1080/23262133.2017.1302216

### Schwann Cell Biology
- **Liu, Z., Jin, Y.Q., Chen, L., Wang, Y., Yang, X., Cheng, J., Wu, W., Qi, Z., & Shen, Z. (2015).** Specific marker expression and cell state of Schwann cells during culture in vitro. *PLOS ONE, 10*(4), e0123278. https://doi.org/10.1371/journal.pone.0123278

- **Wolbert, J., Li, X., Heming, M., Mausberg, A.K., Akkermann, D., Frydrychowicz, C., Fledrich, R., Stettner, M., Alonso Gonzalez, N., Wiendl, H., Meyer zu Hörste, G., & Stassart, R.M. (2020).** Redefining the heterogeneity of peripheral nerve cells in health and autoimmunity. *Proceedings of the National Academy of Sciences, 117*(17), 9466-9476. https://doi.org/10.1073/pnas.1912139117

### Axolotl Genomics
- **NCBI. (2024).** Ambystoma mexicanum genome assembly UKY_AmexF1_1 (GCF_040938575.1). NCBI Datasets. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040938575.1/

- **Axolotl-omics. (n.d.).** Assemblies. https://www.axolotl-omics.org/assemblies

### Software and Methods
- **Cell Ranger**: 10X Genomics. (2023). Cell Ranger Software. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

- **10X Genomics Chromium**: Zheng, G.X.Y., Terry, J.M., Belgrader, P., Ryvkin, P., Bent, Z.W., Wilson, R., Ziraldo, S.B., Wheeler, T.D., McDermott, G.P., Zhu, J., Grego>

- **Seurat**: Hao, Y., Stuart, T., Kowalski, M.H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology, 42*, 293-304. https://doi.org/10.1038/s41587-023-01767-y

- **SCTransform**: Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology, 20*, 296. https://doi.org/10.1186/s13059-019-1874-1

- **future package**: Bengtsson, H. (2021). A unifying framework for parallel and distributed processing in R using futures. *The R Journal, 13*(2), 208-227. https://doi.org/10.32614/RJ-2021-048

- **Monocle3**: Cao, J., Spielmann, M., Qiu, X., Huang, X., Ibrahim, D.M., Hill, A.J., Zhang, F., Mundlos, S., Christiansen, L., Steemers, F.J., Trapnell, C., & Shendure, J. (2019). The single-cell transcriptional landscape of mammalian organogenesis. *Nature, 566*(7745), 496-502. https://doi.org/10.1038/s41586-019-0969-x

- **CellChat**: Jin, S., Guerrero-Juarez, C.F., Zhang, L., Chang, I., Ramos, R., Kuan, C.H., Myung, P., Plikus, M.V., & Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications, 12*(1), 1088. https://doi.org/10.1038/s41467-021-21246-9

### Cell Type Annotation Tools
- **webCSEA**: Cell-type Specific Enrichment Analysis of Genes. UTHealth School of Biomedical Informatics. https://bioinfo.uth.edu/webcsea/

### Tutorials and Documentation
- **Satija Lab. (2023).** Seurat - Guided Clustering Tutorial. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

- **Babraham Bioinformatics. (n.d.).** Seurat workflow for 10X single-cell RNA-seq analysis. https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html

### AI Tools
- **Anthropic. (2025).** Claude (4.5 Sonnet). https://claude.ai/

*Note: AI tools were used for gene annotation assistance and text refinement. All content has been reviewed by the author to ensure accuracy and originality.*

## Contact

**Ehda Dolatabadi**
Master of Science in Bioinformatics, Northeastern University
Email: [dolatabadi.e@northeastern.edu](mailto:dolatabadi.e@northeastern.edu)

---

**Last Updated**: November 20, 2025
**Status**: Active Development
**Version**: 1.0.0
