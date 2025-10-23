# Bioinformatics Capstone

This repository is for the capstone project as a part of my Master’s degree in Bioinformatics at Northeastern University.
The project reproduces and critically re-analyzes Li et al. (2021): *Dynamic cell transition and immune response landscapes of axolotl limb regeneration revealed by single-cell analysis*, using the axolotl (*Ambystoma mexicanum*) as the model organism for regeneration studies.
The aim is to use single-cell transcriptomic data to reconstruct Schwann cell trajectories and map their ligand-receptor interactions with immune cells across the regenerative timeline.

## Data

- Single-cell RNA-seq dataset from Li et al. (2021).
- Public scRNA-seq dataset (NCBI SRA accession PRJNA589484 / CNSA CNP0000706).
- ~41,000 single cells across regeneration timepoints.
- Eight regeneration timepoints: 0h, 3h, 1d, 3d, 7d, 14d, 21d, and 33d post-amputation.
- Sequenced using 10X Genomics v2 chemistry.

## Reference Genomes and Annotations

- https://www.axolotl-omics.org/dl/AmexG_v6.0-DD.fa.gz
- https://www.axolotl-omics.org/dl/AmexT_v47-AmexG_v6.0-DD.gtf.gz
- https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/938/575/GCF_040938575.1_UKY_AmexF1_1/GCF_040938575.1_UKY_AmexF1_1_genomic.fna.gz
- https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/938/575/GCF_040938575.1_UKY_AmexF1_1/GCF_040938575.1_UKY_AmexF1_1_genomic.gtf.gz

## Project Directions

- Reanalyze Li et al. (2021) scRNA-seq dataset using the updated genome (UKY_AmexF1_1)
- Implement a reproducible HPC workflow with Cell Ranger, Seurat, Monocle3, and CellChat
- Improve Schwann-cell identification and annotation using newer datasets
- Map Schwann-cell trajectories and neuro-immune ligand–receptor interactions across time points
- Interpret how Schwann-immune signaling contributes to axolotl limb regeneration

## Technical Skills Inventory

- Bioinformatics: scRNA-seq preprocessing, clustering, trajectory inference.
- Workflow management: Snakemake, Nextflow
- Environments: Conda, Docker, Singularity
- Programming: Python, R, SQL, Bash

## Learning Goal

- Apply computational methods to address a biological research question
- Gain skills in research methodologies and scientific communication
- Build experience in reproducing published analyses

## References

- Axolotl‐omics. (n.d.). *Assemblies.* https://www.axolotl-omics.org/assemblies
- Farkas, J. E., & Monaghan, J. R. (2017). A brief history of the study of nerve dependent regeneration. *Neurogenesis, 4*(1), e1302216. https://doi.org/10.1080/23262133.2017.1302216
- Gerber, T., Murawala, P., Knapp, D., Masselink, W., Schuez, M., Hermann, S., Gac-Santel, M., Nowoshilow, S., Kageyama, J., Khattak, S., Currie, J. D., Camp, J. G., Tanaka, E. M., & Treutlein, B. (2018). Single-cell analysis uncovers convergence of cell identities during axolotl limb regeneration. *Science (New York, N.Y.), 362*(6413), eaaq0681. https://doi.org/10.1126/science.aaq0681
- Leigh, N. D., Dunlap, G. S., Johnson, K., Mariano, R., Oshiro, R., Wong, A. Y., Bryant, D. M., Miller, B. M., Ratner, A., Chen, A., Ye, W. W., Haas, B. J., & Whited, J. L. (2018). Transcriptomic landscape of the blastema niche in regenerating adult axolotl limbs at single-cell resolution. *Nature Communications, 9*(1), 5153. https://doi.org/10.1038/s41467-018-07604-0
- Li, H., Wei, X., Zhou, L., Zhang, W., Wang, C., Guo, Y., Li, D., Chen, J., Liu, T., Zhang, Y., Ma, S., Wang, C., Tan, F., Xu, J., Liu, Y., Yuan, Y., Chen, L., Wang, Q., Qu, J., Shen, Y., … Xu, X. (2021). Dynamic cell transition and immune response landscapes of axolotl limb regeneration revealed by single-cell analysis. *Protein & cell, 12*(1), 57–66. https://doi.org/10.1007/s13238-020-00763-1
- McCusker, C., Bryant, S. V., & Gardiner, D. M. (2015). The axolotl limb blastema: cellular and molecular mechanisms driving blastema formation and limb regeneration in tetrapods. *Regeneration (Oxford, England), 2*(2), 54–71. https://doi.org/10.1002/reg2.32
- NCBI. (2024). Ambystoma mexicanum genome assembly UKY_AmexF1_1 (GCF_040938575.1) [*Genome assembly*]. NCBI Datasets. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040938575.1/
- Qin, T., Han, J., Fan, C., Sun, H., Rauf, N., Wang, T., Yin, Z., & Chen, X. (2024). Unveiling axolotl transcriptome for tissue regeneration with high-resolution annotation via long-read sequencing. *Computational and structural biotechnology journal, 23*, 3186–3198. https://doi.org/10.1016/j.csbj.2024.08.014
- Rodgers, A. K., Smith, J. J., & Voss, S. R. (2020). Identification of immune and non-immune cells in regenerating axolotl limbs by single-cell sequencing. *Experimental cell research, 394*(2), 112149. https://doi.org/10.1016/j.yexcr.2020.112149
