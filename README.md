# Bioinformatics Capstone

This repository is for the capstone project as a part of my Master’s degree in Bioinformatics at Northeastern University.
The project reproduces and critically re-analyzes Li et al. (2021): *Dynamic cell transition and immune response landscapes of axolotl limb regeneration revealed by single-cell analysis*, using the axolotl (Ambystoma mexicanum) as the model organism for regeneration studies.

## Data

- Single-cell RNA-seq dataset from Li et al. (2021).
- Public scRNA-seq dataset (NCBI SRA accession PRJNA589484 / CNSA CNP0000706).
- ~41,000 single cells across regeneration timepoints.
- Eight regeneration timepoints: 0h, 3h, 1d, 3d, 7d, 14d, 21d, and 33d post-amputation.
- Sequenced using 10X Genomics v2 chemistry.

## Project Possible Directions

- Address gaps in the current axolotl reference genome and many missing annotations.  
- Consider the newer genome assembly for improved analysis.
- Perform trajectory analysis for connective tissues or other cell populations.  
- Examine Schwann cell responses to event over time.
- Map fibroblast state changes across regeneration. 
- Explore additional signaling responses of interest. 
- Investigate how the immune system interacts with the trauma event.
- Investigate the large migration of macrophages at the timepoint.

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

- Li, H., Wei, X., Zhou, L., Zhang, W., Wang, C., Guo, Y., Li, D., Chen, J., Liu, T., Zhang, Y., Ma, S., Wang, C., Tan, F., Xu, J., Liu, Y., Yuan, Y., Chen, L., Wang, Q., Qu, J., Shen, Y., … Xu, X. (2021). Dynamic cell transition and immune response landscapes of axolotl limb regeneration revealed by single-cell analysis. *Protein & cell, 12*(1), 57–66. https://doi.org/10.1007/s13238-020-00763-1
- McCusker, C., Bryant, S. V., & Gardiner, D. M. (2015). The axolotl limb blastema: cellular and molecular mechanisms driving blastema formation and limb regeneration in tetrapods. *Regeneration (Oxford, England), 2*(2), 54–71. https://doi.org/10.1002/reg2.32
