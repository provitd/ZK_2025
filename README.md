This directory contains scripts of the bioinformatic analysis of the paper (...).

- [Description of the directory](#description-of-the-directory)
- [R and Package versions](#r-and-package-versions)
- [Explore scRNA-seq dataset interactively](#explore-scrna-seq-dataset-interactively)

# Description of the directory


The analysis is separated into different folders:
- [RNA_seq](RNA_seq/): analysis of femur and kidney RNA-seq datasets using [DESeq2](https://github.com/thelovelab/DESeq2)
- [scRNA_seq](scRNA_seq/): analysis of kidney scRNA-seq dataset using [Seurat](https://github.com/satijalab/seurat), [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) and [CellChat](https://github.com/jinworks/CellChat)
  

The files in [utils](utils/) directory contain R scripts with general functions for the analysis:
- fct_pathway.R: functions for pathway enrichment analysis based on [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) 
- sctype_wrapper.R: a copy of [sc-type](https://github.com/IanevskiAleksandr/sc-type) code saved for reproducibility reasons

# R and Package versions

The version of R used is 4.3.2 (2023-10-31) -- "Eye Holes".

The versions of the packages are described in the `renv.lock` file.
