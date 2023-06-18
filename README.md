# Bioinformatics of short variant calling

This repository contains data and code used to prepare figures for the review on variant calling methods for medical genetics. The `bmc_review.R` script contains all the code used for figure generations. Additional data files in the repository include:

* Summary data files used for figure genration (`rma_counts.tsv`, `clinvar_counts.tsv`, `MF_assemblies.tsv`, `giab_allregions.tsv`)
* BED files of the genomic regions used in the analysis (hard-to-call CDS regions from GIAB v.4.2., BED files with GENCODE v43 CDS remapped to b37 and T2T
* Scripts used to count RMA sites (`count_rma.sh`, `liftover_t2t.sh`) and evaluate call concordance for different sets of regions (`rtg_nonmaster.sh`, `make_matrix.py`).
