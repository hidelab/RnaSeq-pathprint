# RnaSeq-pathprint
New pipeline for pathprinting RNA-seq data

## Status
Operation | Status | Note
--- | --- | ---
Experiment download | `-` | -
Experiment Fingerprinting | `-` | -
Post-processing step 1 | `-` | -
Post-processing step 2 | `-` | -
Post-processing step 3 | `-` | -
Metadata matrix creation | `-` | -
Fingerprint matrix creation | `-` | -

## Scripts

### Pre-processing
1. rna_pathprint_download.R : Script to download the RNA-seq data (gene TPMs) from Gtex portal (GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct) and split them in a number of chunks (adjustable by number_of_chunks variable). Before saving to an equal number of matrices, the Ensemble gene IDs are mapped to Entrez IDs and updated using the entrezUpdate.R script from the package GMAfunctions.

### Post-processing


## Notes
