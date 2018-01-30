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
1. **rna_pathprint_download.R** : Script to download the RNA-seq data (gene TPMs) from Gtex portal (GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct). Then it loads the Ensembl-Entrez ID mapping table (from org.Hs.eg.db), filters it by keeping the pathprint genes and then the Gtex matrix is also filtered based on the mapping table. Finally, is splits the Gtex matrix in a number of chunks (adjustable by number_of_chunks variable). Before saving to an equal number of matrices.

### Post-processing


## Notes
