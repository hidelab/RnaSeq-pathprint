# RnaSeq-pathprint
New pipeline for pathprinting RNA-seq data

## Status
Operation | Status | Note
--- | --- | ---
Experiment download | `Complete` | `Less experiments than expected (418) during run-time`
Experiment Fingerprinting | `Complete` | -
Post-processing step 1 | `Complete` | -
Post-processing step 2 | `Complete` | -
Post-processing step 3 | `Complete` | -
Metadata matrix creation | `Complete` | -
Fingerprint matrix creation | `-` | -

## Scripts

### Pre-processing
1. **rna_pathprint_download.R** : Script to download the RNA-seq data (gene TPMs) from Gtex portal (GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct). Then it loads the Ensembl-Entrez ID mapping table (from BioMart), filters it by keeping the pathprint genes and then the Gtex matrix is also filtered based on the mapping table. Finally, is splits the Gtex matrix in a number of chunks (adjustable by number_of_chunks variable). Before saving to an equal number of matrices. It's worth noting that duplicate ensembl IDs are removed and duplicate entrez IDs are averaged.

2. **rna_pathprint_runner.R** : Script that initiates a number of cores, each one utilizing a number of experiment expression files to run a *rna_parallel_x.R* script.

3. **rna_parallel_x.R** : Parallel scripts(14). Each initiates the experiment-fingerprinting process of multiple experiments by running the *gctx2fingerprint.R* script.

4. **gtex2fingerprint.R** : Script that accepts a experiment name, loads it's matrix and produces the experiment fingerprints. It's a merge of geo2fingerprint.R and exprs2fingerprint.R. It utilizes *custom.single.chip.enrichment.R* to produce the experiment fingerprints.

5. **custom.single.chip.enrichment.R** : Slightly tweaked version of *single.chip.enrichment.R* that omits some unrelated parameters.

### Post-processing

1. **Fingerprint_post_processing_step1_rna.R** :  Creates the *sq_.frame.xxxx-xx-xx.RData* file (SCE dataframe) and the *platform.frame.xxxx-xx-xx.RData* file (platform dataframe). It also creates the *sq_.pathway.SCE.x.RData* files for each pathway.

2. **Fingerprint_post_processing_step2_rna.R** :  Creates the POE (probability of expression) file for each pathway. 

3. **Fingerprint_post_processing_step3_rna.R** :  Creates the *sq_.matrix.xxxx-xx-xx.RData* POE matrix file (one file for all pathways). 

4. **matrix_metadata_rna.R** : Creates the *LINCS.metadata.matrix.RData* file.
 
5. **constructing POE thresholds.R** : Obtains POE threshold values and creates *LINCS.fingerprint.matrix.RData* file.

## Notes

1. We are using ensembl75 from BioMart for the ensembl-entrez ID mapping

2. **gtex2fingerprint.R** is used because the data structure is identical after the preprocessing step in **rna_pathprint_download.R**
