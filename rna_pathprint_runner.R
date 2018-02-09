# Script's Job: Initiates the scripts that produce experiment-fingerprints
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 07.11.2017
# Script that initiates a number of cores, each one utilizing a number of 
# experiments to run the gctx2fingerprint.R script and produce the 
# experiment-fingerprints.

# Script start
# load libraries
library(GMAfunctions)
library(parallel)
library(foreach)
library(doMC)

# Go to the parallel path
setwd("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/rna_parallel")

# Select the number of core
number_of_cores <- 14
registerDoMC(cores = number_of_cores)

# Load the sample set
load("/shared/hidelab2/shared/Sokratis/pathprint_lincs/sampleset.RData")

# The experiment, Entrez history and gtex matrices file paths
gtex_exp_path <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/gtex_experiments/"
history_file_path <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/EntrezHistoryList.RData"
matrices_path <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/rnaseq_Folder/"

# These scripts run for single chip enrichment

a = mclapply(source("rna_parallel_1.R"))
b = mclapply(source("rna_parallel_2.R"))
c = mclapply(source("rna_parallel_3.R"))
d = mclapply(source("rna_parallel_4.R"))
e = mclapply(source("rna_parallel_5.R"))
f = mclapply(source("rna_parallel_6.R"))
g = mclapply(source("rna_parallel_7.R"))
h = mclapply(source("rna_parallel_8.R"))
i = mclapply(source("rna_parallel_9.R"))
j = mclapply(source("rna_parallel_10.R"))
k = mclapply(source("rna_parallel_11.R"))
l = mclapply(source("rna_parallel_12.R"))
m = mclapply(source("rna_parallel_13.R"))
n = mclapply(source("rna_parallel_14.R"))

# Collecting processes/threads
print("Collecting threads...")
Sys.sleep(10)
x=do.call(rbind,list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))
print("Finished creating experiment fingerprints successfully!")
