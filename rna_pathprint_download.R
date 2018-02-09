# Script's Job: Get the .gtex matrix and start threads which create fingerprint 
# files per experiment 
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 29.1.2018
# Script to download the RNA-seq data (gene TPMs) from Gtex portal 
# (GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct) and split them in a 
# number of chunks (adjustable by number_of_chunks variable). Before saving to 
# an equal number of matrices, the Ensemble gene IDs are mapped to Entrez IDs 
# and updated using the entrezUpdate.R script from the package GMAfunctions. 
# Finally, the script initiates a number of cores each one handling running the
# gctx2fingerprint.R script to produce the experiment fingerprints.

# For this script we need:
# > (optional) The .gctx matrix downloaded here "~/LINCS project/LINCS specific scripts/gctx_Folder"

# == Load ready things ==
expr <- readRDS("~/RNA-seq pathprint/pathprint_rnaseq/rnaseq_Folder/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.RDS")
pathprint.Hs.gs <- readRDS("~/RNA-seq pathprint/pathprint_rnaseq/rnaseq_Folder/pathprint.Hs.gs20170427.RDS")
ENSEMBL_ENTREZ_mapping <- readRDS("~/RNA-seq pathprint/pathprint_rnaseq/rnaseq_Folder/ensembl75toentrezgene_for_RNAseq_pathprint.RDS")
# == End loading ==


# Script start
# load library tha handles .GTEx files
library(CePa)
library(GMAfunctions)
library(pathprint)
library(org.Hs.eg.db)

# Select the number of chunks
number_of_chunks <- 15


# Check if the GTEx matrix is in the rnaseq_Folder. If not download manually and add it.
if(!file.exists("~/RNA-seq pathprint/pathprint_rnaseq/rnaseq_Folder/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct",showWarnings = TRUE)[1])
    stop("Cannot find .GTEx file.")

# Read the  gene tpm data and save to .RDS format
expr<-read.gct("rnaseq_Folder/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")
saveRDS(expr, file = "rnaseq_Folder/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.RDS")

load("EntrezHistoryList.RData")
source("entrezGeneIDUpdate.R")

# ---------------------------------------------------------------------------- #
# ------ Create updated pathprint gene list pathprint.Hs.gs20170427.RDS ------ #
# ---------------------------------------------------------------------------- #
data("pathprint.Hs.gs")

preIDs<-unique(unlist(pathprint.Hs.gs))
print(paste("Entrez gene IDs in pathprint: ", length(preIDs), sep = ""))
# Entrez gene IDs in pathprint: 10903

# Update with recent Entrez IDs
pathprint.Hs.gs<-lapply(pathprint.Hs.gs, entrezGeneIDUpdate, 
                        EntrezHistoryList=EntrezHistoryList)

# Get the updated entrez IDs in pathprint
pathprint.Hs.gs<-lapply(pathprint.Hs.gs, setdiff, y="-")
postIDs<-unique(unlist(pathprint.Hs.gs))
print(paste("Updated Entrez gene IDs in pathprint: ", length(postIDs), sep = ""))

saveRDS(pathprint.Hs.gs, file="rnaseq_Folder/pathprint.Hs.gs20170427.RDS")
# Updated Entrez gene IDs in pathprint: 10819
# ---------------------------------------------------------------------------- #
# -------------------- End of updated pathprint gene list -------------------- #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# ----------------------- Create ENSEMBL_ENTREZ_mapping ---------------------- #
# ---------------------------------------------------------------------------- #

# Ensembl gene IDs
ensgs=rownames(expr)
ensgs<-sub("[[:punct:]][[:digit:]]*", "", ensgs)


library("biomaRt")
listMarts(host = 'Feb2014.archive.ensembl.org')
ensembl75 <- useMart(host='Feb2014.archive.ensembl.org',
                     biomart='ENSEMBL_MART_ENSEMBL',
                     dataset='hsapiens_gene_ensembl')
result1<-getBM(attributes=c('ensembl_gene_id', 'entrezgene'),
               filters = 'ensembl_gene_id',
               values = ensgs,
               mart = ensembl75)

length(intersect(result1$entrezgene, unique(unlist(pathprint.Hs.gs))))
# Number of Entrez IDs in pathprint gene IDs: 10583

result1$entrezgene<-entrezGeneIDUpdate(result1$entrezgene, EntrezHistoryList)
ENSEMBL_ENTREZ_mapping <- result1

# Get the Ensembl to Entrez IDs table (alternative way)
# ENSEMBL2EG<-toTable(org.Hs.egENSEMBL2EG)

# Number of Ensembl IDs that match Ensembl IDs in the table: 10615/11059
length(intersect(ensgs, ENSEMBL_ENTREZ_mapping$ensembl_gene_id))
nrow(ENSEMBL_ENTREZ_mapping)

# Dublicate Entrez and Ensembl IDs
dent <- sum(duplicated(ENSEMBL_ENTREZ_mapping$entrezgene))    # 441
dens <- sum(duplicated(ENSEMBL_ENTREZ_mapping$ensembl_gene_id)) # 502

# Limit to ENSGs measured by GTEx
ENSEMBL_ENTREZ_mapping<-subset(ENSEMBL_ENTREZ_mapping, subset=(entrezgene %in% unique(unlist(pathprint.Hs.gs)) ))
length(unique(ENSEMBL_ENTREZ_mapping$entrezgene)) #10572
length(unique(ENSEMBL_ENTREZ_mapping$ensembl_gene_id)) #10615
sum(duplicated(ENSEMBL_ENTREZ_mapping)) #0
ENSEMBL_ENTREZ_mapping<-ENSEMBL_ENTREZ_mapping[!duplicated(ENSEMBL_ENTREZ_mapping),]
length(unique(ENSEMBL_ENTREZ_mapping$entrezgene)) #10572
length(unique(ENSEMBL_ENTREZ_mapping$ensembl_gene_id)) #10615

# Update entrez IDs
ENSEMBL_ENTREZ_mapping$entrezgene<-entrezGeneIDUpdate(ENSEMBL_ENTREZ_mapping$entrezgene, EntrezHistoryList)
length(intersect(unique(unlist(pathprint.Hs.gs)), unique(ENSEMBL_ENTREZ_mapping$entrezgene))) #10569

saveRDS(ENSEMBL_ENTREZ_mapping, file="rnaseq_Folder/ensembl75toentrezgene_for_RNAseq_pathprint.RDS")

# ---------------------------------------------------------------------------- #
# ----------------------- End of ENSEMBL_ENTREZ_mapping ---------------------- #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# ------------------ Create mapped [genes x samples] matrix ------------------ #
# ---------------------------------------------------------------------------- #

dim(expr) # GTEx dimensions ensembl IDs:56202 and samples:11688

# After that we must be careful to only keep the genes we have in both IDs
# Removing .x as this is a vrsion number and not stable
rownames(expr) <- sub("[[:punct:]][[:digit:]]*", "", rownames(expr))

# Removing duplicate ensembl IDs (441 duplicate Ensembl IDs removed)
ENSEMBL_ENTREZ_mapping = ENSEMBL_ENTREZ_mapping[!duplicated(ENSEMBL_ENTREZ_mapping$ensembl_gene_id),]

active_ensembl_ids <- ENSEMBL_ENTREZ_mapping$ensembl_gene_id
active_entrez_ids <- ENSEMBL_ENTREZ_mapping$entrezgene

active_expr <- subset(expr, rownames(expr) %in% active_ensembl_ids)
new_rownames <- rownames(active_expr)

dim(active_expr) # new GTEx dimensions ensembl IDs:10645 and samples:11688

# Convert Ensembl to Entrez IDs (Entrez dublicates still in)
new_rownames

# Get the entrez ID for any ensembl ID
mapIDs <- function(ensg_code) {
    ENSEMBL_ENTREZ_mapping[2]$entrezgene[which(grepl(ensg_code, ENSEMBL_ENTREZ_mapping$ensembl_gene_id))]
}

# Actual convertion
entrez_row<-lapply(new_rownames, mapIDs)
rownames(active_expr) <- entrez_row

saveRDS(entrez_expr, file="rnaseq_Folder/entrez_expr.RDS") # 157 duplicated Entrez IDs

# Average the rows that have duplicate entrez IDs
library(data.table)

# 1. Create first column with the rownames and aggregate based on that 
temp_mat <- cbind(rownames(entrez_expr),entrez_expr)
colnames(temp_mat)[1] <- "Genes"
dim(temp_mat) # 10645 11689

# 2. Average duplicate rows
# Find duplicates
entrez_expr_numeric<-sapply(entrez_expr, as.numeric)
dups<-rownames(entrez_expr)[duplicated(rownames(entrez_expr))]

# Frst this
while(TRUE){
a<-which(rownames(entrez_expr) == dups[1]) 

if(length(a) ==2 ) {
    #2
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), as.numeric(entrez_expr[a[2], ])))
    entrez_expr <- entrez_expr[-c(a[2]), ]
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")
}
if (length(a) ==3 ) {
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), as.numeric(entrez_expr[a[2], ]), as.numeric(entrez_expr[a[3], ])))
    entrez_expr <- entrez_expr[-c(a[2],a[3]), ]
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")    
}
if(length(a) ==1) {
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")
}
if (length(a) ==4 ) {
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), as.numeric(entrez_expr[a[2], ]), as.numeric(entrez_expr[a[3], ]), as.numeric(entrez_expr[a[4], ])))
    entrez_expr <- entrez_expr[-c(a[2],a[3],a[4]), ]
    dims<-dim(entrez_expr)
    dups <- dups[-(1)]
    print("DONE")
}
if (length(a) ==5 ) {
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), as.numeric(entrez_expr[a[2], ]), as.numeric(entrez_expr[a[3], ]), as.numeric(entrez_expr[a[4], ]), as.numeric(entrez_expr[a[5], ])))
    entrez_expr <- entrez_expr[-c(a[2],a[3],a[4],a[5]), ]
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")
}
if (length(a) ==6 ) {
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), as.numeric(entrez_expr[a[2], ]), as.numeric(entrez_expr[a[3], ]), as.numeric(entrez_expr[a[4], ]), 
                                          as.numeric(entrez_expr[a[5], ]), as.numeric(entrez_expr[a[6], ])))
    entrez_expr <- entrez_expr[-c(a[2],a[3],a[4],a[5],a[6]), ]
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")
}
if (length(a) ==14 ) {
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), 
                                          as.numeric(entrez_expr[a[2], ]), 
                                          as.numeric(entrez_expr[a[3], ]), 
                                          as.numeric(entrez_expr[a[4], ]), 
                                          as.numeric(entrez_expr[a[5], ]),
                                          as.numeric(entrez_expr[a[6], ]), 
                                          as.numeric(entrez_expr[a[7], ]), 
                                          as.numeric(entrez_expr[a[8], ]), 
                                          as.numeric(entrez_expr[a[9], ]),
                                          as.numeric(entrez_expr[a[10], ]), 
                                          as.numeric(entrez_expr[a[11], ]), 
                                          as.numeric(entrez_expr[a[12], ]), 
                                          as.numeric(entrez_expr[a[13], ]),
                                          as.numeric(entrez_expr[a[14], ])
                                          ))
    entrez_expr <- entrez_expr[-c(a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13],a[14]), ]
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")
}
if (length(a) ==10 ) {
    entrez_expr[a[1], ] <- colMeans(rbind(as.numeric(entrez_expr[a[1], ]), 
                                          as.numeric(entrez_expr[a[2], ]), 
                                          as.numeric(entrez_expr[a[3], ]), 
                                          as.numeric(entrez_expr[a[4], ]), 
                                          as.numeric(entrez_expr[a[5], ]),
                                          as.numeric(entrez_expr[a[6], ]), 
                                          as.numeric(entrez_expr[a[7], ]), 
                                          as.numeric(entrez_expr[a[8], ]), 
                                          as.numeric(entrez_expr[a[9], ]),
                                          as.numeric(entrez_expr[a[10], ])
    ))
    entrez_expr <- entrez_expr[-c(a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10]), ]
    dups <- dups[-(1)]
    dims<-dim(entrez_expr)
    print("DONE")
}
}

saveRDS(entrez_expr, file="rnaseq_Folder/entrez_expr_final.RDS")
# dim(entrez_expr) 10488 11688

# ---------------------------------------------------------------------------- #
# ------------------ End of mapped [genes x samples] matrix ------------------ #
# ---------------------------------------------------------------------------- #

# column/experiment and gene count
experiments_count <- dim(entrez_expr_final)[2] # 11688
gene_count <- dim(entrez_expr_final)[1] # 10488

# The experiment files are here
gtex_exp_path <- "~/RNA-seq pathprint/pathprint_rnaseq/gtex_experiments/"

# Locate the Entrez history file to be used for updating
historyFile <- "~/RNA-seq pathprint/pathprint_rnaseq/Entrez Gene History/EntrezHistoryList.RData"

# We have to split the matrix column-wise/experiment-wise since its 1.5 Gb altogether
chunk <- as.integer(experiments_count / number_of_chunks)

# gtex matrices are located here
exp_path <- "~/RNA-seq pathprint/pathprint_rnaseq/rnaseq_Folder/gtex_matrix_"

sampleset <- vector("list", number_of_chunks)
sampleset_count <- 1

for (matrix_chunk in seq(1,experiments_count, by = chunk)){
    current_start <- matrix_chunk
    
    # Check if its last loop
    if (matrix_chunk + chunk > experiments_count) current_end <- experiments_count
    else current_end <- matrix_chunk + chunk - 1
    
    partial_gtex_matrix <- entrez_expr_final[,current_start:current_end]

    # Updating the Entrez IDs
    rownames(partial_gtex_matrix)<-entrezUpdate(rownames(partial_gtex_matrix), historyFile = historyFile)
    print("Entrez IDs updated.")
    
    # Fill the sample set
    sampleset[[sampleset_count]] <- colnames(partial_gtex_matrix)
    
    save_name <- paste(exp_path,current_start,"-",current_end,".RData", sep="")
    save(partial_gtex_matrix, file = save_name)
    print(paste("Downloaded and saved experiments ", current_start, " to ", current_end, " in ", save_name))
    sampleset_count <- sampleset_count + 1
    
    for (e in 1:length(colnames(partial_gtex_matrix))) {
        
        expression_column <- matrix(partial_gtex_matrix[,e])
        gene_column <- matrix(rownames(partial_gtex_matrix))
        exp_matrix <- cbind(gene_column,expression_column)
        
        experiment_name <- colnames(partial_gtex_matrix)[e]
        experiment_name <- gsub(":", "_", experiment_name)
        full_experiment_name <- paste(gtex_exp_path,experiment_name,".RData", sep = "")
        save(exp_matrix, file = full_experiment_name)
    }
    
}

save(sampleset, file = "sampleset.RData")
print("Sample set saved.")
