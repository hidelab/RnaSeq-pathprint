# Script's Job: Produces experiment-fingerprints
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 08.11.2017
# Script that accepts a experiment name, loads it's matrix and produces the
# experiment fingerprints. It's a merge of 
# geo2fingerprint-exprs2fingerprint.R. It uses custom.single.chip.enrichment.R 
# to produce the experiment fingerprints.

# Load geneset related files
load("/shared/hidelab2/shared/Sokratis/pathprint_lincs/genesets.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_lincs/pathprint.Hs.gs.rda")


# Load the custom SCE script
source("/shared/hidelab2/shared/Sokratis/pathprint_lincs/custom.single.chip.enrichment.R")

# This is a function to produce a normalized fingerprint directly from an expression 
# matrix. EnrichmentMethod - algorithm used to calculate gene enrichment score and
# it can only be  single.chip.enrichment (SCE). Transformation, statistic, 
# normalizedScore and progressBar are relevant to SCE.

# EXP              : experiment name {string}
# EXPath           : path to the experiment files {string}
# GEOthreshold     : 
# geneset          : geneset/pathway name {string}
# enrichmentMethod : the enrichment method to be used {string}
# transformation   : the transformation to be used {string}
# statistic        : the statistic to be used {string}
# progressBar      : show progress bar for not {boolean}


# NOTE: before we run the whole pipeline hardcode the platform name here!!!

gctx2fingerprint<-function(	EXP, EXPath = tempdir(), GEOthreshold = TRUE,
							geneset = "KEGG and Wikipathways and static",
							enrichmentMethod = "SCE", transformation = "squared.rank",
							statistic = "mean",
							progressBar = FALSE
							)
	{
    # Hardcoded platform that corresponds to the GSE we use: GSE70138
    platform <- "GPL20573"
    
    # Check if the experiment matrix is in the gctx_experiments folder. If not 
    # stop fingerprinting execution and inform the user, otherwise load it.
    file_path <- paste(EXPath, EXP, ".RData", sep="")
    file_path_corrected <- gsub(":", "_", file_path)
    if(!file.exists(file_path_corrected, showWarnings = TRUE)[1]) {
        stop(paste("Cannot find experiment file: "), EXP, sep = "")
    }else {
	load_name <- gsub(":", "_", file_path)     
	load(load_name)
    }
    exp_matrix<-as.data.frame(exp_matrix)
    colnames(exp_matrix)[1] <- "EntrezID"
    
    # Convert to numeric as occasionally factors can complicate matters
    exp_matrix[,-1] <- apply(exp_matrix[,-1, drop = FALSE], 2, function(x){as.numeric(as.character(x))})
    
    # Take care of the duplicate gene entrez IDs (which occured after updating them)
    exp_matrix_unique <- aggregate(exp_matrix[,-1], list(EntrezID = exp_matrix$EntrezID), mean)
    
    # Create the exprs object here to feed to scripts further down. First make 
    # it a data.frame, name it "data" and assign column(ID_REF) and row(gene
    # entrez IDs) names.
    exprs<-as.data.frame(exp_matrix_unique[,2])
    colnames(exprs)[1] <- "ID_REF"
    exprs$ID_REF <- tolower(exprs$ID_REF)
    rownames(exprs) <- exp_matrix_unique[,1]

    #########							
    # Run pathway enrichment	
    #########

    gsdb <- get(genesets["human"])
    
    print("Running fingerprint")
    
    # Use SCE to calculate a score for each pathway based on the mean or median   
   
    SCE <- custom.single.chip.enrichment( 
        exprs = exprs,
        geneset = gsdb,
        transformation = transformation,
        statistic = statistic,
        progressBar = progressBar
    )
    
    #########							
    # Threshold according to GEO corpus background 	
    #########   			

	return(list(SCG = SCE, species = "human", platform = platform))
}

