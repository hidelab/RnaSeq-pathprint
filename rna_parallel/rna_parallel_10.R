# Fingerprinting worker script - running single chip enrichment for a group of
# experiments. Delay allows time to monitor process starting up

# Loading the neccessary scripts and files
source("/shared/hidelab2/shared/Sokratis/pathprint/fingerprint/package/GMAfunctions/R/GSMtable2exprs.R")
source("/shared/hidelab2/shared/Sokratis/pathprint_lincs/gctx2fingerprint.R")
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/genesets.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/pathprint.Hs.gs.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/sampleset.RData")

print("preparing for script launch, 5s synch delay")
Sys.sleep(5)

# Load the sample set in
current_samples <- sampleset[[10]];
total <- length(current_samples)

# load contents of directory once
# this used to be within the loop
e_path = "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/sq_gtex_experiments/"
gctx_e_path <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/gtex_experiments/"
dir <- dir(path = e_path)

for (i in 1:length(current_samples)){
	print(paste(i, "of", total, sep = " "))
    filename<-paste("sq", current_samples[i], sep = "_")
    
    # Only fingerprint experiments when you can't find their sq_file (experiment
    # fingerprint)
	if (filename %in% dir){
		print ("File already analyzed")
	}
	else if (!(filename %in% dir(path = e_path))){
	    #try-catching to skip erroneous experiments to avoid crushing
	    temp<-tryCatch(gctx2fingerprint(
    	    EXP = current_samples[i], 
    	    GEOthreshold = FALSE,
    	    EXPath = gctx_e_path,
    	    geneset = "KEGG and Wikipathways and static",
    	    enrichmentMethod = "SCE",
    	    transformation = "squared.rank",
    	    statistic = "mean",
    	    progressBar = FALSE
    	# Do nothing if an error occurs [add handlers later]*
        ),error=function(e){
		print(e)  
          })
	    
	    temp1<-list("name" = temp)
	    names(temp1)<-current_samples[i]
	    save(temp1, file = paste(e_path, filename, sep =""))
	}
}
