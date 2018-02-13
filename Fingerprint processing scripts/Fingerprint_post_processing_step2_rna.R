# Post processing fingerprint step 2
# Author: Sokratis Kariotis
# Description: Creates the POE (probability of expression) file for each pathway. 
##########
# PART 2 #
##########

# define paths and headers
header<-"sq_"
setwd("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/")

# Define the location of the pathway frame
dir()[grep("platform", dir())]
filename <- "platform.frame.2018-02-12.RData" 
print(paste("Loading platform data: ", filename, sep = ""))
load(filename)

# create directory to store POE files
Pathwaypath <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/pathwayFrames/"
POEpath <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/pathwayFrames/POE/"
try(system(paste("mkdir ", POEpath, header, sep = "")))

# this script is run across 14 nodes

library(metaArray)
library(foreach)
# library(multicore)
library(doMC)
registerDoMC(cores = 14)

print(table(platform.frame$Platform))
platforms<-as.character(unique(platform.frame$Platform))

# load list of pathways
library(pathprint)
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/genesets.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/pathprint.Hs.gs.rda")

pathwayNames<-names(get(genesets[[1]]))

print("processing POE for each platform")
print (length(platforms))
print(platforms)

for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("calculating POE for ", array, sep = ""))
  # for every pathway perform POE
  foreach (i = 1:length(pathwayNames)) %dopar% {
    tname <- paste(POEpath, header, "/pathway.SCE.", i, ".", array, ".RData", sep = "")
    print(tname)
    if (!file.exists(tname)) {
    # load the SCG values for each pathway
    load(paste(Pathwaypath, header, ".pathway.SCE.", i, ".RData", sep = ""))
    print("loaded")
    # subset for the array to be assessed
    sample <- pathway.SCE[grep(array, platform.frame$Platform)]
    print("sample")
    # setup POE vector
    pathway.POE<-sample
    pathway.POE[]<-NA
    # subset to avoid NAs
    sample.valid<-sample[!(is.na(sample))]
    fit<-fit.em(sample.valid,
            cl = rep(0,length(sample.valid)),TRUE
            )
    # write to POE file
    pathway.POE[!(is.na(sample))]<-fit$expr
    save(pathway.POE, file = paste(POEpath, header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    }
    else {
    print("found one")
    }
    }
  }

# END OF PART 2
######################
