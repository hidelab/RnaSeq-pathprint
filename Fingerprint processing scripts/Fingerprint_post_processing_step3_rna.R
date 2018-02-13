# Post processing fingerprint step 3
# Author: Sokratis Kariotis
# Description: Creates the *sq_.matrix.xxxx-xx-xx.RData* POE matrix file (one file for all pathways). 
##########
# PART 3 #
##########

# Define paths and headers
header<-"sq_"
setwd("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/")

Pathwaypath <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/pathwayFrames/"
POEpath <- "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/pathwayFrames/POE/"

# Define the location of the pathway frame
dir()[grep("platform", dir())]
filename <- "platform.frame.2018-02-12.RData" 
print(paste("Loading platform data: ", filename, sep = ""))
load(filename)

table(platform.frame$Platform)
platforms<-as.character(unique(platform.frame$Platform))

print("Got platforms")

# Load list of pathways
library(pathprint)
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/genesets.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/pathprint.Hs.gs.rda")

pathwayNames<-names(get(genesets[[1]]))

# Setup POE.matrix
POE.matrix<-matrix(nrow = length(pathwayNames), ncol = length(platform.frame$Platform))
rownames(POE.matrix)<-pathwayNames
colnames(POE.matrix)<-platform.frame$GEO

print(platforms)
print(length(platforms))

for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("Inserting POE values for ", array, sep = ""))
  for (i in 1:length(pathwayNames)){
    print(paste(POEpath, header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    load(paste(POEpath, header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    POE.matrix[i,names(pathway.POE)]<-pathway.POE
    }
  }

print(colnames(POE.matrix))

save(POE.matrix, file = paste(header, ".POE.matrix.PRE.", Sys.Date(), ".RData", sep = ""))

# Remove arrays that are NA for the first pathway - this is usually due to a mis-match between the platform species and the annotated species
POE.matrix<-POE.matrix[,!(is.na(POE.matrix[1,]))]

save(POE.matrix, file = paste(header, ".POE.matrix.", Sys.Date(), ".RData", sep = ""))


POE.matrix.0.001<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
# How many NAs? - These are probably due to the species or chips where a pathway is not represented properly
(sum(is.na(POE.matrix.0.001)))/(dim(POE.matrix.0.001)[1]*dim(POE.matrix.0.001)[2])
# How many zeros?
(sum((POE.matrix.0.001 == 0), na.rm = TRUE))/(dim(POE.matrix.0.001)[1]*dim(POE.matrix.0.001)[2])

save(POE.matrix.0.001, file = paste("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/", header, ".POE.matrix.0.001.",
                              Sys.Date(), ".RData",
                              sep = ""
                              )
    )

# END OF PART 3
#################
