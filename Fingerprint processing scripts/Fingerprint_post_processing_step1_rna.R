# Post processing fingerprint step 1
# Author: Sokratis Kariotis
# Description: Creates the *sq_.frame.xxxx-xx-xx.RData* file (SCE dataframe) and 
# the *platform.frame.xxxx-xx-xx.RData* file (platform dataframe). It also creates 
# the *sq_.pathway.SCE.x.RData* files for each pathway.
##########
# PART 1 #
##########

# set directory holding fingerprint collection
fingerpath<-"/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/sq_gtex_experiments/"

# define file identifier
header<-"sq_"
setwd(fingerpath)

load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/genesets.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/pathprint.Hs.gs.rda")

print("Loading files")

files<-dir(fingerpath)

probe.data <- as.data.frame(matrix(nrow = length(files), ncol = 3))
colnames(probe.data) <- c("Experiment","GPL","nrow")
probe.data$Experiment <- gsub(header, "", files)
probe.data$GPL <- " GTEx_Analysis_V7"
probe.data$nrow <- length(files)

# for each platform build a table of the data table nrows and find the max (i.e. complete array)
platforms<-levels(as.factor(probe.data$GPL))
platforms.nrow<-lapply(platforms, function(x){table(probe.data$nrow[probe.data$GPL == x])})
names(platforms.nrow)<-platforms

save(probe.data, file = "/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/probe.Rdata")
print("Saved probe.RData")

setwd(fingerpath)

print("Normalizing using the POE approach")
print("constructing dataframes")

# construct list of fingerprints - SCGs (single chip gene set enrichments) - and platforms
# This displays a running total of the arrays processed for tracking
SCE<-vector("list", length(files))
platform<-vector("list", length(files))

# Go through the experiment files and update the probeset removing those that didn't get downloaded
for (i in 1:length(files)){
    if(file.exists(files[i])) {
        load(files[i])
        print(i)
	    print(files[i])
        SCE[[i]]<-temp1[[1]]$SCG
        platform[[i]]<-temp1[[1]]$platform
    }
    else
        print(paste("Experiment file not in sq_gtex_experiments folder: ", files[i]))
 }

# load geneset from data repository to give pathway names
library(pathprint)

# create SCE dataframe
print(paste("saving SCE dataframe with timestamp ", Sys.Date(), sep = " "))
gsub(header, "", files)
names(SCE)<-gsub(header, "", files)
SCE.frame<-t(as.data.frame(SCE))
colnames(SCE.frame)<-names(get(genesets[[1]]))
rownames(SCE.frame)<-gsub(header, "", files)
save(SCE.frame, file = paste("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/", header, ".frame.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

# create platform dataframe with timestamp
print(paste("saving platform dataframe with timestamp ", Sys.Date(), sep = " "))
names(platform)<-gsub(header, "", files)
platform.frame<-data.frame(GEO = names(platform), Platform = unlist(platform))
save(platform.frame, file = paste("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/platform.frame",
                                  Sys.Date(), "RData",
                                  sep = "."
                                  )
                              )

####
print("combining equivalent platforms")

table(platform.frame$Platform)
platforms<-as.character(unique(platform.frame$Platform))
#####
print("Saving pathway-specific data")
library(metaArray)

for (i in 1:ncol(SCE.frame)){
assign("pathway.SCE", SCE.frame[,i])
save(pathway.SCE, file = paste("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/pathwayFrames/", header, 
                  ".pathway.SCE.", i, ".RData", sep = ""))
}

# END OF PART 1
###################

