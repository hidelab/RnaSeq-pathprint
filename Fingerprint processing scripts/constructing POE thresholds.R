# Post processing construction of the POE  thresholds for the fingerprint step 5
# Author: Sokratis Kariotis
# Description: Script to obtain threshold values for the fingerprint. These values will 
# be passed to a dataframe, thresholds.RData. This will be included into the 
# fingerprinting R package. The background parameters are determined by the POE method.
# The threshold was optimized against a panel of tissue samples

# Define path of repository
definePath<-function(){
    readline("define pathprint repository, or blank for default (/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/) : ")
}
pathprintRepository<-definePath()

if (pathprintRepository == ""){
    pathprintRepository<-"/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/"
}

# load raw SCE and POE matrix, need to select appropriate dates
dataPath<-"/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/"

# To make it more versatile
# print("Relevant available files:")
# dir(path = dataPath, pattern = "POE.matrix.20")
# print("select POE file")

POEfile<-"sq_.POE.matrix.2018-02-13.RData"

# dir(path = dataPath, pattern = ".20")
# print("select SCE file")

SCEfile<-"sq_.frame.2018-02-12.RData"

# dir(path = dataPath,
# pattern = "platform.frame.20")
# print("select platform file")

platformfile<-"platform.frame.2018-02-12.RData"

load(paste(dataPath, POEfile, sep = ""))
load(paste(dataPath, SCEfile, sep = ""))
load(paste(dataPath, platformfile, sep = ""))

# N.B. the POE.matrix is slightly smaller as a few arrays were removed
dim(SCE.frame)
dim(POE.matrix)
dim(platform.frame)

# re-organize
SCE.frame<-t(SCE.frame)
colnames(SCE.frame)<-gsub("SCE_", "", colnames(SCE.frame))

platform.frame$GEO<-gsub("SCE_", "", platform.frame$GEO)
rownames(platform.frame)<-platform.frame$GEO

platform.frame<-platform.frame[rownames(platform.frame) %in% colnames(POE.matrix),]

SCE.frame<-SCE.frame[,colnames(SCE.frame) %in% colnames(POE.matrix)]

# check dimensions and orders
dim(SCE.frame)
dim(POE.matrix)
dim(platform.frame)
sum(colnames(SCE.frame) == colnames(POE.matrix))
sum(platform.frame$GEO == colnames(POE.matrix))

# Now need to define thresholds for each platform
# N.B. this produces warnings and NaN values where pathway scores are NA
# That's okay for now
platforms<-unique(as.character(platform.frame$Platform))
print(platforms)
# set threshold as 0.001
threshold<-0.001
platform.thresholds<-vector("list", length(platforms))
names(platform.thresholds)<-platforms
print("Check that the outputs are all close to defined threshold level")
for (i in 1:length(platform.thresholds)){
  # define submatrix for each platform
  POE.subMatrix<-POE.matrix[,platform.frame$Platform == names(platform.thresholds)[i]]
  SCE.subMatrix<-SCE.frame[,platform.frame$Platform == names(platform.thresholds)[i]]
  # setup matrix to contain threshold values for each pathway
  thresholdTable<-matrix(nrow = nrow(POE.subMatrix), ncol = 2)
  rownames(thresholdTable)<-rownames(POE.subMatrix)
  colnames(thresholdTable)<-c("high", "low")

# find closest reference to threshold
  for (j in 1:nrow(POE.subMatrix)){
    # obtain locations of points closest to chosen threshold
    highThresholdRef<-which(abs(threshold - POE.subMatrix[j,]) == min(abs(threshold - POE.subMatrix[j,]), na.rm = TRUE))
    lowThresholdRef<-which(abs((-threshold) - POE.subMatrix[j,]) == min(abs((-threshold) - POE.subMatrix[j,]), na.rm = TRUE))
    # use the mean value of these points as the threshold
    # print for error checking
    print(c(mean(POE.subMatrix[j,highThresholdRef]), mean(POE.subMatrix[j,highThresholdRef])))
    # insert into matrix
    thresholdTable[j,1]<-mean(SCE.subMatrix[j,highThresholdRef])
    thresholdTable[j,2]<-mean(SCE.subMatrix[j,lowThresholdRef])
  }
  platform.thresholds[[i]]<-thresholdTable
  }

print(i)  
save(platform.thresholds, file = paste(dataPath, "POE_derived_thresholds_0.001.",
    Sys.Date(), ".RData",
    sep = ""
    )
)

# also save into fingerprinting repository
try(system(paste("cp ",
             pathprintRepository,
             "platform.thresholds.RData ",
             dataPath,
             "platform.thresholds.RData.old",
             sep = "")))
try(system(paste("cp ",
             pathprintRepository,
             "RNA.fingerprint.matrix.RData ",
             dataPath,
             "RNA.fingerprint.matrix.RData.old",
             sep = "")))
save(platform.thresholds, file = paste(pathprintRepository, "platform.thresholds.RData", sep = ""))

RNA.fingerprint.matrix<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
save(RNA.fingerprint.matrix, file = paste(pathprintRepository, "RNA.fingerprint.matrix.RData", sep = ""))

print(names(platform.thresholds))
# END
