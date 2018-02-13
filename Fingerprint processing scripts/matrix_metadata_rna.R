# Post processing metadata matrix creation step 4
# Author: Sokratis Kariotis
# Description: Script to compile a table of metadata for the fingerprint matrix.
# Due to the nature of the metadata found in GSE70138_Broad_LINCS_inst_info.txt, 
# the LINCS.metadata.matrix will differ from the correspodning GEO file. In addition
# to the present columns, we have GPL and GSE columns. Subject to change!

# Define path of repository
definePath<-function(){
readline("define pathprint repository, or blank for default (/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/) : ")
}
pathprintRepository<-definePath()

if (pathprintRepository == ""){
	pathprintRepository<-"/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/"
	}

# Load fingerprint matrix
defineFile<-function(){
readline("enter filename: ")
}
dataPath<-"/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/post_processing_general_files/"

dir(path = dataPath, pattern = "platform.frame.20")
print("select platform file")
platformfile<-"platform.frame.2018-02-12.RData"
load(paste(dataPath, platformfile, sep = ""))

# Read metadata table
RNA.metadata.matrix = read.delim("/shared/hidelab2/shared/Sokratis/pathprint_rnaseq/rnaseq_Folder/GTEx_v7_Annotations_SampleAttributesDS.txt", fill = TRUE, header = TRUE, sep = "\t")

# Construct metadata table
RNA.metadata.matrix <- as.data.frame(append(RNA.metadata.matrix, list(Species = "Homo sapiens"), after = 1))
RNA.metadata.matrix <- as.data.frame(append(RNA.metadata.matrix, list(Analysis = "GTEx_Analysis_V7"), after = 1))

# Filtering out incomplete lines
#LINCS.metadata.matrix <- LINCS.metadata.matrix[!(is.na(LINCS.metadata.matrix$inst_id == "") | LINCS.metadata.matrix$cell_id == "" | LINCS.metadata.matrix$det_plate == ""), ]

# Convert to ascii strings for better package compatibility
RNA.metadata.matrix<-RNA.metadata.matrix

# first convert "ï" to "i"
# RNA.metadata.matrix.temp[,5] <- gsub("ï", "i", RNA.metadata.matrix.temp[,5])
# RNA.metadata.matrix.temp[,6] <- gsub("ï", "i", RNA.metadata.matrix.temp[,6])
# RNA.metadata.matrix.temp[,7] <- gsub("ï", "i", RNA.metadata.matrix.temp[,7])
# 
# LINCS.metadata.matrix.temp[,5] <- iconv(LINCS.metadata.matrix.temp[,5], "latin1", "ASCII", "byte")
# LINCS.metadata.matrix.temp[,6] <- iconv(LINCS.metadata.matrix.temp[,6], "latin1", "ASCII", "byte")
# LINCS.metadata.matrix.temp[,7] <- iconv(LINCS.metadata.matrix.temp[,7], "latin1", "ASCII", "byte")

# Show changes
#all.equal(LINCS.metadata.matrix.temp, LINCS.metadata.matrix)

# LINCS.metadata.matrix<-LINCS.metadata.matrix.temp

# Save file into fingerprint package
try(system(paste("cp ",
             pathprintRepository,
             "RNA.metadata.matrix.RData ",
             dataPath,
             "RNA.metadata.matrix.RData.old",
             sep = "")))

save(RNA.metadata.matrix, file = paste(pathprintRepository, "RNA.metadata.matrix.RData", sep = ""), compress = TRUE)

# end
