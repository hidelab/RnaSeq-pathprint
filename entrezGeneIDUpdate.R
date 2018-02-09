entrezGeneIDUpdate<-function(genelist, EntrezHistoryList = EntrezHistoryList) {
    temp <- as.character(genelist)
    discontinued<-intersect(temp, names(EntrezHistoryList))
    temp1 <- unlist(EntrezHistoryList[discontinued])
    for (i in 1:length(temp1)) {
        temp[temp %in% names(temp1)[i]] <- temp1[i]
    }
    return(temp)
}
