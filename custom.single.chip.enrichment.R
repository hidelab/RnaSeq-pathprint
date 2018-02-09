custom.single.chip.enrichment <-
function( 
   exprs,
   geneset,
   transformation = "rank",
   statistic = "mean",
   normalizedScore = FALSE,
   progressBar = TRUE
   )

###################
# custom.single.chip.enrichment - function to assess erichment of gene sets using various statistics
# exprs - an expression matrix, rownames correspond to gene ids in the geneset
# geneset - list of pathways or geneset over which to assess statistic
# transformation - "rank", "squared.rank" or "log.rank" - applied to each column of exprs
# statistic - "mean" or "median" - summary statistic to be applied
# progressBar - shows progress of script, good to check running okay, set to FALSE for possible faster running
##################

{
  # check parameters
  if(!(transformation %in% c("rank", "squared.rank", "log.rank"))) stop("transformation should be rank, squared.rank or log.rank")
  if(!(statistic %in% c("mean", "median"))) stop("transformation should be mean or median")
  
  # Sample normalization
  Ns <- ncol(exprs)
  Ng <- nrow(exprs)
  gene.names<-rownames(exprs)
  geneset.names<-names(geneset)
  
  # rank
  exprs<-apply(exprs, 2, rank, ties.method = "average")
  if (transformation == "log.rank") {
      exprs <- log(exprs) 
      }
  else if (transformation == "squared.rank") {
      exprs <- exprs^2
      }

  # Loop statistic over gene sets
  if(progressBar == TRUE){
    pb<-txtProgressBar(min = 0, max = length(geneset), style = 3)
    }
  score.matrix <- matrix(0, nrow = length(geneset), ncol = Ns)
  for (i in 1:length(geneset)) {
      overlap <- intersect(geneset[[i]], gene.names)      
      # print(paste(i, "gene set:", geneset.names[i], " overlap=", length(overlap))) # for testing
      if (length(overlap) == 0) {
         score.matrix[i, ] <- NA
      } else {
          if (statistic == "mean"){
          #score.matrix[i, ] <- apply(exprs, 2, function(x){mean(x[match(overlap, gene.names)])})
          score.matrix[i, ] <- apply(exprs, 2, function(x){mean(x[overlap])}) # this runs slightly faster          
          if(normalizedScore == TRUE){
            n <- length(overlap)
            # print(i) # for testing
            # print(n) # for testing
            # print(Ng) # for testing
            # expected sample mean = population mean
            if (transformation == "rank"){
              E.mean <- mean(1:Ng)
              E.sd <- ((sd(1:Ng)/(n^0.5)))*(((Ng-n)/(Ng-1))^0.5)
              }
            else if (transformation == "log.rank"){
              E.mean <- mean(log(1:Ng))
              E.sd <- ((sd(log(1:Ng))/(n^0.5)))*(((Ng-n)/(Ng-1))^0.5)
              }
            else if (transformation == "squared.rank"){
              E.mean <- mean((1:Ng)^2)
              E.sd <- ((sd((1:Ng)^2)/(n^0.5)))*(((Ng-n)/(Ng-1))^0.5)
              }
            
            # population sd needs to be corrected for selection without replacement
            # e.g. http://courses.wcupa.edu/rbove/Berenson/10th%20ed%20CD-ROM%20topics/section7_3.pdf
            # print(E.mean) # for testing        
            # print(E.sd) # for testing
            # Use these parameters to normalize to a score between -0.5 and +0.5
            score.matrix[i, ] <- sapply(score.matrix[i, ], pnorm, mean = E.mean, sd = E.sd) - 0.5
              }          
            }
          else if (statistic == "median"){
            #score.matrix[i, ] <- apply(exprs, 2, function(x){median(x[match(overlap, gene.names)])})       
            score.matrix[i, ] <- apply(exprs, 2, function(x){median(x[overlap])}) # this runs slightly faster
            }
        }         
        if(progressBar == TRUE){
          setTxtProgressBar(pb, i)
          }
      }
    colnames(score.matrix)<-colnames(exprs)
    rownames(score.matrix)<-geneset.names
    
    return(score.matrix)
    
    } # end of single.chip.enrichment

