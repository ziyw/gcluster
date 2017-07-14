
# all cluster functions return data frame of genes and cluster

require(BHC)

cluster_bhc <- function(dataset, itemLabels, timePoints, dataType = "time-course"){
  # Avgs:
  #  dataset: dataset to do clustering
  #  itemLabels: name for all items
  #  timePoints: time points for datasets 
  #  dataType: default is 'time-course'
  
  # Return:
  #  cluster: a list of cluster the first row is genes, the second row is cluster
  
  bhc_cluster <- bhc(dataset,itemLabels,timePoints,dataType = dataType,verbose=TRUE)
  tmpFile <- 'tmp.txt'
  WriteOutClusterLabels(bhc_cluster, tmpFile, verbose=TRUE)
  
  cluster <- matrix( rep( 0, len=length(itemLabels)))
  rownames(cluster) <- as.vector(itemLabels)
  # read thefile back, get cluster list then delete the temp file
  
  
  count_cluster <- 0
  content <- file(tmpFile, "r")
  
  
  while (TRUE) {
    line <- readLines(content, n=1) 
    if (length(line) == 0) break
    
    if (!is.na(pmatch("---", line)) ) {
      count_cluster = count_cluster + 1
      next
    }
    
    index <- match(line, itemLabels)
    if (is.na(index)) next
    cluster[index] <- count_cluster
  }
  
  close(content)
  
  if (file.exists(tmpFile)) {
    file.remove(tmpFile)
  }

  return(cluster)
  
}

#filename <- "rand_data/polya_rand_1000.csv"

#dataset <- read.csv(filename)
#dataset <- dataset[1:10,]

#itemLabels <- t(dataset[1])
#dataset <- dataset[,2:6]
#timePoints <- c(1.5,3,6,9,12)

#cluster_bhc(dataset, itemLabels, timePoints)




