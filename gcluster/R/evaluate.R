#' Evaluation function
#' return the evalution of the cluster results
#' @param dataset The dataset to run cluster on
#' @param genes of the cluster 
#' @keywords evalution
#' @export
#' @examples
#' evalution(dataset, genes, cluster, measure = "all")
#' 

# The purpose of evaluation function is to run cluster method directly and get evaluation results

library("clValid")
library("org.Hs.eg.db")
library("GOSim")

# bhi_value
get_bhi_value <- function(genes, cluster) {
  # Get the GO term overlap results for cluster of genes
  # Avgs:
  #    genes: list of genes
  #    cluster: list of the assignment of cluster for each genes
  # Return:
  #    value : BHI score of the cluster
  new_genes <- vector()
  
  # set the first annotation as the result annotation
  for (i in 1:length(genes)){
    annotation <- unname(unlist(mget(x=genes[i],envir=org.Hs.egALIAS2EG)))
    new_genes[i] <- annotation[1]
  }
  
  if(require("Biobase") && require("annotate") && require("GO.db") &&
     require("org.Hs.eg.db")) {
    value <- BHI(cluster, annotation="org.Hs.eg.db", names=new_genes, category="all")
  }
  return(value)
}

# Cluster evaluation according to euclidean distance of the clusters
get_dunn_index <- function(dataset, cluster){	
  # need the dataset to evaluate the cluster
  Dist <- dist(dataset, method = "euclidean")	
  dunn_index <- dunn(Dist, cluster)
  return (dunn_index)
}

# Cluster evaluation according to gene annotation
get_goto_term <- function(genes, cluster) {
  # Get the GO term overlap results for cluster of genes
  # Avgs:
  #    genes: list of genes
  #    cluster: list of the assignment of cluster for each genes
  # Return:
  #    g_value : average silhouette value of the cluster
  
  # change gene names into annotation 
  new_genes <- vector()
  # set the first annotation as the result annotation
  for (i in 1:length(genes)){
    annotation <- unname(unlist(mget(x=genes[i],envir=org.Hs.egALIAS2EG)))
    new_genes[i] <- annotation[1]
  }
  
  # get the similarity matrix of genes
  sim <- getGeneSim(new_genes,similarity="dot",method="Tanimoto", verbose=FALSE, normalization = FALSE)
  # get rid of the NA/NaN/Inf elements in simliarity matrix 
  sim[is.na(sim)] <- 0
  sim[is.infinite(sim)] <- 0
  sim[is.nan(sim)] <- 0
  
  kept_genes <- rownames(sim)
  kept_cluster <- cluster[is.element(new_genes, kept_genes)]
  
  ev <- evaluateClustering(kept_cluster, sim)
  g_value <- mean(ev$clustersil)
  return(g_value)
}

evaluate <- function(dataset, genes, cluster, measure = "all") {

  
  denn_index <- get_dunn_index(dataset, cluster)
  bhi_value <- get_bhi_value(genes, cluster)
  goto_term <- get_goto_term(genes, cluster)
  
  eval <- data.frame('goto' = goto_term, "denn" = denn_index, "bhi" = bhi_value)
  return(eval)
}




