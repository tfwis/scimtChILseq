GetGraphDistance_ <- function(cell.embedding,dims,K) {
  message("Constructing kNN graph")
  knn.result <- FNN::get.knn(cell.embedding[,dims], k = K)
  KNN.adj.mat <- matrix(0, nrow = nrow(cell.embedding), ncol = nrow(cell.embedding))
  for (i in 1:nrow(KNN.adj.mat)) {
    KNN.adj.mat[i, knn.result$nn.index[i, ]] <- 1
  }
  message("Constructing graph distance matrix")
  knn.graph <- igraph::graph_from_adjacency_matrix(KNN.adj.mat, 
                                                   mode = "undirected")
  graph.dist.mat <- igraph::distances(knn.graph)
  if (is.infinite(max(graph.dist.mat))) 
    stop("The cell-cell kNN graph has disconnected components. Please increase K.")
  message(sprintf("The largest graph distance is %s", max(graph.dist.mat)))
  graph.dist.mat
}

GetGeneEmbedding_ <- function (dist.mat, K = 10, sigma = NULL, nEV = 30, t = 1) {
  nEV <- min(nEV,ncol(dist.mat)-1)
  dm.res <- GeneTrajectory::diffusion.map(dist.mat, K = K, sigma = sigma, nEV = nEV + 1, t = t)
  dm.res[["diffu.emb"]] <- dm.res[["diffu.emb"]][, 2:(nEV + 1)]
  colnames(dm.res[["diffu.emb"]]) <- paste0("DM_", 1:ncol(dm.res[["diffu.emb"]]))
  dm.res[["eigen.vals"]] <- dm.res[["eigen.vals"]][2:(nEV + 1)]
  dm.res
}

pre.coarse.grain <- function(cell.embedding, graph.dist, N = 1000, random.seed = 1) {
  message("Run k-means clustering")
  set.seed(random.seed)
  km.res <- stats::kmeans(cell.embedding, N,iter.max = 100)
  message("Coarse-grain matrices")
  KNN.membership.mat <- matrix(0, nrow = N, ncol = nrow(cell.embedding))
  for (i in 1:ncol(KNN.membership.mat)) {
    KNN.membership.mat[km.res$cluster[i], i] <- 1
  }
  KNN.membership.mat <- KNN.membership.mat/apply(KNN.membership.mat, 1, sum)
  graph.dist.updated <- KNN.membership.mat %*% graph.dist %*% t(KNN.membership.mat)
  colnames(KNN.membership.mat) <- rownames(cell.embedding)
  output <- list()
  output[["KNN.membership.mat"]] <- KNN.membership.mat
  output[["graph.dist"]] <- graph.dist.updated
  return(output)
}

exp.coarse.grain <- function (pre.cg.output, gene.expression) {
  KNN.membership.mat <- biclust::binarize(pre.cg.output[['KNN.membership.mat']],0)
  int <- base::intersect(colnames(KNN.membership.mat),colnames(gene.expression))
  #gene.expression <- scale(t(gene.expression[,int]),center = F)
  gene.expression.updated <- KNN.membership.mat[,int] %*% t(gene.expression[,int])
  
  output <- list()
  output[["gene.expression"]] <- gene.expression.updated
  output[["graph.dist"]] <- pre.cg.output[["graph.dist"]]
  output[["features"]] <- colnames(gene.expression.updated)
  output
}
