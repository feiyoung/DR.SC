## Add external data to DR.SC
# usethis::use_data(HCC1)
# Run to build the website
# pkgdown::build_site()
### Generate data without spatial dependence.
gendata_noSp <- function(n=100, p =100, q=15, K = 8,  alpha=10, sigma2=1, seed=1){
  
  q <- 2; K <- 10
  if(q <2) stop("error:gendata_noSp::  q must be greater than 2!")
  ## fixed after generation
  set.seed(1)
  Lambda = sigma2*(abs(rnorm(p, sd=1)))
  W1 <- matrix(rnorm(p*q), p, q)
  W <- qr.Q(qr(W1))
  
  # set different random seed for each repeat.
  set.seed(seed)
  Pi <- rep(1, K)/K
  Y <- t(rmultinom(n, size=1, Pi))
  cluster <- apply(Y, 1, which.max)
  
  mu <- matrix(0, q,  K)
  diagmat = array(0, dim = c(q, q, K))
  if(q > K){
    q1 <- floor(K/2)
    for(j in 1:q1){
      if(j <= (q1/2)) mu[j,j] <- alpha
      if(j > (q1/2)) mu[j,j] <- -alpha
    }
    mu[(q1+1):q, K] <- -alpha
    
  }else if(q <= K){
    for(k in 1:K)
       mu[,k] <- rep(alpha/8 *k, q) #
  }
  for(k in 1:K){
    tmp  <- rep(1, q)
    if(k <= K/2){
      tmp[q] <- alpha
    }
    diag(diagmat[,,k]) <- tmp
  }
  
  
  Mu <- t(mu)
  Ez <- 0
  for(k in 1:K){
    
    Ez <- Ez + Pi[k]* Mu[k,]
  }
  Mu <- Mu - matrix(Ez, K, q, byrow=T)
  
  Sigma <- diagmat
  Z <- matrix(0, n, q)
  for(k in 1:K){
    nk <- sum(Y[,k])
    Z[Y[,k]==1, ] <- MASS::mvrnorm(nk, Mu[k,], Sigma[,,k])
  }
  
  X = Z %*% t(W) + MASS::mvrnorm(n, rep(0,p), diag(Lambda))
  
  svd_Sig <- svd(cov(Z))
  W12 <- W %*% svd_Sig$u %*% diag(sqrt(svd_Sig$d))
  signal <- sum(svd(W12)$d^2)
  snr <- sum(svd(W12)$d^2) / (sum(svd(W12)$d^2)+ sum(Lambda))
  cat("SNR=", round(snr,4), '\n')
  return(list(X=X, Z=Z, cluster=cluster, W=W, Mu=Mu, Sigma=Sigma, Pi=Pi, Lam_vec=Lambda, snr=snr))
}



# dat1 <- gendata_noSp()
# dat1$Mu
### Generate data with spatial dependence.
gendata_sp <- function(height=30, width=30, p =100, q=10, K=7,  G=4, beta=1, sigma2=1, tau=8, seed=1, view=T){
  # height <- 70
  # width <- 70
  # G <- 4
  # beta <- 1.0
  # K <- 7
  # q <- 10
  # p <- 1000
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  require(GiRaF)
  require(MASS)
  n <- height * width # # of cell in each indviduals 
  
  
  ## generate deterministic parameters, fixed after generation
  set.seed(1)
  # sigma2 <- 1
  Lambda <- sigma2*abs(rnorm(p, sd=1))
  W1 <- matrix(rnorm(p*q), p, q)
  W <- qr.Q(qr(W1))
  mu <- matrix(0, q,  K)
  diagmat = array(0, dim = c(q, q, K))
  if(q > K){
    q1 <- floor(K/2)
    for(j in 1:q1){
      if(j <= (q1/2)) mu[j,j] <- tau
      if(j > (q1/2)) mu[j,j] <- -tau
    }
    mu[(q1+1):q, K] <- -tau
    
  }else if(q <= K){
    for(k in 1:K)
      mu[,k] <- rep(tau/8 *k, q) #
  }
  for(k in 1:K){
    tmp  <- rep(1, q)
    if(k <= K/2){
      tmp[q] <- tau
    }
    diag(diagmat[,,k]) <- tmp
  }
  
  
  Mu <- t(mu)
  Sigma <- diagmat
  set.seed(seed)
  # generate the spatial dependence for state variable y, a hidden Markov RF
  y <- sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width, ncolors = K, nei = G, param = beta,
                   initialise = FALSE, view = view)
  y <- c(y) + 1
  
  Z <- matrix(0, n, q)
  for(k in 1:K){
    nk <- sum(y==k)
    Z[y==k, ] <- MASS::mvrnorm(nk, Mu[k,], Sigma[,,k])
  }
  Ez <- colMeans(Z)
  Mu <- Mu - matrix(Ez, K, q, byrow=T) # center Z
  X <- Z %*% t(W) + MASS::mvrnorm(n, rep(0,p), diag(Lambda))
  
  svd_Sig <- svd(cov(Z))
  W12 <- W %*% svd_Sig$u %*% diag(sqrt(svd_Sig$d))
  signal <- sum(svd(W12)$d^2)
  snr <- sum(svd(W12)$d^2) / (sum(svd(W12)$d^2)+ sum(Lambda))
  
  cat("SNR=", round(snr,4), '\n')
  
  # make position
  pos <- cbind(rep(1:height, width), rep(1:height, each=width))
  
  return(list(X=X, Z=Z, cluster=y, W=W, Mu=Mu, Sigma=Sigma, Lam_vec=Lambda, beta=beta,pos=pos, snr=snr))
}


#### Generate Spatial data with ST platform
gendata_spatial <- function(height=30, width=30, p =100, q=10, K=7,  G=4, beta=1, sigma2=1, tau=8, seed=1, view=F){
  
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  require(GiRaF)
  require(MASS)
  n <- height * width # # of cell in each indviduals 
  
  
  ## generate deterministic parameters, fixed after generation
  set.seed(1)
  # sigma2 <- 1
  Lambda <- sigma2*abs(rnorm(p, sd=1))
  W1 <- matrix(rnorm(p*q), p, q)
  W <- qr.Q(qr(W1))
  mu <- matrix(0, q,  K)
  diagmat = array(0, dim = c(q, q, K))
  if(q > K){
    q1 <- floor(K/2)
    for(j in 1:q1){
      if(j <= (q1/2)) mu[j,j] <- tau
      if(j > (q1/2)) mu[j,j] <- -tau
    }
    mu[(q1+1):q, K] <- -tau
    
  }else if(q <= K){
    for(k in 1:K)
      mu[,k] <- rep(tau/8 *k, q) #
  }
  for(k in 1:K){
    tmp  <- rep(1, q)
    if(k <= K/2){
      tmp[q] <- tau
    }
    diag(diagmat[,,k]) <- tmp
  }
  
  
  Mu <- t(mu)
  Sigma <- diagmat
  set.seed(seed)
  # generate the spatial dependence for state variable y, a hidden Markov RF
  y <- sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width, ncolors = K, nei = G, param = beta,
                   initialise = FALSE, view = view)
  y <- c(y) + 1
  
  Z <- matrix(0, n, q)
  for(k in 1:K){
    nk <- sum(y==k)
    Z[y==k, ] <- MASS::mvrnorm(nk, Mu[k,], Sigma[,,k])
  }
  Ez <- colMeans(Z)
  Mu <- Mu - matrix(Ez, K, q, byrow=T) # center Z
  X <- Z %*% t(W) + MASS::mvrnorm(n, rep(0,p), diag(Lambda))
  
  svd_Sig <- svd(cov(Z))
  W12 <- W %*% svd_Sig$u %*% diag(sqrt(svd_Sig$d))
  signal <- sum(svd(W12)$d^2)
  snr <- sum(svd(W12)$d^2) / (sum(svd(W12)$d^2)+ sum(Lambda))
  
  cat("SNR=", round(snr,4), '\n')
  
  # make position
  pos <- cbind(rep(1:height, width), rep(1:height, each=width))
  #  make BayesSpace metadata used in BayesSpace-------------------------------------------------
  counts <- t(X) - min(X)
  p <- ncol(X); n <- nrow(X)
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- paste0("spot_", seq_len(n))
  
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  row.names(cdata) <- colnames(counts)
  library(Seurat)
  ## Make SCE
  seu <-  CreateSeuratObject(counts= exp(counts)-1, meta.data=cdata)
  seu$true_clusters <- y
  return(seu)
}


# Data reading ------------------------------------------------------------


read10XVisium <- function (dirname) {
  
  require(Seurat)
  uniquifyFeatureNames <- function (ID, names)
  {
    if (length(ID) != length(names)) {
      stop("lengths of 'ID' and 'names' must be equal")
    }
    names <- as.character(names)
    ID <- as.character(ID)
    missing.name <- is.na(names)
    names[missing.name] <- ID[missing.name]
    dup.name <- names %in% names[duplicated(names)]
    names[dup.name] <- paste0(names[dup.name], "_", ID[dup.name])
    return(names)
  }
  spatial_dir <- file.path(dirname, "spatial")
  matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_dir))
    stop("Matrix directory does not exist:\n  ", matrix_dir)
  if (!dir.exists(spatial_dir))
    stop("Spatial directory does not exist:\n  ", spatial_dir)
  colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"),
                      header = FALSE)
  colnames(colData) <- c("spot", "in_tissue", "row",
                         "col", "imagerow", "imagecol")
  rownames(colData) <- colData$spot
  colData <- colData[colData$in_tissue > 0, ]
  rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"),
                        header = FALSE)
  
  colnames(rowData) <- c("gene_id", "gene_name",
                         "feature_type")
  rowData <- rowData[, c("gene_id", "gene_name")]
  rownames(rowData) <- uniquifyFeatureNames(rowData$gene_id,
                                            rowData$gene_name)
  counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"),
                         header = FALSE)
  colnames(counts) <- barcodes$V1
  rownames(counts) <- rownames(rowData)
  counts <- counts[, rownames(colData)]
  
  seu <- CreateSeuratObject(counts = counts, meta.data = colData)
  seu[['RNA']] <- AddMetaData(seu[['RNA']], metadata=rowData)
  seu@tools$platform <- "Visium"
  seu
}


# Data prepocessing -------------------------------------------------------
## use SPARK to choose spatially variable genes
FindSVGs <- function(seu, nfeatures=2000,num_core=1, verbose=T){
  require(SPARK)
  require(Seurat)
  sp_count <- seu@assays$RNA@counts
  if(nrow(sp_count)>5000){
    seu <- FindVariableFeatures(hcc1, nfeatures = 5000, verbose=verbose)
    sp_count <- seu@assays$RNA@counts[seu@assays$RNA@var.features,]
  }
  location <- as.data.frame(cbind(seu$row, seu$col))
  if(verbose){
    cat("Find the spatially variables genes by SPARK-X...")
  }
  sparkX <- sparkx(sp_count,location,numCores=num_core,option="mixture", verbose=verbose)
  if(nfeatures > nrow(sp_count)) nfeatures <- nrow(sp_count)
  genes <- row.names(sp_count)[order(sparkX$res_mtest$adjustedPval)[1:nfeatures]]
  is.SVGs <- rep(F, nrow(seu))
  names(is.SVGs) <- row.names(seu)
  is.SVGs[genes] <- T
  seu[['RNA']]@meta.features$is.SVGs <- is.SVGs
  seu
}

# Define DR.SC S3 function ------------------------------------------------


DR.SC <- function(object, ...) UseMethod("DR.SC")
DR.SC.fit <- function(X,Adj_sp=NULL, q=15, K= NULL,
                            error.heter= T, K_set = seq(2, 10), beta_grid=seq(0.5, 5, by=0.5),
                            maxIter=25, epsLogLik=1e-5, verbose=F, maxIter_ICM=6,pen.const=1,
                            wpca.int=F, parallel='doSNOW', num_core=5){
  # if (!inherits(X, "dgCMatrix"))
  #   stop("X must be dgCMatrix objects or Seurat objects")
  if(is.null(K)){
    cat("Start chooing number of clusters...\n")
    cat(" The candidate set is: ", K_set, '\n')
    icMat <- selectClustNumber(X, Adj_sp, q, K_set= K_set, parallel=parallel, num_core = num_core, pen.const=pen.const)
    K_best <- K_set[which.min(icMat[,"mbic"])]
    cat("The best number of cluster is: ", K_best, '\n')
  }else{
    K_best <- K
  }
  resList <- simulDRcluster(X,Adj_sp = Adj_sp, q, K_best, error.heter= error.heter, 
                            beta_grid=beta_grid,
                            maxIter=maxIter, epsLogLik=epsLogLik, verbose=verbose, maxIter_ICM=maxIter_ICM,pen.const=pen.const,
                            alpha=F, wpca.int=wpca.int, diagSigmak=FALSE)
  if(is.null(K)){
    resList$K_best <- K_best
    resList$icMat <- icMat
  }
  return(resList)
}

DR.SC.Seurat <- function(seu, q=15, K=NULL, platform= c("Visium", "ST", 'scRNAseq'),
                         nfeatures=2000,K_set = seq(2, 10), variable.type=c("HVGs", "SVGs"),...){
  require(Seurat)
  if (!inherits(seu, "Seurat"))
    stop("method is only for Seurat objects")
  if(platform == 'scRNAseq'){
    Adj_sp <- NULL
  }else{
    Adj_sp <- getAdj(seu,  platform)
  }
  if(variable.type=='HVGs'){
    if(is.null(seu@assays$RNA@var.features)){
      seu <- FindVariableFeatures(seu, nfeatures = nfeatures)
    }
    var.features <- seu@assays$RNA@var.features
  
  }else{
    var.features <- row.names(seu)[seu[["RNA"]]@meta.features$is.SVGs]
  }
  
 
  
  X <- Matrix::t(LogNormalize(seu@assays$RNA@counts[var.features,]))
  
  resList <- DR.SC.fit(X,Adj_sp = Adj_sp, q=q, K=K,K_set =K_set, ...)
  if(is.null(K)){
    K <- resList$K_best
  }
  hZ <-resList$hZ
  row.names(hZ) <- colnames(seu)
  colnames(hZ) <- paste0('DR-SC', 1:q)
  seu@reductions$"dr-sc" <- CreateDimReducObject(embeddings = hZ, key='DRSC_', assay='RNA')
  seu$spatial.drsc.cluster <- resList$cluster
  Idents(seu) <- factor(paste0("cluster", seu$spatial.drsc.cluster), levels=paste0('cluster',1:K))
  seu@tools <- resList[-c(1,2)]
  return(seu)
  
}




### This function includes the main methods: first is simultaneous dimension reduction and 
### clustering with homo variance error and no sptial information, second is the 
### simulDRcluster with heter variance and no spatial information, and third is the
### simulDRcluster with heter variance and spatial information
simulDRcluster <- function(X,Adj_sp = NULL, q, K, error.heter= T, beta_grid=seq(0.5, 5, by=0.5),
                           maxIter=30, epsLogLik=1e-5, verbose=F, maxIter_ICM=6,pen.const=0.5,
                           alpha=F, wpca.int=T, diagSigmak=FALSE){
  
  n <- nrow(X); p <- ncol(X)
  X <- scale(X, scale=F)
  if(verbose){
    cat("-------------------Calculate inital values------------- \n")
  }
  
  require(mclust)
  tic <- proc.time()
  princ <- wpca(X, q, weighted=wpca.int)
  if(error.heter){
    Lam_vec0 <- princ$Lam_vec
  }else{
    Lam_vec0 <- rep(mean(princ$Lam_vec), p)
  }
  
  W0 <- princ$loadings
  hZ <- princ$PCs
  set.seed(1)
  mclus2 <- Mclust(hZ, G=K)
  toc_gmm <- proc.time() - tic
  
  y <- mclus2$classification
  if(alpha){
    alpha0 <- mclus2$parameters$pro
  }else{
    alpha0 <- rep(0, K)
  }
  
  Mu0 <- t(mclus2$parameters$mean)
  Sigma0 <- mclus2$parameters$variance$sigma
  
  if(verbose){
    cat("-------------------Finish computing inital values------------- \n")
  }
  
  
  
  
  if(verbose){
    verbose <- 1
  }else{
    verbose <- 0
  }
  
  if(verbose)
    cat("-------------------Starting  EM algortihm------------- \n")
  if((!is.null(Adj_sp))){
    
    resList <- icmem_heterCpp(X, Adj_sp, y,  Mu0, W0, Sigma0,  Lam_vec0,
                              alpha=alpha0,  beta_int=1.5, beta_grid=beta_grid, maxIter_ICM, maxIter, 
                              epsLogLik, verbose, !error.heter, diagSigmak)
    resList$aic <- -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* 2* log(log(p+n))*pen.const # adjusted  bic and aic for high dimension
    resList$bic <-  -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* log(n)* log(log(p+n))*pen.const 
    
  }else if(is.null(Adj_sp)){
    alpha0 <- mclus2$parameters$pro
    resList <- EMmPCpp_heter(X, alpha0, Mu0, W0,Sigma0, Lam_vec0,maxIter, epsLogLik, 
                             verbose, !error.heter, diagSigmak)
    resList$aic <- -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* 2* log(log(p+n))*pen.const # adjusted  bic and aic for high dimension
    resList$bic <-  -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* log(n)* log(log(p+n)) *pen.const
  }
  
  toc_heter <- proc.time() - tic
  if(verbose) cat("-------------------Complete!------------- \n")
  cat("elasped time is :", toc_heter, '\n')
  resList$cluster_init <- y
  time_used <- c(toc_gmm[3], toc_heter[3])
  names(time_used) <- c("pcgmm", "simul")
  resList$time <- time_used
  return(resList)
}


find_neighbors <- function(sce, platform='ST') {
  
  if (platform == "Visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (platform == "ST") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  spot.positions <- colData(sce)[, c("col", "row")]
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  
  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  
  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions), 
                     as.data.frame(spot.positions), 
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)
  
  ## Shift to zero-indexing for C++
  #neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
  
  ## Group neighbor indices by spot 
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  
  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  ## df_j <- map(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  df_j <- lapply(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  
  ## Log number of spots with neighbors
  n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
          ncol(sce), " spots.")
  
  n <- length(df_j) 
  
  D <- matrix(0,  nrow = n, ncol = n)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  ij <- which(D != 0, arr.ind = T)
  ij
}

getAdj <- function(seu, platform = c('Visium', 'ST')){
  
  if (!inherits(seu, "Seurat"))
    stop("method is only for Seurat objects")
  require(Matrix)
  if (platform == "Visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (platform == "ST") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  spot.positions <- as.data.frame(cbind(row=seu$row, col=seu$col))
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  
  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  
  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions), 
                     as.data.frame(spot.positions), 
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)
  
  ## Shift to zero-indexing for C++
  #neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
  
  ## Group neighbor indices by spot 
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  
  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  ## df_j <- map(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  df_j <- lapply(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  
  ## Log number of spots with neighbors
  n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
          ncol(seu), " spots.")
  
  n <- length(df_j) 
  
  D <- matrix(0,  nrow = n, ncol = n)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  ij <- which(D != 0, arr.ind = T)
  
  Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1)
  return(Adj_sp)
}

runAdj <- function(X, pos, platform='ST'){
  require(purrr)
  require(SingleCellExperiment)
  
  # make sce structure
  n <- nrow(X)
  p <- ncol(X)
  counts <- t(X)
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- paste0("spot_", seq_len(n))
  
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  ## Make SCE
  ## note: scater::runPCA throws warning on our small sim data, so use prcomp
  sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
  ij <- find_neighbors(sce, platform)
  library(Matrix)
  Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1)
  return(Adj_sp)
}


# select cluster number K -------------------------------------------------



selectClustNumber <- function(X,Adj_sp, q, K_set= 3:10, parallel="parallel", num_core = 1,...){
  
  
  nK <- length(K_set)
  if(!is.null(parallel)){
    if (num_core > 1) {
      if (num_core > detectCores()) {
        warning("selectClustNumber:: the number of cores you're setting is larger than detected cores!")
        num_core = detectCores()
      }
    }
    if(parallel=='doSNOW'){
      require(doSNOW)
      cl <- makeSOCKcluster(num_core)
      registerDoSNOW(cl)
      
      ## set Prgogress Bar
      pb <- txtProgressBar(min=1, max=nK, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      k <- 1
      icMat <- foreach(k = 1:nK,.packages="MixPPCA" ,.options.snow=opts,
                       .combine='rbind') %dopar% {
                         reslist <- simulDRcluster(X,Adj_sp = Adj_sp, q=q, K=K_set[k],  ...) 
                         
                         c(reslist$bic, reslist$aic)
                       }
      close(pb)
      stopCluster(cl)
    }else if(parallel=='parallel'){
      library(parallel)
      
      cl <- makeCluster(num_core)
      clusterExport(cl, list("EMmPCpp_heter", "icmem_heterCpp", "simulDRcluster"))
      cat("Starting parallel computing...")
      # clusterCall(cl, function() library(MixPPCA))
      # Run
      icMat <- parSapply(cl, X=K_set, parafun1, XX=X, Adj_sp=Adj_sp, q=q, ...)
      stopCluster(cl)
      icMat <- t(icMat)
    }
      
  }else{
    icMat <- matrix(NA, nK, 2)
    pb <- txtProgressBar()
    for(k in 1:nK){
      reslist <- simulDRcluster(X,Adj_sp = Adj_sp, q=q, K=K_set[k],  ...) 
      setTxtProgressBar(pb, k)
      icMat[k, ] <- c(reslist$bic, reslist$aic)
    }
    close(pb)
  }
  
  
  
  
  icMat <- cbind(K_set, icMat)
  colnames(icMat) <- c("K", 'mbic', 'aic')
  row.names(icMat) <- as.character(K_set)
  return(icMat)
}
parafun1 <- function(K, XX, Adj_sp, q, ...){
  reslist <- simulDRcluster(XX,Adj_sp = Adj_sp, q=q, K=K,  ...) 
  
  return(c(reslist$bic, reslist$aic))
}


# select factor number q --------------------------------------------------

selectFacNumber <- function(X, qmax=15){
  mnlamjFun <- function(eigvals, j){
    p <- length(eigvals)
    lamj <- eigvals[j]
    Sum <- 0
    for(l in (j+1):p){
      Sum <- Sum + 1/(eigvals[l] - lamj)
    }
    res <- Sum + 1/ ((3*lamj + eigvals[j+1])/4 - lamj)
    return(res/(p-j))
  }
  mtnlamjFun <- function(n, eigvals, j){
    p <- length(eigvals)
    rhojn <-  (p-j)/(n-1)
    res <- -(1-rhojn)/ eigvals[j] + rhojn * mnlamjFun(eigvals, j)
    return(res)
  }
  ##Reference: Fan, J., Guo, J., & Zheng, S. (2020). Estimating number of factors by adjusted eigenvalues thresholding. Journal of the American Statistical Association, 1-10.
  n <- nrow(X)
  p <- ncol(X)
  corMat <- cor(X)
  evalues <- eigen(corMat)$values
  hq1 <- sum(evalues>1+sqrt(p/(n-1)))
  if(hq1 < 15){
    hq <- hq1
  }else{ # ajdust the eigvalues
    adj.eigvals <- sapply(1:(p-1), function(j) -1/mtnlamjFun(n, evalues, j))
    hq <- sum(adj.eigvals >1) # overselect
  }
  if(hq > qmax || hq < 5) hq <- qmax
  
  
  
  
  
  propvar <- sum(evalues[1:hq]) / sum(evalues)
  res <- list()
  res$q <- hq
  res$propvar <- sum(evalues[1:hq]) / sum(evalues)
  
  return(res)
}



# Weighted PCs ------------------------------------------------------------
#  @Wei Liu
#  This function includes  functions to impliment Weighted PCA in 
#  reference Bai, J., & Liao, Y. (2013). Statistical inferences using large estimated covariances for panel data and factor models. arXiv preprint arXiv:1307.2662.
#  and Inferences in panel data with interactive effects using largecovariance matrices
#  It considers the heterogeneous error term in approximated factor model.
wpca <- function(X, q, weighted=T){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  out <- wpcaCpp(X, q, weighted)
  return(out)
}
wpca2 <- function(X, q, weighted=T){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  n <- nrow(X)

  svdX <- svd(X, nu=q, nv=q)
  PCs <- svdX$u %*% diag(svdX$d[1:q])
  loadings <- svdX$v
  dX <- PCs %*% t(loadings) - X
  Lam_vec = colSums(dX^2)/ n
  if(weighted){
    svdXw <- svd(X %*% diag(1/sqrt(Lam_vec)), nu=q, nv=q)
    PCs <- svdXw$u %*% diag(svdXw$d[1:q])
    loadings <- diag(sqrt(Lam_vec)) %*% svdXw$v
    dX <- PCs %*% t(loadings) - X
    Lam_vec = colSums(dX^2)/ n
  }
  return(list(PCs=PCs, loadings=loadings, Lam_vec=Lam_vec))
}
