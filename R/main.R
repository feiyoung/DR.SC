## Add external data to DR.SC
# usethis::use_data(HCC1)
# Run to build the website
# pkgdown::build_site()
# pkgdown::build_reference()
# build_home()
# R CMD check --as-cran DR.SC_3.0.tar.gz

### Generate data without spatial dependence.
gendata_noSp <- function(n=100, p =100, q=15, K = 8,  alpha=10, sigma2=1, seed=1){
  
  q <- 2; K <- 10
  if(q <2) stop("error:gendata_noSp::  q must be greater than 2!")
  ## fixed after generation
  
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
  message("SNR=", round(snr,4), '\n')
  return(list(X=X, Z=Z, cluster=cluster, W=W, Mu=Mu, Sigma=Sigma, Pi=Pi, Lam_vec=Lambda, snr=snr))
}

### Generate data with spatial dependence.
gendata_sp <- function(height=30, width=30, p =100, q=10, K=7,  G=4, beta=1, sigma2=1, tau=8, seed=1, view=FALSE){
  # height <- 70
  # width <- 70
  # G <- 4
  # beta <- 1.0
  # K <- 7
  # q <- 10
  # p <- 1000
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  # require(GiRaF)
  # require(MASS)
  n <- height * width # # of cell in each indviduals 
  
  
  ## generate deterministic parameters, fixed after generation
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
  
  message("SNR=", round(snr,4), '\n')
  # make position
  pos <- cbind(rep(1:height, width), rep(1:height, each=width))
  
  return(list(X=X, Z=Z, cluster=y, W=W, Mu=Mu, Sigma=Sigma, Lam_vec=Lambda, beta=beta,pos=pos, snr=snr))
}


#### Generate Spatial data with ST platform
gendata_RNAExp <- function(height=30, width=30, platform="ST", p =100, q=10, K=7, 
                            G=4,sigma2=1, tau=8, seed=1, view=FALSE){
  
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  #require(GiRaF)
  #require(MASS)
  n <- height * width # # of cell in each indviduals 
  
  if(platform=="ST"){
    beta= 1
  }else if(platform=='scRNAseq'){
    beta = 0
  }
  ## generate deterministic parameters, fixed after generation
  
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
  X <- Z %*% t(W) + MASS::mvrnorm(n, mu=rep(0,p), Sigma=diag(Lambda))
  
  
  
  
  # make position
  pos <- cbind(rep(1:height, width), rep(1:height, each=width))
  
  counts <- t(X) - min(X)
  p <- ncol(X); n <- nrow(X)
  rownames(counts) <- paste0("gene", seq_len(p))
  colnames(counts) <- paste0("spot", seq_len(n))
  counts <- as.data.frame(exp(counts)-1)
  ## Make array coordinates - filled rectangle
  
  if(platform=="ST"){
    cdata <- list()
    cdata$row <- pos[,1]
    cdata$col <- pos[,2]
    cdata <- as.data.frame(do.call(cbind, cdata))
    cdata$imagerow <- cdata$row
    cdata$imagecol <- cdata$col 
    row.names(cdata) <- colnames(counts)
    
    #library(Seurat)
    ## Make SCE
    seu <-  CreateSeuratObject(counts= counts, meta.data=cdata) #
  }else if(platform=='scRNAseq'){
    # library(Seurat)
    ## Make SCE
    seu <-  CreateSeuratObject(counts= counts)
  }else{
    stop("gendata_RNAExp: Unsupported platform \"", platform, "\".")
  }
  
  seu$true_clusters <- y
  return(seu)
}


# Data reading ------------------------------------------------------------


read10XVisium <- function (dirname) {
  
  #require(Seurat)
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

readscRNAseq <- function(mtx, cells, features, ...){
  #require(Seurat)
  spMat <- ReadMtx(mtx,cells,features, ...)
  seu <- CreateSeuratObject(counts = spMat)
  seu@tools$platform <- "scRNAseq"
  seu
}



# Data prepocessing -------------------------------------------------------
## use SPARK to choose spatially variable genes
FindSVGs <- function(seu, nfeatures=2000, covariates=NULL, num_core=1, verbose=TRUE){
  
  if (!inherits(seu, "Seurat"))
    stop("method is only for Seurat objects")
  # require(SPARK)
  # require(Seurat)
  assy <- DefaultAssay(seu)
  sp_count <- seu[[assy]]@counts
  if(nrow(sp_count)>5000){
    seu <- FindVariableFeatures(seu, nfeatures = 5000, verbose=verbose)
    sp_count <- seu[[assy]]@counts[seu[[assy]]@var.features,]
  }
  location <- as.data.frame(cbind(seu$row, seu$col))
  if(verbose){
    message("Find the spatially variables genes by SPARK-X...\n")
  }
  sparkX <- sparkx(sp_count,location, X_in = covariates, numCores=num_core, option="mixture",  verbose=verbose)
  if(nfeatures > nrow(sp_count)) nfeatures <- nrow(sp_count)
  
  ## Find top nfeatures smallest adjusted p-values
  order_nfeatures <- order(sparkX$res_mtest$adjustedPval)[1:nfeatures]
  genes <- row.names(sp_count)[order_nfeatures]
  
  ## Access the gene based on gene name 
  is.SVGs <- rep(FALSE, nrow(seu))
  order.SVGs <- rep(NA, nrow(seu))
  adjusted.pval.SVGs <- rep(NA, nrow(seu))
  names(is.SVGs) <- names(order.SVGs)<- names(adjusted.pval.SVGs) <- row.names(seu)
  
  order.SVGs[genes] <- 1:length(genes)
  is.SVGs[genes] <- TRUE
  adjusted.pval.SVGs[genes] <- sparkX$res_mtest$adjustedPval[order_nfeatures]
  
  seu[[assy]]@meta.features$is.SVGs <- is.SVGs
  seu[[assy]]@meta.features$order.SVGs <- order.SVGs
  seu[[assy]]@meta.features$adjusted.pval.SVGs <- adjusted.pval.SVGs
  seu[[assy]]@var.features <- genes
  seu
}

topSVGs <- function(seu, ntop=5){
  if (!inherits(seu, "Seurat"))
    stop("method is only for Seurat objects")
  if(ntop > nrow(seu)) warning(paste0("Only ", nrow(seu), ' SVGs will be returned since the number of genes is less than ', ntop, '\n'))
  assy <- DefaultAssay(seu)
  if(is.null(seu[[assy]]@meta.features$is.SVGs)) 
    stop("There is no information about SVGs in default Assay. Please use function FindSVGs first!")
  
  SVGs <- row.names(seu)[seu[[assy]]@meta.features$is.SVGs]
  order_features <- seu[[assy]]@meta.features$order.SVGs
  
  idx <- order(order_features[!is.na(order_features)])[1:ntop]
  SVGs[idx]
  
}
# Define DR.SC S3 function ------------------------------------------------



DR.SC_fit <- function(X, K, Adj_sp=NULL, q=15,
                            error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                            maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                            wpca.int=FALSE,int.model="EEE", approxPCA=FALSE,  coreNum = 5){
  if (!(inherits(X, "dgCMatrix") || inherits(X, "matrix")))
    stop("X must be dgCMatrix object or matrix object")
  if(is.null(colnames(X))) colnames(X) <- paste0('gene', 1:ncol(X))
  ## Check whether X include zero-variance genes
  sd_zeros <- apply(X, 2, sd)
  if(sum(sd_zeros==0)>0){
    warning(paste0('There are ', sum(sd_zeros==0), ' zero-variance genes that will be removed!\n'))
    X <- X[, sd_zeros !=0]
  }
  if(length(K)==1) coreNum <- 1 

  message("Fit DR-SC model...\n")
  resList <- drsc(X,Adj_sp = Adj_sp, q=q, K=K, error.heter= error.heter, 
                  beta_grid=beta_grid,maxIter=maxIter, epsLogLik=epsLogLik,
                  verbose=verbose, maxIter_ICM=maxIter_ICM,
                  alpha=FALSE, wpca.int=wpca.int, diagSigmak=FALSE, int.model=int.model,
                  approxPCA=approxPCA, coreNum =coreNum)
  message("Finish DR-SC model fitting\n")
  
  return(resList)
}
DR.SC <- function(seu, K, q=15,  platform= c('Visium', "ST", "Other_SRT", "scRNAseq"), ...) {
  UseMethod("DR.SC")
}
  
DR.SC.Seurat <- function(seu, K, q=15,  platform= c('Visium', "ST", "Other_SRT", "scRNAseq"), ...){
  # require(Seurat)
  if (!inherits(seu, "Seurat"))
    stop("method is only for Seurat objects")
  
  platform <- match.arg(platform)
  if(platform == 'scRNAseq'){
    Adj_sp <- NULL
  }else{
    Adj_sp <- getAdj(seu,  platform)
  }
 
  
  assy <- DefaultAssay(seu)
  
  
  var.features <- seu[[assy]]@var.features
  
  X <- Matrix::t(seu[[assy]]@data[var.features, ])
  
  resList <- DR.SC_fit(X,Adj_sp = Adj_sp, q=q, K=K, ...)
  
  # hZ <-resList$hZ
  # row.names(hZ) <- colnames(seu)
  # colnames(hZ) <- paste0('DR-SC', 1:q)
  # seu@reductions$"dr-sc" <- CreateDimReducObject(embeddings = hZ, key='DRSC_', assay=DefaultAssay(seu))
  # seu$spatial.drsc.cluster <- resList$cluster
  # Idents(seu) <- factor(paste0("cluster", seu$spatial.drsc.cluster), levels=paste0('cluster',1:K))
  # seu@tools <- resList[-c(1,2)]
  seu[[assy]]@misc[['dr-scInf']] <- extractInfModel(resList)
  seu <- selectModel(seu,  criteria = 'MBIC', pen.const=0.5)
  
  return(seu)
  
}


extractInfModel <- function(resList){
  
  
  if(!inherits(resList, 'drscObject')) stop('extractInfModel: resList must be "drscObject" class!\n')
  
  nObj <- length(resList$Objdrsc)
  if(nObj<1) stop("extractInfModel: the length of 'resList' is zero!")
  if(nObj >= 1){
    loglik_vec <- rep(NA, nObj)
    dmat <- matrix(NA, nObj, 3)
    colnames(dmat) <- c("n", "p", "df")
    clusterList <- list()
    hZList <- list()
    for(i in 1:nObj){
      loglik_vec[i] <- resList$Objdrsc[[i]]$loglik
      dmat[i, ] <- resList$out_param[i, c(3,4,6)]
      clusterList[[i]] <- resList$Objdrsc[[i]]$cluster
      hZList[[i]] <- resList$Objdrsc[[i]]$hZ
    }
    
    KinfMat <- cbind(K=resList$K_set,loglik=loglik_vec,dmat)
    
    return(list(clusterList=clusterList, hZList=hZList, KinfMat=KinfMat))
  }
  
}

### This function includes the main methods: first is simultaneous dimension reduction and 
### clustering with homo variance error and no sptial information, second is the 
### drsc with heter variance and no spatial information, and third is the
### drsc with heter variance and spatial information
### This version use parallel package to evalute initial values, then use multi-thread in C++
### to evaluate Z and y in multi-K.
drsc <- function(X,Adj_sp = NULL, q, K, error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                           maxIter=30, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                           alpha=FALSE, wpca.int=TRUE, diagSigmak=FALSE, int.model="EEE",
                 approxPCA=FALSE, coreNum = 5){
  
  n <- nrow(X); p <- ncol(X)
  X <- scale(X, scale=FALSE)
  if(verbose){
    message("-------------------Calculate inital values-------------")
  }
  
  
  tic <- proc.time()
  if(approxPCA){
    message("Using approxmated PCA to obtain initial values")
    princ <- approxPCA(X, q)
  }else{
    message("Using accurate PCA to obtain initial values")
    princ <- wpca(X, q, weighted=wpca.int)
  }
  
  if(error.heter){
    Lam_vec0 <- princ$Lam_vec
  }else{
    Lam_vec0 <- rep(mean(princ$Lam_vec), p)
  }
  W0 <- princ$loadings
  hZ <- princ$PCs
  
  rm(princ) ## delete the temporary variable.
  
  nK <- length(K)
  if(nK>1){
    message("Starting parallel computing intial values...")
    ## set the number of cores in evaluating initial values
    num_core <- coreNum
    if(Sys.info()[1]=="Windows"){
      cl <- makeCluster(num_core)
    }else if(Sys.info()[1]=='Linux'){
      cl <- makeCluster(num_core, type="FORK")
    }else{
      cl <- makeCluster(num_core)
    }
    intList <- parLapply(cl, X=K, parfun_int, Z=hZ, alpha=alpha,int.model=int.model)
    stopCluster(cl)
  }else{
    intList <- list(parfun_int(K, hZ, alpha, int.model=int.model))
  }
  
  
  
  
  alpha0List = list()
  Mu0List = list()
  Sigma0List = list()
  ymat = matrix(0, n, nK)
  Pi0List = list()
  
  for (kk in 1: nK){
    
    ymat[,kk] <- intList[[kk]]$yveck
    alpha0List[[kk]] <-  intList[[kk]]$alpha0k
    Mu0List[[kk]] <- intList[[kk]]$Mu0k
    Sigma0List[[kk]] <- intList[[kk]]$Sigma0k
    
    Pi0List[[kk]] <- intList[[kk]]$Pi0k
  }
  rm(intList) ## delete the temporary variable.
  
  if(verbose){
    message("-------------------Finish computing inital values------------- ")
  }
  
  if(verbose){
    verbose <- 1
  }else{
    verbose <- 0
  }
  
  if(verbose)
    message("-------------------Starting  ICM-EM algortihm-------------")
  if((!is.null(Adj_sp))){
    
    beta0 = matrix(1, length(K), 1)*1.5
    resList <- icmem_heterCpp(X, Adj_sp, ymat,Mu0List, W0, Sigma0List,  Lam_vec0,
                              alpha0List,  beta_int=beta0, beta_grid=beta_grid, maxIter_ICM, maxIter, 
                              epsLogLik, verbose, !error.heter, diagSigmak, max(K), min(K), coreNum)
    
  }else if(is.null(Adj_sp)){
    
    resList <- EMmPCpp_heter(X, Pi0List, Mu0List, W0, Sigma0List, Lam_vec0, maxIter, epsLogLik, 
                             verbose, !error.heter, diagSigmak, max(K), min(K), coreNum)
  }
  
  toc_heter <- proc.time() - tic
  if(verbose) {
    message("-------------------Complete!-------------")
    message("elasped time is :", round(toc_heter[3], 2))
  }
  resList$K_set <- K
  
  class(resList) <- "drscObject"
  return(resList)
}



mycluster <- function(Z, G, int.model='EEE'){
  
  mclus2 <- Mclust(Z, G=G, modelNames=int.model,verbose=FALSE)
  return(mclus2)
}

parfun_int <- function(k, Z,  alpha, int.model='EEE'){
  
 
  mclus2 <- mycluster(Z, k, int.model)
   
  yveck <- mclus2$classification
  
  if(alpha){
    alpha0k <- mclus2$parameters$pro
  }else{
    alpha0k <- rep(0, k)
  }
  
  Mu0k <- t(mclus2$parameters$mean)
  Sigma0k <- mclus2$parameters$variance$sigma
  
  Pi0k <- mclus2$parameters$pro
  
  return(list(yveck=yveck, alpha0k=alpha0k, Mu0k = Mu0k, Sigma0k=Sigma0k, Pi0k=Pi0k))
}


selectModel <- function(obj, criteria = 'MBIC', pen.const=1){
  UseMethod("selectModel")
}

selectModel.drscObject <- function(obj, criteria = 'MBIC', pen.const=1){
  
  # select the best model based on the returned values from SimulDRcluster
  if(!inherits(obj, 'drscObject')) 
    stop('selectModel: method is only for Seurat or drscObject object!\n')
  
 
  reslist <- extractInfModel(obj)
  dfInf <- reslist$KinfMat
  K_set <- reslist$KinfMat[,"K"]
  nK <- length(K_set)
  icVec <-  rep(Inf, nK)
  for(k in 1:nK){
    # k <- 1
    
    n <- dfInf[k,3]; p <- dfInf[k,4];  dfree <- dfInf[k,5]
    loglik <- dfInf[k, 2]
    icVec[k] <- switch(criteria, 
                         MAIC = -2.0* loglik +dfree * 2 * log(log(p+n))*pen.const,
                         AIC = -2.0* loglik +dfree * 2,
                         MBIC =  -2.0* loglik +dfree * log(n) * log(log(p+n))*pen.const, 
                         BIC =  -2.0* loglik +dfree * log(n))
      
    
    
  }
  min_indx <- which.min(icVec)
  bestK <- K_set[min_indx]
  icMat <- cbind(K=K_set, IC=icVec)

  cluster_PCList <- list(bestK= bestK, cluster=as.vector(reslist$clusterList[[min_indx]]),
                         hZ = reslist$hZList[[min_indx]], icMat=icMat)
  return(cluster_PCList)
  
}

selectModel.Seurat <- function(obj,  criteria = 'MBIC', pen.const=1){
  if (!inherits(obj, "Seurat"))
    stop("selectModel: method is only for Seurat or drscObject object")
  
  assy <- DefaultAssay(obj)
  
  reslist <- obj[[assy]]@misc[['dr-scInf']] 
  dfInf <- reslist$KinfMat
  K_set <- reslist$KinfMat[,"K"]
  nK <- length(K_set)
  
  icMat <- matrix(Inf, nK, 3)
  colnames(icMat) <- toupper(c("MBIC", "BIC", "AIC"))
  for(k in 1:nK){
    # k <- 1
    
    n <- dfInf[k,3]; p <- dfInf[k,4];  dfree <- dfInf[k,5]
    loglik <- dfInf[k, 2]
    icMat[k, ] <-c(MBIC =  -2.0* loglik +dfree * log(n) * log(log(p+n))*pen.const, 
                     BIC =  -2.0* loglik +dfree * log(n),
                       AIC = -2.0* loglik +dfree * 2)
    
    
    
  }
  criteria <- toupper(criteria)
  icVec <- icMat[,criteria]
  min_indx <- which.min(icVec)
  bestK <- K_set[min_indx]
  icMat <- cbind(K=K_set, icMat)
  
  hZ <- reslist$hZList[[min_indx]]
  row.names(hZ) <- colnames(obj)
  colnames(hZ) <- paste0('DR-SC', 1: ncol(hZ))
  obj@reductions$"dr-sc" <- CreateDimReducObject(embeddings = hZ, key='DRSC_', assay=assy)
  obj$spatial.drsc.cluster <- as.vector(reslist$clusterList[[min_indx]])
  Idents(obj) <- factor(paste0("cluster", obj$spatial.drsc.cluster), levels=paste0('cluster',1: bestK))
  obj[[assy]]@misc[['icMat']] <- icMat 
  
  return(obj)
}

getAdj <- function(obj, platform = c('Visium', "ST", "Other_SRT"), ...) UseMethod("getAdj")

getAdj.Seurat <- function(obj, platform = c('Visium', "ST", "Other_SRT"), ...){
  
  if (!inherits(obj, "Seurat"))
    stop("method is only for Seurat or matrix objects")
  # require(Matrix)
  
  platform <- match.arg(platform)
  if (tolower(platform) == "visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (tolower(platform) == "st") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  }else if(tolower(platform) == 'other_srt'){
    pos <- as.matrix(cbind(row=obj$row, col=obj$col))
    Adj_sp <- getAdj_auto(pos, ...)
    
    return(Adj_sp)
  }else{
    stop("getAdj: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  spot.positions <- as.data.frame(cbind(row=obj$row, col=obj$col))
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
          ncol(obj), " spots.")
  
  n <- length(df_j) 
  
  D <- matrix(0,  nrow = n, ncol = n)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  ij <- which(D != 0, arr.ind = T)
  
  Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1, dims=c(n, n))
  return(Adj_sp)
}

## Bisection method to search the optimal radius to make the  median of neighborhoods between 4~6.
find_neighbors <- function(pos, platform=c('ST', "Visium")) {
  # require(purrr)
  # require()
  if (tolower(platform) == "visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (tolower(platform) == "st") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop("find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  colnames(pos) <- c("row", "col")
  pos <- DataFrame(pos)
  spot.positions <- pos
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
          nrow(pos), " spots.")
  
  n <- length(df_j) 
  
  D <- matrix(0,  nrow = n, ncol = n)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  ij <- which(D != 0, arr.ind = T)
  ij
}
getAdj_reg <- function(pos, platform ='Visium'){
    # require(Matrix)
    ij <- find_neighbors(pos, platform)
    n <- nrow(pos)
    Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1, dims=c(n, n))
    return(Adj_sp)
}

getAdj_auto <- function(pos, lower.med=4, upper.med=6, radius.upper= NULL){
  if (!inherits(pos, "matrix"))
    stop("method is only for  matrix object!")
  
  
  ## Automatically determine the upper radius
  n_spots <- nrow(pos)
  idx <- sample(n_spots, min(100, n_spots))
  dis <- dist(pos[idx,])
  if(is.null(radius.upper)){
    #radius.upper <- max(dis)
    radius.upper <- sort(dis)[20] ## select the nearest 20 spots.
  }
  radius.lower <- min(dis[dis>0])
  Adj_sp <- getneighborhood_fast(pos, radius=radius.upper)
  Med <- summary(Matrix::rowSums(Adj_sp))['Median']
  if(Med < lower.med) stop("The radius.upper is too smaller that cannot find median neighbors greater than 4.")
  start.radius <- 1
  Med <- 0
  message("Find the adjacency matrix by bisection method...")
  maxIter <- 30
  k <- 1
  while(!(Med >= lower.med && Med <=upper.med)){ # ensure that each spot has about 4~6 neighborhoods in median.
    
    Adj_sp <- getneighborhood_fast(pos, radius=start.radius)
    Med <- summary(Matrix::rowSums(Adj_sp))['Median']
    if(Med < lower.med){
      radius.lower <- start.radius
      start.radius <- (radius.lower + radius.upper)/2
    }else if(Med >upper.med){
      radius.upper <- start.radius
      start.radius <- (radius.lower + radius.upper)/2
    }
    message("Current radius is ", round(start.radius, 2))
    message("Median of neighborhoods is ", Med)
    if(k > maxIter) {
      message("Reach the maximum iteration but can not find a proper radius!")
      break;
    }
    k <- k + 1
  }
  
  return(Adj_sp)
}



getAdj_manual <- function(pos, radius){
  if (!inherits(pos, "matrix"))
    stop("pos must be a matrix!")
  if(radius <=0){
    stop('radius must be a positive real!')
  }
  
  Adj_sp <- getneighborhood_fast(pos, radius=radius)
  return(Adj_sp)
}





# Approximated PCA for fast computation--------------------------------------------------------------

approxPCA <- function(X, q){ ## speed the computation for initial values.
  # require(irlba) 
  n <- nrow(X)
  svdX  <- irlba(A =X, nv = q)
  PCs <- svdX$u %*% diag(svdX$d[1:q])
  loadings <- svdX$v
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}


# Weighted PCs ------------------------------------------------------------
#  @Wei Liu
#  This function includes  functions to impliment Weighted PCA in 
#  reference Bai, J., & Liao, Y. (2013). Statistical inferences using large estimated covariances for panel data and factor models. arXiv preprint arXiv:1307.2662.
#  and Inferences in panel data with interactive effects using largecovariance matrices
#  It considers the heterogeneous error term in approximated factor model.
wpca <- function(X, q, weighted=TRUE){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  out <- wpcaCpp(X, q, weighted)
  return(out)
}
wpca2 <- function(X, q, weighted=TRUE){
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

RunWPCA <- function(object,q=15) UseMethod("RunWPCA")

RunWPCA.Seurat <- function(object, q=15){
  
  if (!inherits(object, "Seurat"))
    stop("method is only for Seurat, dgCMatrix or matrix objects")
  if(length(object@assays[[DefaultAssay(object)]]@scale.data) == 0)
    stop("The slot scale.data is empty! Please use ScaleData function first then use RunWPCA!")
  scale.data <- object@assays[[DefaultAssay(object)]]@scale.data
  hZ <- wpca(t(scale.data), q=q, weighted = T)$PCs
  row.names(hZ) <- colnames(object)
  colnames(hZ) <- paste0('WPCA', 1:q)
  object@reductions$"wpca" <- CreateDimReducObject(embeddings = hZ, key='WPCA_', assay=DefaultAssay(object))
  return(object)
}

RunWPCA.matrix <- function(object, q=15){
  if (!inherits(object, "matrix"))
    stop("method is only for Seurat, dgCMatrix or matrix objects")
  hZ <- wpca(t(object), q=q, weighted = T)
  if(is.null(colnames(object))){
    warning('colnames(object) is null, so the colnames are assigned with spot 1:ncol(object)!')
    colnames(object) <- paste0('spot', 1:ncol(object))
  }
  row.names(hZ) <- colnames(object)
  colnames(hZ) <- paste0('WPCA', 1:q)
  return(hZ)
}
RunWPCA.dgCMatrix<- function(object, q=15){
  if (!inherits(object, "dgCMatrix"))
    stop("method is only for Seurat, dgCMatrix or matrix objects")
  hZ <- wpca(Matrix::t(object), q=q, weighted = T)
  if(is.null(colnames(object))){
    warning('colnames(object) is null, so the colnames are assigned with spot 1:ncol(object)!')
    colnames(object) <- paste0('spot', 1:ncol(object))
  }
  row.names(hZ) <- colnames(object)
  colnames(hZ) <- paste0('WPCA', 1:q)
  return(hZ)
}







# Previous core functions  Code in previous package-------------------------------------------------


simulDRcluster <- function(X,Adj_sp = NULL, q, K, error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                           maxIter=30, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,pen.const=0.5,
                           alpha=FALSE, wpca.int=TRUE, diagSigmak=FALSE){
  
  n <- nrow(X); p <- ncol(X)
  X <- scale(X, scale=FALSE)
  if(verbose){
    message("-------------------Calculate inital values-------------")
  }
  
  # require(mclust)
  tic <- proc.time()
  princ <- wpca(X, q, weighted=wpca.int)
  if(error.heter){
    Lam_vec0 <- princ$Lam_vec
  }else{
    Lam_vec0 <- rep(mean(princ$Lam_vec), p)
  }
  
  W0 <- princ$loadings
  hZ <- princ$PCs
  mclus2 <- Mclust(hZ, G=K, verbose=FALSE)
  toc_gmm <- proc.time() - tic
  
  y <- mclus2$classification
  if(alpha){
    alpha0 <- mclus2$parameters$pro
  }else{
    alpha0 <- rep(0, K)
  }
  
  Mu0 <- t(mclus2$parameters$mean)
  Sigma0 <- mclus2$parameters$variance$sigma
  
  alpha0List = list()
  Mu0List = list()
  Sigma0List = list()
  
  ymat = matrix(0, n, 1)
  Pi0List = list()
  
  for (kk in 1: 1){
    
    ymat[,kk] <- y
    alpha0List[[kk]] <-  alpha0
    Mu0List[[kk]] <- Mu0
    Sigma0List[[kk]] <- Sigma0
    
    Pi0List[[kk]] <- mclus2$parameters$pro
  }
  
  if(verbose){
    message("-------------------Finish computing inital values------------- ")
  }
  
  
  
  
  if(verbose){
    verbose <- 1
  }else{
    verbose <- 0
  }
  
  if(verbose)
    message("-------------------Starting  ICM-EM algortihm-------------")
  if((!is.null(Adj_sp))){
    
    resList <- icmem_heterCpp(X, Adj_sp, ymat,Mu0List, W0, Sigma0List,  Lam_vec0,
                              alpha0List,  beta_int=1.5, beta_grid=beta_grid, maxIter_ICM, maxIter, 
                              epsLogLik, verbose, !error.heter, diagSigmak, max(K), min(K), 1)
    resList <- resList[[1]][[1]]
    resList$aic <- -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* 2* log(log(p+n))*pen.const # adjusted  bic and aic for high dimension
    resList$bic <-  -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* log(n)* log(log(p+n))*pen.const 
    
  }else if(is.null(Adj_sp)){
    
    resList <- EMmPCpp_heter(X, Pi0List, Mu0List, W0, Sigma0List, Lam_vec0, maxIter, 
                             epsLogLik, 
                             verbose, !error.heter, diagSigmak, max(K), min(K), 1)
    resList <- resList[[1]][[1]]
    resList$aic <- -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* 2* log(log(p+n))*pen.const # adjusted  bic and aic for high dimension
    resList$bic <-  -2.0* resList$loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* log(n)* log(log(p+n)) *pen.const
  }
  
  toc_heter <- proc.time() - tic
  if(verbose) {
    message("-------------------Complete!-------------")
    message("elasped time is :", round(toc_heter[3], 2))
  }
  resList$cluster_init <- y
  time_used <- c(toc_gmm[3], toc_heter[3])
  names(time_used) <- c("pcgmm", "simul")
  resList$time <- time_used
  return(resList)
}
selectClustNumber <- function(X,Adj_sp=NULL, q, K_set= 3:10, parallel="parallel", num_core = 1,...){
  
  
  nK <- length(K_set)
  if(!is.null(parallel)){
    if (num_core > 1) {
      if (num_core > parallel::detectCores()) {
        warning("selectClustNumber:: the number of cores you're setting is larger than detected cores!")
        num_core = parallel::detectCores()
      }
    }
    
    if(parallel=='parallel'){
      #library(parallel)
      
      cl <- parallel::makeCluster(num_core)
      # parallel::clusterExport(cl, list("simulDRcluster"))
      ## "EMmPCpp_heter", "icmem_heterCpp",
      message("Starting parallel computing...")
      # clusterCall(cl, function() library(MixPPCA))
      # Run
      icMat <- parSapply(cl, X=K_set, parafun1, XX=X, Adj_sp=Adj_sp, q=q, ...)
      parallel::stopCluster(cl)
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
