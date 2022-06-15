# Call sparkx -------------------------------------------------------------


sparkx <- function (count_in, locus_in, X_in = NULL, numCores = 1, option = "mixture", 
                    verbose = TRUE){
  if (is(count_in, "matrix")) {
    raw_count_mat <- as(as.matrix(count_in), "sparseMatrix")
  }
  else if (is(count_in, "vector")) {
    raw_count_mat <- as(t(as.matrix(count_in)), "sparseMatrix")
  }
  else if (is(count_in, "sparseMatrix")) {
    raw_count_mat <- count_in
  }
  else {
    stop("counts has to be of following forms: vector,matrix or sparseMatrix")
  }
  rm(count_in)
  locus_in <- as.matrix(locus_in)
  totalcount <- as.vector(sp_sums_Rcpp(raw_count_mat))
  keep_cell_idx <- which(totalcount != 0)
  count_mat <- raw_count_mat[, keep_cell_idx]
  fil_loc <- locus_in[keep_cell_idx, ]
  if (!is.null(X_in)) {
    X_in <- X_in[keep_cell_idx, ]
  }
  rm(raw_count_mat)
  keep_gene_idx <- which(as.vector(sp_sums_Rcpp(count_mat, 
                                                rowSums = TRUE)) != 0)
  count_mat <- count_mat[keep_gene_idx, ]
  numGene <- nrow(count_mat)
  numCell <- ncol(count_mat)
  if(verbose){
    message(paste("## ===== SPARK-X INPUT INFORMATION ===="))
    message(paste("## number of total samples:", numCell))
    message(paste("## number of total genes:", numGene))
    if (numCores > 1) {
      message(paste("## Running with", numCores, "cores \n"))
    }
    else {
      message(paste("## Running with single core, may take some time \n"))
    }
  }
  
  GeneNames <- rownames(count_mat)
  if (sum(is.na(GeneNames)) > 0) {
    GeneNames[is.na(GeneNames)] <- "NAgene"
  }
  sparkx_list <- list()
  icount = 0
  if (option == "mixture") {
    if(verbose)
       message(paste0("## Testing With Projection Kernel"))
    icount = icount + 1
    final_location <- fil_loc
    sparkx_list[[icount]] <- sparkx.sk(count_mat, final_location, 
                                       X_mat = X_in, mc_cores = numCores, verbose = verbose)
    rm(final_location)
    for (iker in 1:5) {
      icount = icount + 1
      if(verbose)
         message(paste0("## Testing With Gaussian Kernel ", iker))
      final_location <- apply(fil_loc, 2, transloc_func, 
                              lker = iker, transfunc = "gaussian")
      sparkx_list[[icount]] <- sparkx.sk(count_mat, final_location, 
                                         X_mat = X_in, mc_cores = numCores, verbose = verbose)
      rm(final_location)
    }
    for (iker in 1:5) {
      icount = icount + 1
      if(verbose)
      message(paste0("## Testing With Cosine Kernel ", iker))
      final_location <- apply(fil_loc, 2, transloc_func, 
                              lker = iker, transfunc = "cosine")
      sparkx_list[[icount]] <- sparkx.sk(count_mat, final_location, 
                                         X_mat = X_in, mc_cores = numCores, verbose = verbose)
      rm(final_location)
    }
    names(sparkx_list) <- c("projection", paste0("gaus", 
                                                 1:5), paste0("cos", 1:5))
  }
  else if (option == "single") {
    if(verbose)
       message(paste0("## Testing With Projection Kernel"))
    icount = icount + 1
    final_location <- fil_loc
    sparkx_list[[icount]] <- sparkx.sk(count_mat, final_location, 
                                       X_mat = X_in, mc_cores = numCores, verbose = verbose)
    rm(final_location)
    names(sparkx_list) <- c("projection")
  }
  else {
    stop("option should be one of following: single or mixture")
  }
  allstat <- sapply(sparkx_list, function(x) {
    x$stat
  })
  allpvals <- sapply(sparkx_list, function(x) {
    x$pval
  })
  rownames(allstat) <- rownames(allpvals) <- GeneNames
  comb_pval <- apply(allpvals, 1, ACAT)
  pBY <- p.adjust(comb_pval, method = "BY")
  joint_pval <- cbind.data.frame(combinedPval = comb_pval, 
                                 adjustedPval = pBY)
  final_res <- list(stats = allstat, res_stest = allpvals, 
                    res_mtest = joint_pval)
  return(final_res)
}


sparkx.sk <- function(counts, infomat, X_mat = NULL, mc_cores = 1, verbose = TRUE) 
{
  geneName <- rownames(counts)
  if (sum(is.na(geneName)) > 0) {
    geneName[is.na(geneName)] <- "NAgene"
  }
  if (is.null(X_mat)) {
    Xinfomat <- apply(infomat, 2, scale, scale = FALSE)
    loc_inv <- solve(crossprod(Xinfomat, Xinfomat))
    kmat_first <- Xinfomat %*% loc_inv
    LocDim <- ncol(infomat)
    Klam <- eigen(crossprod(Xinfomat, kmat_first), only.values = T)$values
    EHL <- counts %*% Xinfomat
    numCell <- nrow(Xinfomat)
    adjust_nominator <- as.vector(sp_sums_Rcpp(counts^2, 
                                               TRUE))
    vec_stat <- apply(EHL, 1, function(x) {
      x %*% loc_inv %*% as.matrix(x)
    }) * numCell/adjust_nominator
    vec_ybar <- as.vector(sp_means_Rcpp(counts, TRUE))
    vec_ylam <- unlist(parallel::mclapply(1:nrow(counts), function(x) {
      1 - numCell * vec_ybar[x]^2/adjust_nominator[x]
    }, mc.cores = mc_cores))
    vec_daviesp <- unlist(parallel::mclapply(1:nrow(counts), function(x) {
      sparkx_pval(x, vec_ylam, Klam, vec_stat)
    }, mc.cores = mc_cores))
    res_sparkx <- as.data.frame(cbind(vec_stat, vec_daviesp))
  }
  else {
    if (ncol(counts) < 30000) {
      counts <- as.matrix(counts)
    }
    numCell <- nrow(infomat)
    XTX_inv <- solve(crossprod(X_mat, X_mat))
    Xadjust_mat <- crossprod(infomat, X_mat) %*% crossprod(XTX_inv, 
                                                           t(X_mat))
    Xinfomat <- infomat - t(Xadjust_mat)
    info_inv <- solve(crossprod(Xinfomat, Xinfomat))
    kmat_first <- Xinfomat %*% info_inv
    LocDim <- ncol(Xinfomat)
    Klam <- eigen(crossprod(Xinfomat, kmat_first), only.values = T)$values
    res_sparkx_list <- parallel::mclapply(X = 1:nrow(counts), FUN = sparkx.sksg, 
                                expmat = counts, xmat = X_mat, scaleinfo = Xinfomat, 
                                numDim = LocDim, lambda_K = Klam, loc_inv = info_inv, 
                                mc.cores = mc_cores, verbose = verbose)
    res_sparkx <- as.data.frame(do.call(rbind, res_sparkx_list))
  }
  colnames(res_sparkx) <- c("stat", "pval")
  rownames(res_sparkx) <- geneName
  return(res_sparkx)
}

sparkx.sksg <- function (igene, expmat, xmat, scaleinfo, numDim, lambda_K, loc_inv, 
                         verbose = TRUE) {
  if (verbose) {
    message("gene", igene, "\n")
  }
  single_gene <- expmat[igene, ]
  numCell <- length(single_gene)
  XTX_inv <- solve(crossprod(xmat, xmat))
  GTX <- crossprod(single_gene, xmat)
  Gadjust_mat <- GTX %*% tcrossprod(XTX_inv, GTX)
  adj_nominator <- 1/as.vector(crossprod(single_gene, single_gene))
  lambda_G <- as.vector(crossprod(single_gene, single_gene) - 
                          Gadjust_mat) * adj_nominator
  YHL <- single_gene %*% scaleinfo
  scoredavies <- adj_nominator * (YHL %*% loc_inv %*% t(YHL)) * 
    numCell
  Zsort <- sort(lambda_G * lambda_K, decreasing = TRUE)
  results_score <- try(davies(scoredavies, Zsort))
  if (!inherits(results_score, "try-error")) { # class(results_score) != "try-error"
    pout <- results_score$Qq
    if (pout <= 0) {
      pout <- liu(scoredavies, Zsort)
    }
  }
  else {
    pout <- NA
  }
  return(c(scoredavies, pout))
}

sparkx_pval<-function (igene, lambda_G, lambda_K, allstat) 
{
  Zsort <- sort(lambda_G[igene] * lambda_K, decreasing = TRUE)
  results_score <- try(davies(allstat[igene], Zsort))
  if (!inherits(results_score, "try-error")) { # class(results_score) != "try-error"
    pout <- results_score$Qq
    if (pout <= 0) {
      pout <- liu(allstat[igene], Zsort)
    }
  }
  else {
    pout <- NA
  }
  return(pout)
}

transloc_func <- function (coord, lker, transfunc = "gaussian"){
  coord <- scale(coord, scale = F)
  l <- quantile(abs(coord), probs = seq(0.2, 1, by = 0.2))
  if (transfunc == "gaussian") {
    out <- exp(-coord^2/(2 * l[lker]^2))
  }
  else if (transfunc == "cosine") {
    out <- cos(2 * pi * coord/l[lker])
  }
  return(out)
}
ACAT <- function (Pvals, Weights = NULL) {
  if (sum(is.na(Pvals)) > 0) {
    stop("Cannot have NAs in the p-values!")
  }
  if ((sum(Pvals < 0) + sum(Pvals > 1)) > 0) {
    stop("P-values must be between 0 and 1!")
  }
  is.zero <- (sum(Pvals == 0) >= 1)
  is.one <- (sum(Pvals == 1) >= 1)
  if (is.zero && is.one) {
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero) {
    return(0)
  }
  if (is.one) {
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  if (is.null(Weights)) {
    Weights <- rep(1/length(Pvals), length(Pvals))
  }
  else if (length(Weights) != length(Pvals)) {
    stop("The length of weights should be the same as that of the p-values")
  }
  else if (sum(Weights < 0) > 0) {
    stop("All the weights must be positive!")
  }
  else {
    Weights <- Weights/sum(Weights)
  }
  is.small <- (Pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(Weights * tan((0.5 - Pvals) * pi))
  }
  else {
    cct.stat <- sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(Weights[!is.small] * tan((0.5 - 
                                                           Pvals[!is.small]) * pi))
  }
  if (cct.stat > 1e+15) {
    pval <- (1/cct.stat)/pi
  }
  else {
    pval <- 1 - pcauchy(cct.stat)
  }
  return(pval)
}
