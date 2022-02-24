
mbicPlot <- function(seu, criteria="MBIC"){
  
  if (!inherits(seu, "Seurat"))
    stop("seu must be a Seurat object!")
  
  assy <- DefaultAssay(seu)
  
  if(is.null(seu[[assy]]@misc$icMat)) stop("There is no MBIC-related information in 'seu' object!")
  # library(ggplot2)
  
  
  icMat <- as.data.frame(seu[[assy]]@misc$icMat)
  if(nrow(icMat)<2) stop("mbicPlot: there are less than two candidates in K!")
  
  criteria <- toupper(criteria)
  ggplot(data=icMat,
         aes_string(x='K', y=criteria)) + geom_line(size=1) + cowplot::theme_cowplot() + 
    ylab(paste0(criteria, " value"))
}

spatialPlotClusters <- function(seu){
  
  if (!inherits(seu, "Seurat"))
    stop("seu must be a Seurat object!")
  # require(ggplot2)
  cols <- brewer.pal(12, "Set3")
  K <- length(unique(seu$spatial.drsc.cluster))
  dat <- data.frame(row=seu$row, col=seu$col, clusters=factor(seu$spatial.drsc.cluster))
  clusters <- dat$clusters
  p1 <- ggplot(dat, aes(x=row, y=col, color=clusters)) +
    geom_point(size = 3, alpha=0.7) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text( size=16),
          legend.title = element_text(face='bold', size=18))
  if(K <= 12){
    p1 <- p1 +  scale_color_manual(values=cols[1:K])
  }
  return(p1)
}


drscPlot <- function(seu, dims=1:5, visu.method='tSNE',...){
  
  if (!inherits(seu, "Seurat"))
    stop("seu must be a Seurat object!")
  if(is.null(seu@reductions$'dr-sc')) stop("There is no 'dr-sc' component in slot 'reduction' of 'seu',
                            please perform DR-SC dimension reduction using 'DR.SC' function first!")
  # require(Seurat)
  if(visu.method=='tSNE'){
    seu <- RunTSNE(seu, reduction="dr-sc", dims=dims,verbose = F)
  }else if(visu.method=='UMAP'){
    seu <- RunUMAP(seu, reduction="dr-sc", dims=dims,verbose = F)
  }
  
  p1 <- DimPlot(seu, ...) 
  return(p1)
  
}



