\name{readscRNAseq}
\alias{readscRNAseq}
\title{Read the scRNAseq data measured on scRNA sequencing platform}
\description{
  Read the single cell RNA sequencing data measured on scRNA sequencing platform  as a Seurat object.
}
\usage{
  readscRNAseq(mtx, cells, features, ...)
}
\arguments{
  \item{mtx}{a string, ame or remote URL of the mtx file}
  \item{cells}{a string, Name or remote URL of the cells/barcodes file}
  \item{features}{a string, Name or remote URL of the features/genes file}
  \item{...}{the arguments passing to \link{ReadMtx}}
}
\details{
  Nothing
}
\value{
  return a Seurat object including expression matrix.
}
\references{
  None
}
\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  None
}
\examples{
\dontrun{
  ### set the file directory, then read it.
  seu <- readscRNAseq(mtx="GSM3755564_16_Liver_Treg_matrix.mtx.gz",
                       features='GSM3755564_16_Liver_Treg_genes.tsv.gz',
                       cells='GSM3755564_16_Liver_Treg_barcodes.tsv.gz' )
  seu
  }
}
\keyword{Input of scRNAseq data}
