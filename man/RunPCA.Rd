\name{RunPCA}
\alias{RunPCA}
\alias{RunPCA.SCArrayAssay}
\alias{RunPCA.SC_GDSMatrix}

\title{Run PCA}
\description{
    Performs PCA on a Seurat SCArrayAssay or a DelayedMatrix object.
}

\usage{
# RunPCA(object, ...)
\method{RunPCA}{SCArrayAssay}(object, assay=NULL, features=NULL, npcs=50,
    rev.pca=FALSE, weight.by.var=TRUE, verbose=TRUE, ndims.print=seq_len(5),
    nfeatures.print=30, reduction.key="PC_", seed.use=42, ...)
\method{RunPCA}{SC_GDSMatrix}(object, assay=NULL, npcs=50, rev.pca=FALSE,
    weight.by.var=TRUE, verbose=TRUE, ndims.print=seq_len(5),
    nfeatures.print=30, reduction.key="PC_", seed.use=42, approx=TRUE, ...)
}

\arguments{
    \item{object}{input R object (e.g., a SCArrayAssay object)}
    \item{assay}{NULL for using the active assay, or an assay name}
    \item{features}{if \code{NULL}, PCA will be run on the scaled data;
        otherwise, features to compute PCA on}
    \item{npcs}{# of top PCs to be calculated}
    \item{rev.pca}{By default (\code{FALSE}), perform PCA on the cell x gene
        matrix; otherwise, compute it on gene x cell matrix}
    \item{weight.by.var}{if \code{TRUE}, weight the cell embeddings (when
        \code{rev.pca=FALSE}) or the the gene loadings (when
        \code{rev.pca=TRUE}) by the variance of each PC}
    \item{verbose}{if \code{TRUE}, show information}
    \item{ndims.print}{which PCs to print genes for}
    \item{nfeatures.print}{# of genes to print for each PC}
    \item{reduction.key}{dimensional reduction key}
    \item{seed.use}{a random seed; or NULL for not setting a seed internally}
    \item{approx}{if \code{TRUE}, use IrlbaSVD; otherwise use ExactSVD}
    \item{...}{additional arguments passed to specific methods}
}
\details{
    \code{RunPCA()} computes the covariance matrix of genes (if # of genes <=
# of cells) or the cell covariance matrix for the PCA calculation, which can
reduce the times of accessing on-disk data.
}
\value{
    Return a data frame for reduction data (via \code{CreateDimReducObject}).
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[Seurat]{RunPCA}}, \code{\link{CreateDimReducObject}}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scNewAssayGDS(fn)
d <- Seurat::CreateSeuratObject(a)

d <- NormalizeData(d)
d <- FindVariableFeatures(d, nfeatures=250)
d <- ScaleData(d)

d <- RunPCA(d, ndims.print=1:2)
DimPlot(d, reduction="pca")

remove(a, d)
}

\keyword{methods}
\keyword{GDS}
