\name{SCArrayAssay-methods}
\alias{SCArrayAssay-methods}
\alias{GetAssayData}
\alias{GetAssayData.SCArrayAssay}
\alias{SetAssayData}
\alias{SetAssayData.SCArrayAssay}
\alias{subset}
\alias{subset.SCArrayAssay}

\title{SCArrayAssay S3 methods}
\description{
    Gets and sets data in the Seurat Assay object.
}
\usage{
\method{GetAssayData}{SCArrayAssay}(object,
    slot=c("data", "scale.data", "counts"), ...)
\method{SetAssayData}{SCArrayAssay}(object, layer, new.data,
    slot=c('data', 'scale.data', 'counts'), ...)

\method{subset}{SCArrayAssay}(x, cells=NULL, features=NULL, ...)
}
\arguments{
    \item{object, x}{a SCArrayAssay object (inherited from SeuratObject::Assay)}
    \item{layer}{layer}
    \item{new.data}{a new data matrix (dgCMatrix or SC_GDSMatrix)}
    \item{slot}{data matrix in the Assay object, "data" is used by default}
    \item{cells}{names or indices for selected cells}
    \item{features}{names or indices for selected features}
    \item{...}{further arguments to be passed to or from other methods}
}
\value{
    Return a data matrix or an instance of \code{\link{SCArrayAssay}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[SeuratObject]{Assay}}
}

\keyword{GDS}
\keyword{SingleCell}
