\name{scGetFiles}
\alias{scGetFiles}
\alias{scGetFiles,Assay-method}
\alias{scGetFiles,SCArrayAssay-method}
\alias{scGetFiles,Seurat-method}

\title{File names for on-disk backend}
\description{
    Get a list of file names for DelayedArray with an on-disk backend.
}

\usage{
scGetFiles(object, ...)
\S4method{scGetFiles}{Assay}(object, ...)
\S4method{scGetFiles}{SCArrayAssay}(object, ...)
\S4method{scGetFiles}{Seurat}(object, ...)
}

\arguments{
    \item{object}{input R object (e.g., a Seurat object)}
    \item{...}{additional arguments passed to specific methods}
}
\value{
    Return a character vector storing file names.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[SCArray]{scGetFiles}}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scNewAssayGDS(fn)
d <- Seurat::CreateSeuratObject(a)

scGetFiles(a)
scGetFiles(d)

remove(a, d)
}

\keyword{methods}
\keyword{GDS}
