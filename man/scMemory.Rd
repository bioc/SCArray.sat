\name{scMemory}
\alias{scMemory}
\alias{scMemory,SCArrayAssay-method}
\alias{scMemory,Seurat-method}

\title{Load Data to Memory}
\description{
    Loads the internal data to memory for any on-disk object.
}

\usage{
scMemory(x, ...)
\S4method{scMemory}{SCArrayAssay}(x, slot=NULL, ...)
\S4method{scMemory}{Seurat}(x, assay=NULL, slot=NULL, ...)
}

\arguments{
    \item{x}{input R object (e.g., a Seurat object)}
    \item{assay}{NULL for using the active assay, or a list of assay names}
    \item{slot}{NULL for all \code{"counts"}, \code{"data"} and
        \code{"scale.data"}; or a character vector including \code{"counts"},
        \code{"data"} or \code{"scale.data"}; see details}
    \item{...}{additional arguments passed to specific methods}
}
\details{
    If \code{slot=NULL}, return a \code{Assay} object instead of
\code{SCArrayAssay} object, so it can downgrade a \code{SCArrayAssay} object to
a \code{Assay} object.
}
\value{
    Return an object (it maybe a different type from \code{class(x)}).
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[SCArray]{scMemory}}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

d1 <- scNewSeuratGDS(fn)
is(GetAssay(d1))

d2 <- scMemory(d1)
is(GetAssay(d2))

remove(a, d1, d2)
}

\keyword{methods}
\keyword{GDS}
