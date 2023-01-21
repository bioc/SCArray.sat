\name{scGetAssayFromGDS}
\alias{scGetAssayFromGDS}
\title{Get Assay object}
\description{
    Gets Seurat Assay object.
}
\usage{
scGetAssayFromGDS(fn, verbose=TRUE)
}
\arguments{
    \item{fn}{a file name for the GDS file}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return an instance of \code{\link{SCArrayAssay}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{SCArrayAssay}}
}

\examples{
# raw count data in a GDS file
fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")

assay <- scGetAssayFromGDS(fn)
assay
class(assay)

rm(assay)
}

\keyword{GDS}
\keyword{SingleCell}
