\name{scGetAssayGDS}
\alias{scGetAssayGDS}
\title{Get Assay object}
\description{
    Gets Seurat Assay object.
}
\usage{
scGetAssayGDS(gdsfn, row_data=TRUE, check=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsfn}{a file name for the GDS file}
    \item{row_data}{if \code{TRUE}, add \code{rowData()} to the feature-level
        meta data of the Seurat Assay}
    \item{check}{if \code{TRUE}, check the feature names}
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
fn <- system.file("extdata", "example.gds", package="SCArray")

assay <- scGetAssayGDS(fn)
assay
class(assay)

rm(assay)
}

\keyword{GDS}
\keyword{SingleCell}
