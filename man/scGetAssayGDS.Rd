\name{scGetAssayGDS}
\alias{scGetAssayGDS}
\title{Get Assay object}
\description{
    Gets Seurat Assay object.
}
\usage{
scGetAssayGDS(gdsfile, name="counts", key="rna_", row_data=TRUE, check=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a file name for the GDS file, or a \code{SCArrayFileClass}
        object}
    \item{name}{characters for the name of data matrix in the GDS file}
    \item{key}{characters for the Assay key}
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
    \code{\link{SCArrayAssay}}, \code{\link{SCArrayFileClass}},
    \code{\link{scExperiment}}
}

\examples{
# raw count data in a GDS file
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scGetAssayGDS(fn)
a
class(a)

rm(a)
}

\keyword{GDS}
\keyword{SingleCell}
