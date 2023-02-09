\name{scNewAssayGDS}
\alias{scNewAssayGDS}
\title{Create Assay Object}
\description{
    Creates a new Seurat Assay object (SCArrayAssay) from a GDS file.
}
\usage{
scNewAssayGDS(gdsfile, name="counts", key="rna_", row_data=TRUE, check=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a file name for the GDS file, or a \code{SCArrayFileClass}
        object}
    \item{name}{characters for the name of data matrix in the GDS file;
        if \code{NA_character_}, to use the first assay}
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
    \code{\link{scExperiment}}, \code{\link{scNewSeuratGDS}}
}

\examples{
# raw count data in a GDS file
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scNewAssayGDS(fn)
a
class(a)

d <- Seurat::CreateSeuratObject(a)
d

rm(a, d)
}

\keyword{GDS}
\keyword{SingleCell}
