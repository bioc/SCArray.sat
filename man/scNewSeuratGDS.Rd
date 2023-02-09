\name{scNewSeuratGDS}
\alias{scNewSeuratGDS}
\title{Create Seurat Object}
\description{
    Creates a new Seurat object from a GDS file.
}
\usage{
scNewSeuratGDS(gdsfile, assay.name=NULL, key=c(counts="rna_"), row_data=TRUE,
    col_data=TRUE, check=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a file name for the GDS file, or a \code{SCArrayFileClass}
        object}
    \item{assay.name}{characters for the name of data matrix in the GDS file;
        if \code{NULL}, to use all of the assays}
    \item{key}{a character vector for an assay key map, where its names are
        the GDS node names}
    \item{row_data}{if \code{TRUE}, add \code{rowData()} to the feature-level
        meta data of the first Seurat Assay}
    \item{col_data}{if \code{TRUE}, add \code{colData()} to the cell-level
        meta data of the Seurat object}
    \item{check}{if \code{TRUE}, check the feature names}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    \code{"counts"} must be in the input GDS file and it is used as the raw count data in
the active Seurat assay. If \code{"logcounts"} exists, it is used as normalized data
associated with \code{"counts"}. If there are other data matrices in the GDS file,
they will be added to the assay list.
}
\value{
    Return an instance of \code{\link{Seurat}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{SCArrayAssay}}, \code{\link{SCArrayFileClass}},
    \code{\link{scExperiment}}, \code{\link{scNewAssayGDS}}
}

\examples{
# raw count data in a GDS file
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scNewAssayGDS(fn)
a
class(a)

rm(a)
}

\keyword{GDS}
\keyword{SingleCell}
