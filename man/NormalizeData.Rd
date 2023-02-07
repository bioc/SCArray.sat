\name{NormalizeData}
\alias{NormalizeData}
\alias{NormalizeData.SC_GDSMatrix}

\title{Normalize Count Data}
\description{
    Normalizes the count data in the Seurat assay.
}

\usage{
# NormalizeData(object, ...)
\method{NormalizeData}{SC_GDSMatrix}(object,
    normalization.method="LogNormalize", scale.factor=1e4, margin=1,
    verbose=TRUE, ...)
}

\arguments{
    \item{object}{input R object (e.g., a SC_GDSMatrix object)}
    \item{normalization.method}{"LogNormalize", "CLR" or "RC";
        see \code{NormalizeData.Seurat} for more details}
    \item{scale.factor}{the scale factor for cell-level normalization}
    \item{margin}{only applicable when \code{normalization.method="CLR"},
        normalize across features (\code{margin=1}) or cells (\code{margin=2})}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to specific methods}
}
\details{
    \code{NormalizeData()} does not store the normalized data in a GDS file,
since the calculation is "delayed" until it is needed.
}
\value{
    Returns a SC_GDSMatrix matrix.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[Seurat]{NormalizeData}}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scGetAssayGDS(fn)
d <- Seurat::CreateSeuratObject(a)

d <- NormalizeData(d)

remove(a, d)
}

\keyword{methods}
\keyword{GDS}
