\name{CreateAssayObject2}
\alias{CreateAssayObject2}
\title{Create an Assay object}
\description{
    Create an Assay5 object with DelayedMatrix layers from counts or
prenormalized data.
}
\usage{
CreateAssayObject2(counts, data, min.cells=0, min.features=0,
    key=NULL, check.matrix=FALSE, ...)
}
\arguments{
    \item{counts}{Unnormalized raw counts (matrix, dgCMatrix or DelayedMatrix)}
    \item{data}{Prenormalized data (matrix, dgCMatrix or DelayedMatrix)}
    \item{min.cells}{if > 0, a lower cutoff for filtering cells}
    \item{min.features}{if > 0, a lower cutoff for filtering features}
    \item{check.matrix}{Check counts matrix for NA, NaN, Inf, and non-integer
        values}
    \item{key}{Key name for the assay}
    \item{...}{Arguments passed to \code{\link[SeuratObject]{as.sparse}} when \code{counts}
        or \code{data} is matrix or dgCMatrix}
}
\details{
    Similar to \code{SeuratObject::CreateAssayObject()}, except allowing
DelayedMatrix \code{counts} or \code{data}, and returning an \code{Assay5}
object. \code{counts} and \code{data} should not be provided at the same time.
}
\value{
    Return an instance of \code{\link[SeuratObject:Assay5-class]{Assay5}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[SeuratObject]{CreateAssayObject}}, \code{\link{scGetFiles}}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

x <- scArray(fn, "counts")
colnames(x) <- paste0("c", 1:ncol(x))
rownames(x) <- paste0("g", 1:nrow(x))
x

a <- CreateAssayObject2(x)
a

scGetFiles(x)
scGetFiles(a)

remove(x, a)
}

\keyword{GDS}
\keyword{SingleCell}
