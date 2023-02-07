\name{ScaleData}
\alias{ScaleData}
\alias{ScaleData.SC_GDSMatrix}

\title{Scale and Center the Data}
\description{
    Scales and centers features or residuals in the dataset.
}

\usage{
# ScaleData(object, ...)
\method{ScaleData}{SC_GDSMatrix}(object, features=NULL, vars.to.regress=NULL,
    latent.data=NULL, split.by=NULL, model.use='linear', use.umi=FALSE,
    do.scale=TRUE, do.center=TRUE, scale.max=10, block.size=1000,
    min.cells.to.block=3000, verbose=TRUE, use_gds=TRUE, rm_tmpfile=TRUE, ...)
}

\arguments{
    \item{object}{input R object (e.g., a SC_GDSMatrix object)}
    \item{features}{if \code{NULL}, to use the variable features (found via
        \code{FindVariableFeatures()}); or features names to scale/center}
    \item{vars.to.regress}{NULL or variable names to regress out}
    \item{latent.data}{NULL or a \code{data.frame} to regress out the
        covariates}
    \item{split.by}{variable name in the metadata, or a vector or factor
        defining grouping of cells}
    \item{model.use}{regression model: \code{"linear"} (default),
        \code{"poisson"} or \code{"negbinom"}}
    \item{use.umi}{only applicable when the covariates are given in
        \code{vars.to.regress} for regression; default is \code{FALSE} for
        linear regression, \code{TRUE} for negbinom and poisson models}
    \item{do.scale}{if \code{TRUE}, scale the data}
    \item{do.center}{if \code{TRUE}, center the data}
    \item{scale.max}{max value in the resulting scaled data; see
        \code{ScaleData.Seurat} for more details}
    \item{block.size}{not used}
    \item{min.cells.to.block}{not used}
    \item{verbose}{if \code{TRUE}, show information}
    \item{use_gds}{if \code{TRUE}, use \code{SC_GDSMatrix} for the scaled data;
        if \code{FALSE}, to use a dense in-memory scaled matrix; or a GDS file
        name for storing the resulting data matrix; see details}
    \item{rm_tmpfile}{if \code{TRUE}, remove any temporary GDS file after
        the calculation; the temporary file will be created when calculating
        the residuals}
    \item{...}{additional arguments passed to specific methods}
}
\details{
    \code{ScaleData()} stores the scaled data in a GDS file when
\code{use_gds=TRUE} or an output GDS file name is given via \code{use_gds}.
When \code{vars.to.regress} and \code{split.by} are both \code{NULL}, an output
GDS file is not needed, since the resulting DelayedMatrix can be represented
as common operations on the count matrix. If \code{use_gds=TRUE}, an output
file name "_scale_data.gds" will be used if it does not exists, or
"_scale_data2.gds" (if not exists), "_scale_data3.gds" and so on.
    When \code{vars.to.regress} are given, a temporary GDS file (e.g., 
"_temp_scale_data.gds", \code{use_gds} with a prefix "_temp") will be created
to store the residuals before scaling. This temporary file will be deleted
after the calculation when \code{rm_tmpfile=TRUE}.
}
\value{
    Returns a SC_GDSMatrix matrix if \code{use_gds=TRUE} or \code{use_gds} is an
output file name, otherwise returns an in-memory matrix.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link[Seurat]{ScaleData}}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

a <- scGetAssayGDS(fn)
d <- Seurat::CreateSeuratObject(a)

d <- NormalizeData(d)
d <- FindVariableFeatures(d, nfeatures=50)
d <- ScaleData(d)

GetAssayData(d, "scale.data")  # DelayedMatrix

# scale with split.by
ss <- rep(c(TRUE, FALSE), length.out=ncol(d))
d <- ScaleData(d, split.by=ss)

fn <- scGetFiles(d)
fn[2L]    # the file name storing scaled data

remove(a, d)
unlink(fn[grepl("^_scale", fn)], force=TRUE)
}

\keyword{methods}
\keyword{GDS}
