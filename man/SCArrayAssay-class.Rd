\name{SCArrayAssay-class}
\alias{SCArrayAssay}
\alias{SCArrayAssay-class}

\docType{class}
\title{GDS-specific Assay Class}
\description{
    The SCArrayAssay class extends the Assay class of Seurat with the new slots
\code{counts2}, \code{data2} and \code{scale.data2} replacing \code{counts},
\code{data} and \code{scale.data}.
}

\section{Slots}{
\describe{
\item{\code{counts2}}{Unnormalized raw counts (dgCMatrix or SC_GDSMatrix),
    replacing \code{Assay@counts}}

\item{\code{data2}}{Normalized expression data (dgCMatrix or SC_GDSMatrix),
    replacing \code{Assay@data}}

\item{\code{scale.data2}}{Scaled expression data (NULL, matrix or SC_GDSMatrix),
    replacing \code{Assay@scale.data}}
}}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{Assay-class}}, \code{\link{GetAssayData}},
    \code{\link{SetAssayData}}
}

\keyword{GDS}
\keyword{SingleCell}
