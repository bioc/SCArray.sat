\name{SCArray.sat-package}
\alias{SCArray.sat-package}
\alias{SCArray.sat}
\docType{package}
\title{Large-scale single-cell RNA-seq data analysis using GDS files and Seurat}
\description{
    The package extends the Seurat classes and functions to support GDS files
as a DelayedArray backend for data representation. It introduces a new
SCArrayAssay class (derived from the Seurat Assay), which wraps raw counts,
normalized expressions and scaled data matrix based on DelayedMatrix. It is
designed to integrate seamlessly with the SeuratObject and Seurat packages to
provide common data analysis, with the optimized algorithms for GDS data files.
}
\details{
    \tabular{ll}{
        Package: \tab SCArray.sat \cr
        Type: \tab Package \cr
        License: \tab GPL version 3 \cr
    }
}
\author{Xiuwen Zheng}

\keyword{GDS}
\keyword{SingleCell}
