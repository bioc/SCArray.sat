\name{SCArray.sat-package}
\alias{SCArray.sat-package}
\docType{package}
\title{Large-scale single-cell RNA-seq data analysis using GDS files and Seurat}
\description{
    The package extends Seurat classes and functions to support GDS files as a
DelayedArray backend for data representation. It provides the SCArrayAssay
class (inherited from the Seurat Assay) to wrap the DelayedMatrix-based raw
counts, normalized expression values and scaled data matrix. It is designed to
integrate seamlessly with the SeuratObject and Seurat packages to provide the
downstream data analysis, with the optimized algorithms for GDS data files.
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
