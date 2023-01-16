#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2020-2022    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


# Package-wide variable
.packageEnv <- new.env()


#######################################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")

.plural <- function(num) if (num > 1L) "s" else ""

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

x_check <- function(x, msg) SCArray:::x_check(x, msg)



scGetAssayFromGDS <- function(fn, verbose=TRUE)
{
    # load gds data
    sce <- scExperiment(fn)
    m <- counts(sce)
    # fake matrix
    m0 <- sparseMatrix(i=integer(), j=integer(), x=double(),
        dims=c(nrow(m), ncol(m)))
    rownames(m0) <- rownames(m)
    colnames(m0) <- colnames(m)
    init.meta.features <- data.frame(row.names=rownames(m))
    # output
    rv <- new(Class = "SCArrayAssay",
        counts = m0, data = m0,
        counts2 = m, data2 = m,
        scale.data = new(Class="matrix"),
        meta.features = init.meta.features,
        misc = list())
    return(rv)
}

