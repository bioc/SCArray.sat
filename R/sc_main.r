#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2022-2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
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
    # output
    rv <- new(Class = "SCArrayAssay",
        counts2 = m, data2 = m, scale.data2 = NULL,
        meta.features = data.frame(row.names=rownames(m)),
        misc = list())
    return(rv)
}

