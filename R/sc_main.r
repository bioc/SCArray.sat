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



scGetAssayGDS <- function(gdsfn, check=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfn), length(gdsfn)==1L)
    stopifnot(is.logical(check), length(check)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    verbose <- isTRUE(verbose)
    # load gds data
    if (verbose) .cat("Input: ", gdsfn)
    sce <- scExperiment(gdsfn)
    m <- counts(sce)
    if (verbose)
        .cat("    counts: ", NROW(m), " x ", NCOL(m))
    # check feature IDs
    s <- rownames(m)
    if (isTRUE(check) && any(grepl('_', s)))
    {
        warning(
            "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
            immediate. = TRUE
        )
        rownames(m) <- gsub('_', '-', s)
    }
    # output
    new(Class = "SCArrayAssay",
        counts2 = m, data2 = m, scale.data2 = NULL,
        meta.features = data.frame(row.names=rownames(m)),
        misc = list())
}

