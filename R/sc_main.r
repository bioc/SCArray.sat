#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data analysis using GDS files and Seurat
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

# if getOption("SCArray.verbose")=TRUE, show message
x_check <- function(x, msg) SCArray:::x_check(x, msg)

# if getOption("SCArray.verbose")=TRUE, show message
x_msg <- function(msg) SCArray:::x_check(NULL, msg)



#######################################################################

scGetAssayGDS <- function(gdsfn, row_data=TRUE, check=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfn), length(gdsfn)==1L)
    stopifnot(is.logical(check), length(check)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(row_data) || is.data.frame(row_data))
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
    # meta data
    meta_data <- data.frame(row.names=rownames(m))
    if (isTRUE(row_data))
    {
        v <- rowData(sce)
        if (!is.null(v) && ncol(v)>0L)
        {
            v <- as.data.frame(v)
            if (!identical(rownames(m), rownames(v)))
                stop("The rownames of 'rowData()' should be the same as 'count' matrix.")
            meta_data <- v
        }
    }
    # output
    new(Class = "SCArrayAssay",
        counts2 = m, data2 = m, scale.data2 = NULL,
        meta.features = meta_data, misc = list())
}

