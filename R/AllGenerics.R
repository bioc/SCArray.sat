#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################
# S3 Methods for SCArrayAssay

subset.SCArrayAssay <- function(x, cells=NULL, features=NULL, ...)
{
    CheckDots(...)
    # check cells
    if (is.null(cells)) cells <- colnames(x)
    if (all(is.na(cells)))
    {
        cells <- colnames(x)
    } else if (anyNA(cells))
    {
        warning("NAs passed in cells vector, removing NAs", immediate.=TRUE)
        cells <- na.omit(cells)
    }
    # check features
    if (is.null(features)) features <- rownames(x)
    if (all(is.na(features)))
    {
        features <- rownames(x)
    } else if (anyNA(features))
    {
        warning("NAs passed in the features vector, removing NAs", immediate.=TRUE)
        features <- na.omit(features)
    }
    if (all(sapply(list(features, cells), FUN=length) == dim(x)))
        return(x)
    if (is.numeric(features)) features <- rownames(x)[features]
    features <- gsub(paste0('^', Key(x)), '', features)
    features <- intersect(features, rownames(x))
    if (length(features) == 0L)
        stop("Cannot find features provided")

    s <- GetAssayData(x, 'counts')
    if (ncol(s) == ncol(x))
        slot(x, "counts2") <- s[features, cells, drop=FALSE]
    slot(x, "data2") <- GetAssayData(x, "data")[features, cells, drop=FALSE]

    s <- GetAssayData(x, "scale.data")
    if (!is.null(s))
    {
        cells.scaled <- colnames(s)
        cells.scaled <- cells.scaled[cells.scaled %in% cells]
        cells.scaled <- cells.scaled[na.omit(match(colnames(x), cells.scaled))]
        features.scaled <- rownames(s)
        features.scaled <- features.scaled[features.scaled %in% features]
        if (length(cells.scaled)>0L && length(features.scaled)>0L)
        {
            slot(x, "scale.data2") <- s[features.scaled, cells.scaled, drop=FALSE]
        } else {
            slot(x, "scale.data2") <- NULL
        }
    }

    s <- VariableFeatures(x)
    VariableFeatures(x) <- s[s %in% features]
    slot(x, 'meta.features') <- x[[]][features, , drop=FALSE]

    # output
    return(x)
}

