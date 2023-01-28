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


#######################################################################
# Internal functions

.redirect_slot <- function(nm)
{
    stopifnot(is.character(nm), length(nm)==1L)
    switch(nm, data="data2", counts="counts2", scale.data="scale.data2", nm)
}


#######################################################################
# S3 Methods for SCArrayAssay

# Get the subset of a SCArrayAssay object
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

    s <- GetAssayData(x, "counts")
    if (ncol(s) == ncol(x))
        x@counts2 <- s[features, cells, drop=FALSE]
    x@data2 <- GetAssayData(x, "data")[features, cells, drop=FALSE]

    s <- GetAssayData(x, "scale.data")
    if (!is.null(s))
    {
        cells.scaled <- colnames(s)
        cells.scaled <- cells.scaled[cells.scaled %in% cells]
        cells.scaled <- cells.scaled[na.omit(match(colnames(x), cells.scaled))]
        features.scaled <- rownames(s)
        features.scaled <- features.scaled[features.scaled %in% features]
        if (length(cells.scaled) && length(features.scaled))
        {
            x@scale.data2 <- s[features.scaled, cells.scaled, drop=FALSE]
        } else {
            x@scale.data2 <- NULL
        }
    }

    s <- VariableFeatures(x)
    VariableFeatures(x) <- s[s %in% features]
    x@meta.features <- x[[]][features, , drop=FALSE]

    # output
    return(x)
}


# Get data matrix from a SCArrayAssay object
GetAssayData.SCArrayAssay <- function(object,
    slot=c("data", "scale.data", "counts"), ...)
{
    CheckDots(...)
    slot <- match.arg(slot)
    slot(object, .redirect_slot(slot))  # output
}


# Set data matrix in a SCArrayAssay object
SetAssayData.SCArrayAssay <- function(object,
    slot=c('data', 'scale.data', 'counts'), new.data, ...)
{
    # check
    CheckDots(...)
    slot <- match.arg(slot)
    if (!IsMatrixEmpty(new.data))
    {
        s <- rownames(new.data)
        if (any(grepl('_', s, fixed=TRUE)))
        {
            warning(
                "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
                call.=FALSE, immediate.=TRUE)
            rownames(new.data) <- gsub('_', '-', s, fixed=TRUE)
        }
        if (ncol(new.data) != ncol(object))
        {
            stop("The new data doesn't have the same number of cells as the current data",
                call.=FALSE)
        }
        num.counts <- nrow(object)
        counts.names <- rownames(object)
        if (slot=='scale.data' && nrow(new.data)>num.counts)
        {
            warning("Adding more features than present in current data",
                call.=FALSE, immediate.=TRUE)
        } else if (slot %in% c('counts', 'data') && nrow(new.data)!=num.counts)
        {
            warning(
                "The new data doesn't have the same number of features as the current data",
                call.=FALSE, immediate.=TRUE)
        }
        if (!all(rownames(new.data) %in% counts.names))
        {
            warning("Adding features not currently present in the object",
                call.=FALSE, immediate.=TRUE)
        }
        new.features <- na.omit(match(counts.names, rownames(new.data)))
        new.cells <- colnames(new.data)
        if (!all(new.cells %in% colnames(object)))
        {
            stop("All cell names must match current cell names", call.=FALSE)
        }
        new.data <- new.data[new.features, colnames(object), drop=FALSE]
        if (slot %in% c('counts', 'data') && !all(dim(new.data)==dim(object)))
        {
            stop("Attempting to add a different number of cells and/or features",
                call.=FALSE)
        }
    }
    if (!is.vector(rownames(new.data)))
        rownames(new.data) <- as.vector(rownames(new.data))
    if (!is.vector(colnames(new.data)))
        colnames(new.data) <- as.vector(colnames(new.data))
    # set
    slot(object, .redirect_slot(slot)) <- new.data
    return(object)
}

