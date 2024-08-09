#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data analysis using GDS files and Seurat
#
# Copyright (C) 2022-2024    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
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
    x_msg("Calling subset.SCArrayAssay() ...")
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
        warning("NAs passed in the features vector, removing NAs",
            immediate.=TRUE)
        features <- na.omit(features)
    }
    if (all(vapply(list(features, cells), FUN=length, 0L) == dim(x)))
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
        c_scaled <- colnames(s)
        c_scaled <- c_scaled[c_scaled %in% cells]
        c_scaled <- c_scaled[na.omit(match(colnames(x), c_scaled))]
        f_scaled <- rownames(s)
        f_scaled <- f_scaled[f_scaled %in% features]
        if (length(c_scaled) && length(f_scaled))
        {
            x@scale.data2 <- s[f_scaled, c_scaled, drop=FALSE]
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
SetAssayData.SCArrayAssay <- function(object, layer, new.data,
    slot=c('data', 'scale.data', 'counts'), ...)
{
    # check
    CheckDots(...)
    slot <- match.arg(slot)
    if (!IsMatrixEmpty(new.data))
    {
        s <- rownames(new.data)
        if (any(grepl('_', s, fixed=TRUE)))
        {
            warning("Feature names cannot have underscores ('_'), ",
                "replacing with dashes ('-')", call.=FALSE, immediate.=TRUE)
            rownames(new.data) <- gsub('_', '-', s, fixed=TRUE)
        }
        if (ncol(new.data) != ncol(object))
        {
            stop("The new data doesn't have the same number of cells ",
                "as the current data", call.=FALSE)
        }
        num.counts <- nrow(object)
        counts.names <- rownames(object)
        if (slot=='scale.data' && nrow(new.data)>num.counts)
        {
            warning("Adding more features than present in current data",
                call.=FALSE, immediate.=TRUE)
        } else if (slot %in% c('counts', 'data') && nrow(new.data)!=num.counts)
        {
            warning("The new data doesn't have the same number of features ",
                "as the current data", call.=FALSE, immediate.=TRUE)
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
            stop("Attempting to add a different number of cells or features",
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


# Rename cell names
RenameCells.SCArrayAssay <- function(object, new.names=NULL, ...)
{
    CheckDots(...)
    for (data.slot in c("counts2", "data2", "scale.data2"))
    {
        old.data <- slot(object, name=data.slot)
        if (NCOL(old.data) <= 1L) next
        colnames(slot(object, name=data.slot)) <- new.names
    }
    object
}


# Merge a list of SCArrayAssay
.row_merge_mat <- function(all.mat)
{
    all.colnames <- all.rownames <- vector("list", length(all.mat))
    for (i in seq_along(all.mat))
    {
        all.rownames[[i]] <- rownames(all.mat[[i]])
        all.colnames[[i]] <- colnames(all.mat[[i]])
    }
    use.cbind <- all(duplicated(all.rownames)[-1L])
    if (isTRUE(use.cbind))
    {
        new.mat <- do.call(cbind, all.mat)
    } else {
        stop("All matrices should have the same row names for merging.")
    }
    colnames(new.mat) <- make.unique(unlist(all.colnames))
    new.mat
}

merge.SCArrayAssay <- function(x=NULL, y=NULL, add.cell.ids=NULL,
    merge.data=TRUE, ...)
{
    x_msg("Calling merge.SCArrayAssay() ...")
    CheckDots(...)
    assays <- c(x, y)
    if (!is.null(add.cell.ids))
    {
        for (i in seq_along(assays))
            assays[[i]] <- RenameCells(assays[[i]], new.names=add.cell.ids[i])
    }
    # merge the counts
    counts.mats <- lapply(assays, FUN=Seurat:::ValidateDataForMerge,
        slot="counts")
    keys <- vapply(assays, FUN=Key, "")
    merged.counts <- .row_merge_mat(counts.mats)
    combined.assay <- CreateAssayObject2(merged.counts,
        min.cells=-1, min.features=-1)
    if (length(unique(keys)) == 1L)
        Key(combined.assay) <- keys[1L]
    # merge the data
    if (merge.data)
    {
        data.mats <- lapply(assays, FUN=Seurat:::ValidateDataForMerge,
            slot="data")
        merged.data <- .row_merge_mat(data.mats)
        # only keep cells that made it through counts filtering params
        if (!all.equal(colnames(combined.assay), colnames(merged.data)))
        {
            merged.data <- merged.data[, colnames(combined.assay)]
        }
        combined.assay <- SetAssayData(combined.assay,
            slot="data", new.data=merged.data)
    }
    combined.assay
}

merge.Seurat_g <- function(...)
{
    x_msg("Calling merge.Seurat_g() ...")
    new("Seurat_g", NextMethod())
}


####  scGetFiles()  ####

.scget_sc_assay <- function(object, ...)
{
    s <- character()
    if (is(object@counts2, "SC_GDSMatrix"))
        s <- scGetFiles(object@counts2)
    if (is(object@data2, "SC_GDSMatrix"))
        s <- c(s, scGetFiles(object@data2))
    if (is(object@scale.data2, "SC_GDSMatrix"))
        s <- c(s, scGetFiles(object@scale.data2))
    unique(s)
}

.scget_seurat <- function(object, ...)
{
    s <- lapply(Assays(object), function(nm) scGetFiles(object[[nm]]))
    unique(unlist(s))
}

# Get file names for on-disk backend
setMethod("scGetFiles", "Assay", function(object, ...) NULL)
setMethod("scGetFiles", "SCArrayAssay", .scget_sc_assay)
setMethod("scGetFiles", "Seurat", .scget_seurat)


####  scMemory()  ####

.scmemory_sc_assay <- function(x, slot=NULL, ...)
{
    x_msg("Calling scMemory.SCArrayAssay() ...")
    if (is.null(slot))
    {
        # return Assay instead of SCArrayAssay
        counts <- scMemory(x@counts2, ...)
        data <- scMemory(x@data2, ...)
        scale.data <- scMemory(x@scale.data2, ...)
        if (is.null(scale.data)) scale.data <- new("matrix")
        new("Assay",
            counts = counts, data = data, scale.data = scale.data,
            key = x@key, assay.orig = x@assay.orig,
            var.features = x@var.features,
            meta.features = x@meta.features,
            misc = x@misc)
    } else {
        stopifnot(is.character(slot))
        nm <- c("counts", "data", "scale.data")
        s <- setdiff(slot, nm)
        if (length(s))
            stop("'slot' should be ", paste(nm, collapse=", "), ".")
        for (s in slot)
        {
            s <- .redirect_slot(s)
            slot(x, s) <- scMemory(slot(x, s), ...)
        }
        x
    }
}

.scmemory_seurat <- function(x, assay=NULL, slot=NULL, ...)
{
    if (is.null(assay)) assay <- names(x@assays)
    stopifnot(is.character(assay))
    for (nm in assay)
        x[[nm]] <- scMemory(x[[nm]], slot=slot, ...)
    x
}

.scmemory_seurat_g <- function(x, assay=NULL, slot=NULL, ...)
{
    x <- .scmemory_seurat(x, assay=assay, slot=slot, ...)
    if (is(x, "Seurat_g")) x <- as(x, "Seurat")
    x
}

# Load data to memory for SCArrayAssay
setMethod("scMemory", "SCArrayAssay", .scmemory_sc_assay)
# Load data to memory for Seurat
setMethod("scMemory", "Seurat", .scmemory_seurat)
# Load data to memory for Seurat_g
setMethod("scMemory", "Seurat_g", .scmemory_seurat_g)


#######################################################################
# S3 Methods for DelayedMatrix

as.sparse.DelayedMatrix <- function(x, ...)
{
    as(x, "sparseMatrix")
}


#######################################################################
# S3 Methods for SC_GDSMatrix

CheckMatrix.SC_GDSMatrix <- function(object, checks=NULL, ...)
{
    # check efficient row or column
    if (SCArray:::x_type(object) == 2L)
        gd <- rowAutoGrid(object)
    else
        gd <- colAutoGrid(object)
    # block read
    flag <- blockReduce(function(bk, v)
    {
        if (is(bk, "COO_SparseArray")) bk <- as(bk, "SVT_SparseArray")
        v || anyNA(bk) || any(is.infinite(bk))
    }, object, FALSE, grid=gd, as.sparse=NA, verbose=FALSE, BREAKIF=identity)
    if (flag)
        warning("Input matrix contains NA/NaN or infinite values.")
    invisible()
}
