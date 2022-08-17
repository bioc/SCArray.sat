#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2022    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################
# Internal functions

.redirect_slot <- function(nm)
{
    switch(nm, data="data2", counts="counts2", nm)
}



#######################################################################
# Methods for SCArrayAssay

GetAssayData.SCArrayAssay <- function(object,
    slot=c("data", "scale.data", "counts"), ...)
{
    CheckDots(...)
    slot <- match.arg(slot)
    # output
    slot(object, .redirect_slot(slot))
}


SetAssayData.SCArrayAssay <- function(object,
    slot=c('data', 'scale.data', 'counts'), new.data, ...)
{
  CheckDots(...)
  # print(dim(object))
  # print(dim(new.data))
  slot <- slot[1]
  slot <- match.arg(arg = slot)
  if (!IsMatrixEmpty(x = new.data)) {
    if (ncol(x = new.data) != ncol(x = object)) {
      stop(
        "The new data doesn't have the same number of cells as the current data",
        call. = FALSE
      )
    }
    num.counts <- nrow(x = object)
    counts.names <- rownames(x = object)
    if (slot == 'scale.data' && nrow(x = new.data) > num.counts) {
      warning(
        "Adding more features than present in current data",
        call. = FALSE,
        immediate. = TRUE
      )
    } else if (slot %in% c('counts', 'data') && nrow(x = new.data) != num.counts) {
      warning(
        "The new data doesn't have the same number of features as the current data",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!all(rownames(x = new.data) %in% counts.names)) {
      warning(
        "Adding features not currently present in the object",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    new.features <- na.omit(object = match(
      x = counts.names,
      table = rownames(x = new.data)
    ))
    new.cells <- colnames(x = new.data)
    if (!all(new.cells %in% colnames(x = object))) {
      stop(
        "All cell names must match current cell names",
        call. = FALSE
      )
    }
    new.data <- new.data[new.features, colnames(object), drop = FALSE]
    # print(dim(new.data))
    # print(dim(object))
    if (slot %in% c('counts', 'data') && !all(dim(new.data) == dim(object))) {
      stop(
        "Attempting to add a different number of cells and/or features",
        call. = FALSE
      )
    }
  }
  if (!is.vector(x = rownames(x = new.data))) {
    rownames(x = new.data) <- as.vector(x = rownames(x = new.data))
  }
  if (!is.vector(x = colnames(x = new.data))) {
    colnames(x = new.data) <- as.vector(x = colnames(x = new.data))
  }
  slot(object = object, name = .redirect_slot(slot)) <- new.data
  return(object)
}


####  Methods -- NormalizeData()  ####

.log_norm <- function(mat, scale.factor, verbose)
{
    stopifnot(inherits(mat, "DelayedArray"))
    s <- scale.factor / DelayedArray::colSums(mat)
    m <- log1p(DelayedArray::sweep(mat, 2L, s, `*`))
    as(m, "SC_GDSMatrix")
}

NormalizeData.DelayedMatrix <- function(object,
    normalization.method="LogNormalize", scale.factor=10000, margin=1,
    verbose=TRUE, ...)
{
    # check
    CheckDots(...)
    if (!is.null(normalization.method))
    {
        stopifnot(is.character(normalization.method),
            length(normalization.method)==1L)
    }
    stopifnot(margin==1)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is.null(normalization.method)) return(object)
    switch(normalization.method,
        "LogNormalize" = .log_norm(object, scale.factor, verbose),
        # "CLR", "RC"
        stop("Unknown normalization method: ", normalization.method)
    )
}


####  Methods -- FindVariableFeatures()  ####

.row_var_std <- function(mat, mu, sd, vmax, verbose)
{
    stopifnot(inherits(mat, "DelayedArray"))
    # block read
    v <- blockReduce(function(bk, v, mu, sd, vmax)
    {
        b <- pmin(vmax, (bk - mu) / sd)^2L
        dim(b) <- dim(bk)
        v + rowSums(b)
    }, mat, init=0, grid=colAutoGrid(mat), mu=mu, sd=sd, vmax=vmax)
    v[!is.finite(v)] <- 0
    # output
    v / (ncol(mat)-1L)
}

FindVariableFeatures.DelayedMatrix <- function(object,
    selection.method="vst", loess.span=0.3, clip.max="auto",
    mean.function=FastExpMean, dispersion.function=FastLogVMR,
    num.bin=20, binning.method="equal_width", verbose=TRUE, ...)
{
    # check
    CheckDots(...)

    if (selection.method == "vst")
    {
        if (clip.max == "auto")
            clip.max <- sqrt(ncol(object))
        if (verbose)
            cat("Calculating gene variances\n")
        hvf.info <- data.frame(mean=rowMeans(object))
        hvf.info$variance <- rowVars(object)
        hvf.info$variance[is.na(hvf.info$variance)] <- 0
        hvf.info$variance.expected <- 0
        hvf.info$variance.standardized <- 0
        not.const <- hvf.info$variance > 0
        fit <- loess(log10(variance) ~ log10(mean),
            hvf.info[not.const, ], span=loess.span)
        hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted

        # get variance after feature standardization
        if (verbose)
            cat("Calculating feature variances of standardized and clipped values\n")
        hvf.info$variance.standardized <- .row_var_std(
            object, hvf.info$mean, sqrt(hvf.info$variance.expected),
            clip.max, verbose)

        colnames(hvf.info) <- paste0('vst.', colnames(hvf.info))
    } else {
        stop("selection.method!=vst, not implemented yet.")
    }
  return(hvf.info)
}







