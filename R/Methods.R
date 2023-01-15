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
    stopifnot(is.character(nm), length(nm)==1L)
    switch(nm, data="data2", counts="counts2", scale.data="scale.data2", nm)
}



#######################################################################
# Methods for SCArrayAssay

GetAssayData.SCArrayAssay <- function(object,
    slot=c("data", "scale.data", "counts"), ...)
{
    CheckDots(...)
    slot <- match.arg(slot)
    slot(object, .redirect_slot(slot))  # output
}


SetAssayData.SCArrayAssay <- function(object,
    slot=c('data', 'scale.data', 'counts'), new.data, ...)
{
    # check
    CheckDots(...)
    slot <- match.arg(slot)
    if (!IsMatrixEmpty(new.data))
    {
        if (any(grepl('_', rownames(new.data), fixed=TRUE)))
        {
            warning("Feature names cannot have underscores ('_'), replacing with dashes ('-')",
                call.=FALSE, immediate.=TRUE)
            rownames(new.data) <- gsub('_', '-', rownames(new.data), fixed=TRUE)
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
            warning("The new data doesn't have the same number of features as the current data",
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
    {
        rownames(new.data) <- as.vector(rownames(new.data))
    }
    if (!is.vector(colnames(new.data)))
    {
        colnames(new.data) <- as.vector(colnames(new.data))
    }
    # set
    slot(object, .redirect_slot(slot)) <- new.data
    return(object)
}



####  Methods -- CreateSeuratObject()  ####

CreateSeuratObject.SCArrayAssay <- function(counts, project='SeuratProject',
    assay='RNA', names.field=1, names.delim='_', meta.data=NULL, ...)
{
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = counts)))) {
      warning("Some cells in meta.data not present in provided counts matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data), y = colnames(x = counts)), , drop = FALSE]
    }
    if (is.data.frame(x = meta.data)) {
      new.meta.data <- data.frame(row.names = colnames(x = counts))
      for (ii in 1:ncol(x = meta.data)) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
  # Check assay key
  if (!length(x = Key(object = counts)) || !nchar(x = Key(object = counts))) {
    Key(object = counts) <- SeuratObject:::UpdateKey(key = tolower(x = assay))
  }
  assay.list <- list(counts)
  names(x = assay.list) <- assay
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    X = colnames(x = counts),
    FUN = SeuratObject:::ExtractField,
    field = names.field,
    delim = names.delim
  )))
  if (any(is.na(x = idents))) {
    warning(
      "Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = counts))
  }
  names(x = idents) <- colnames(x = counts)
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    meta.data = data.frame(row.names = colnames(x = counts)),
    active.assay = assay,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'SeuratObject')
  )
  object[['orig.ident']] <- idents
  # Calculate nCount and nFeature
  # n.calc <- CalcN(object = counts)
  # if (!is.null(x = n.calc)) {
  #   names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
  #   object[[names(x = n.calc)]] <- n.calc
  # }
  # Add metadata
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}



####  Methods -- NormalizeData()  ####

.log_norm <- function(mat, scale.factor, verbose)
{
    stopifnot(is(mat, "DelayedArray"))
    s <- scale.factor / DelayedArray::colSums(mat)
    m <- log1p(DelayedArray::sweep(mat, 2L, s, `*`))
    # as(m, "SC_GDSMatrix")
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
        stop("Unknown or not implemented normalization method: ",
            normalization.method)
    )
}


####  Methods -- ScaleData()  ####

x_row_scale <- function(mat, center=TRUE, scale=TRUE)
{
    if (center) r_m <- rowMeans(mat)
    if (scale)
    {
        if (center)
            rsd <- rowSds(mat, center=r_m)
        else
            rsd <- rowSds(mat, center=rep(0, nrow(mat)))
    }
    if (center) mat <- mat - r_m
    if (scale) mat <- mat/rsd
    return(mat)
}


ScaleData.DelayedMatrix <- function(object, features=NULL, vars.to.regress=NULL,
    latent.data=NULL, split.by=NULL, model.use='linear', use.umi=FALSE,
    do.scale=TRUE, do.center=TRUE, scale.max=10, block.size=1000,
    min.cells.to.block=3000, use_gds=FALSE, verbose=TRUE, ...)
{
    # check
    CheckDots(...)
    features <- features %||% rownames(object)
    features <- intersect(features, rownames(object))
    object <- object[features, , drop=FALSE]
    object.names <- dimnames(object)
    # min.cells.to.block <- min(min.cells.to.block, ncol(object))
    # suppressWarnings(expr = Parenting(
    #     parent.find = "ScaleData.Assay",
    #     features = features,
    #     min.cells.to.block = min.cells.to.block
    # ))
    split.by <- split.by %||% TRUE
    split.cells <- split(colnames(object), split.by)
    CheckGC()

    if (is.logical(use_gds))
    {
        if (isTRUE(use_gds))
        {
            use_gds <- "_scale_data.gds"
            while (file.exists(use_gds)) use_gds <- paste0("_", use_gds)
        }
    } else if (!is.character(use_gds))
        stop("'use_gds' should be FALSE, TRUE or a file name.")

    if (!is.null(vars.to.regress))
    {
        stop("Not implemented yet.")
        if (is.null(latent.data))
        {
            latent.data <- data.frame(row.names=colnames(object))
        } else {
            latent.data <- latent.data[colnames(object), , drop = FALSE]
            rownames(latent.data) <- colnames(object)
        }
        if (any(vars.to.regress %in% rownames(object)))
        {
            x <- object[intersect(vars.to.regress, rownames(object)), , drop=FALSE]
            latent.data <- cbind(latent.data, realize(t(x), BACKEND=NULL))
            remove(x)
        }
        notfound <- setdiff(vars.to.regress, colnames(latent.data))
        if (length(notfound) == length(vars.to.regress))
        {
            stop("None of the requested variables to regress are present in the object.",
                call.=FALSE)
        } else if (length(notfound) > 0L)
        {
            warning("Requested variables to regress not in object: ",
                paste(notfound, collapse = ", "),
                call.=FALSE, immediate.=TRUE)
            vars.to.regress <- colnames(latent.data)
        }
        if (verbose)
            message("Regressing out ", paste(vars.to.regress, collapse=', '))

        object <- lapply(names(split.cells), FUN=function(x)
        {
            if (verbose && length(split.cells) > 1L)
                message("Regressing out variables from split ", x)
            RegressOutMatrix(
                data.expr = object[, split.cells[[x]], drop = FALSE],
                latent.data = latent.data[split.cells[[x]], , drop = FALSE],
                features.regress = features,
                model.use = model.use,
                use.umi = use.umi,
                verbose = verbose
            )
        })
        object <- do.call(cbind, object)
        dimnames(object) <- object.names
        CheckGC()
    }

    if (verbose && (do.scale || do.center))
    {
        msg <- paste(na.omit(c(
            ifelse(do.center, 'centering', NA_character_),
            ifelse(do.scale, 'scaling', NA_character_))), collapse = ' and ')
        msg <- paste0(
            toupper(substr(msg, 1L, 1L)), substring(msg, 2L),
            ' data matrix (', class(object)[1L],
            ' [', paste(dim(object), collapse=','), '])')
        message(msg)
    }

    # output variable
    if (is.character(use_gds))
    {
        # use DelayedMatrix
        if (file.exists(use_gds))
        {
            if (verbose) message("Open ", sQuote(use_gds))
            outf <- openfn.gds(use_gds, readonly=FALSE)
        } else {
            if (verbose) message("Create ", sQuote(use_gds))
            outf <- createfn.gds(use_gds)
        }
        on.exit(closefn.gds(outf))
        out_nd <- add.gdsn(outf, "scale.data", storage="double",
            valdim=dim(object), replace=TRUE)
    } else {
        scaled.data <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=object.names)
    }
    # scale
    for (x in names(split.cells))
    {
        if (verbose)
        {
            if (length(split.cells)>1 && (do.scale || do.center))
                message(gsub('matrix', 'from split ', msg), x)
        }
        ii <- match(split.cells[[x]], colnames(object))
        m <- object[, ii, drop=FALSE]
        # center & scale, m is DelayedMatrix
        m <- x_row_scale(m, do.center, do.scale)
        # save
        if (!is.character(use_gds))
        {
            for (k in seq_along(ii))
            {
                v <- m[, k, drop=TRUE]
                v[v > scale.max] <- scale.max  # set a bound
                v[is.na(v)] <- 0
                scaled.data[, ii[k]] <- v
            }
        } else {
            for (k in seq_along(ii))
            {
                v <- m[, k, drop=TRUE]
                v[v > scale.max] <- scale.max  # set a bound
                v[is.na(v)] <- 0
                write.gdsn(out_nd, v, start=c(1L, ii[k]), count=c(length(v), 1L))
            }
        }
        # CheckGC()
    }
    if (is.character(use_gds))
    {
        on.exit()
        closefn.gds(outf)  # close the file first
        scaled.data <- scArray(use_gds, "scale.data")
        dimnames(scaled.data) <- object.names
        scaled.data <- scObj(scaled.data)
    }

    return(scaled.data)
}


####  Methods -- FindVariableFeatures()  ####

.row_var_std <- function(mat, mu, sd, vmax, verbose)
{
    stopifnot(is(mat, "DelayedMatrix"))
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

        rownames(hvf.info) <- rownames(object)
        colnames(hvf.info) <- paste0('vst.', colnames(hvf.info))
    } else {
        stop("selection.method!=vst, not implemented yet.")
    }
    return(hvf.info)
}


# as.Seurat


