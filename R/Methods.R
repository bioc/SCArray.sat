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
    # switch(nm, data="data2", counts="counts2", scale.data="scale.data2", nm)
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

ScaleData.DelayedMatrix <- function(object, features=NULL, vars.to.regress=NULL,
    latent.data=NULL, split.by=NULL, model.use='linear', use.umi=FALSE,
    do.scale=TRUE, do.center=TRUE, scale.max=10, block.size=1000,
    min.cells.to.block=3000, verbose=TRUE, ...)
{
    CheckDots(...)
    features <- features %||% rownames(object)
    features <- intersect(features, rownames(object))
    object <- object[features, , drop=FALSE]
    object.names <- dimnames(object)
    min.cells.to.block <- min(min.cells.to.block, ncol(object))
    # suppressWarnings(expr = Parenting(
    #     parent.find = "ScaleData.Assay",
    #     features = features,
    #     min.cells.to.block = min.cells.to.block
    # ))
    split.by <- split.by %||% TRUE
    split.cells <- split(colnames(object), split.by)
    CheckGC()

    if (!is.null(vars.to.regress))
    {
        if (is.null(x = latent.data))
        {
            latent.data <- data.frame(row.names = colnames(x = object))
        } else {
            latent.data <- latent.data[colnames(x = object), , drop = FALSE]
            rownames(x = latent.data) <- colnames(x = object)
        }
        if (any(vars.to.regress %in% rownames(x = object)))
        {
            latent.data <- cbind(latent.data,
                t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], , drop=FALSE])
            )
        }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
    notfound <- setdiff(x = vars.to.regress, y = colnames(x = latent.data))
    if (length(x = notfound) == length(x = vars.to.regress)) {
      stop(
        "None of the requested variables to regress are present in the object.",
        call. = FALSE
      )
    } else if (length(x = notfound) > 0) {
      warning(
        "Requested variables to regress not in object: ",
        paste(notfound, collapse = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
      vars.to.regress <- colnames(x = latent.data)
    }
    if (verbose) {
      message("Regressing out ", paste(vars.to.regress, collapse = ', '))
    }
    chunk.points <- ChunkPoints(dsize = nrow(x = object), csize = block.size)
    if (nbrOfWorkers() > 1) { # TODO: lapply
      chunks <- expand.grid(
        names(x = split.cells),
        1:ncol(x = chunk.points),
        stringsAsFactors = FALSE
      )
      object <- future_lapply(
        X = 1:nrow(x = chunks),
        FUN = function(i) {
          row <- chunks[i, ]
          group <- row[[1]]
          index <- as.numeric(x = row[[2]])
          return(RegressOutMatrix(
            data.expr = object[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
            latent.data = latent.data[split.cells[[group]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = FALSE
          ))
        }
      )
      if (length(x = split.cells) > 1) {
        merge.indices <- lapply(
          X = 1:length(x = split.cells),
          FUN = seq.int,
          to = length(x = object),
          by = length(x = split.cells)
        )
        object <- lapply(
          X = merge.indices,
          FUN = function(x) {
            return(do.call(what = 'rbind', args = object[x]))
          }
        )
        object <- do.call(what = 'cbind', args = object)
      } else {
        object <- do.call(what = 'rbind', args = object)
      }
    } else {
      object <- lapply(
        X = names(x = split.cells),
        FUN = function(x) {
          if (verbose && length(x = split.cells) > 1) {
            message("Regressing out variables from split ", x)
          }
          return(RegressOutMatrix(
            data.expr = object[, split.cells[[x]], drop = FALSE],
            latent.data = latent.data[split.cells[[x]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = verbose
          ))
        }
      )
      object <- do.call(what = 'cbind', args = object)
    }
    dimnames(x = object) <- object.names
    CheckGC()
  }
    if (verbose && (do.scale || do.center))
    {
        msg <- paste(
      na.omit(object = c(
        ifelse(test = do.center, yes = 'centering', no = NA_character_),
        ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
      )),
      collapse = ' and '
        )
        msg <- paste0(
            toupper(x = substr(x = msg, start = 1, stop = 1)),
      substr(x = msg, start = 2, stop = nchar(x = msg)),
      ' data matrix (', class(object)[1L],
      ' [', paste(dim(object), collapse=','), '])'
      )
        message(msg)
    }

    scale_fc_sp <- is_sparse(object)
    if (scale_fc_sp)
    {
        scale.function <- Seurat:::FastSparseRowScale
    } else {
        scale.function <- Seurat:::FastRowScale
    }

    scaled.data <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
        dimnames=object.names)
    max.block <- ceiling(length(features) / block.size)
    for (x in names(split.cells))
    {
        if (verbose)
        {
            if (length(split.cells)>1 && (do.scale || do.center))
                message(gsub('matrix', 'from split ', msg), x)
            pb <- txtProgressBar(0, max.block, style=3L, file=stderr())
        }
        for (i in 1:max.block)
        {
            my.inds <- ((block.size*(i-1L)):(block.size*i - 1L)) + 1L
            my.inds <- my.inds[my.inds <= length(features)]
            m <- object[features[my.inds], split.cells[[x]], drop=FALSE]
            arg.list <- list(
                mat = if (scale_fc_sp) as(m, "sparseMatrix") else as.matrix(m),
                scale = do.scale,
                center = do.center,
                scale_max = scale.max,
                display_progress = FALSE
            )
            arg.list <- arg.list[intersect(names(arg.list), names(formals(scale.function)))]
            # call
            data.scale <- do.call(scale.function, arg.list)
            dimnames(data.scale) <- dimnames(object[features[my.inds], split.cells[[x]]])
            scaled.data[features[my.inds], split.cells[[x]]] <- data.scale
            rm(data.scale, m, arg.list)
            CheckGC()
            if (verbose) setTxtProgressBar(pb, i)
        }
        if (verbose) close(pb)
    }

    dimnames(scaled.data) <- object.names
    scaled.data[is.na(scaled.data)] <- 0
    CheckGC()
    return(scaled.data)
}


####  Methods -- FindVariableFeatures()  ####

.row_var_std <- function(mat, mu, sd, vmax, verbose)
{
    stopifnot(is(mat, "DelayedArray"))
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






