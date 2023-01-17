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

.log_norm <- function(mat, scale.factor=10000, verbose=TRUE)
{
    stopifnot(is(mat, "DelayedArray"))
    s <- scale.factor / colSums(mat)
    m <- log1p(sweep(mat, 2L, s, `*`))
}

NormalizeData.SC_GDSMatrix <- function(object,
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

x_write_gdsn <- function(mat, gdsn, st=1L, scale.max=10, verbose=TRUE)
{
    stopifnot(is(mat, "DelayedMatrix"))
    stopifnot(is(gdsn, "gdsn.class"))
    if (verbose)
        pb <- txtProgressBar(min=st, max=st+ncol(mat), style=3L, file=stderr())
    # block write
    blockReduce(function(bk, i, gdsn, vmax)
    {
        bk[bk > vmax] <- vmax  # set a bound
        bk[is.na(bk)] <- 0
        write.gdsn(gdsn, bk, start=c(1L, i), count=c(-1L, ncol(bk)))
        i <- i + ncol(bk)
        if (verbose) setTxtProgressBar(pb, i)
        i
    }, mat, init=st, grid=colAutoGrid(mat), gdsn=gdsn, vmax=scale.max)
    # finally
    if (verbose) close(pb)
    invisible()
}

ScaleData.SC_GDSMatrix <- function(object, features=NULL, vars.to.regress=NULL,
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
            i <- 1L
            while (file.exists(use_gds))
            {
                i <- i + 1L
                use_gds <- sprintf("_scale_data%d.gds", i)
            }
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
            latent.data <- cbind(latent.data, as.matrix(t(x)))
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
        if (verbose) message("Writing to ", sQuote(use_gds))
        if (file.exists(use_gds))
        {
            outf <- openfn.gds(use_gds, readonly=FALSE)
        } else {
            outf <- createfn.gds(use_gds)
        }
        on.exit(closefn.gds(outf))
        out_nd <- add.gdsn(outf, "scale.data", storage="double",
            valdim=dim(object), replace=TRUE)
        gds_colnm <- NULL
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
        if (is.character(use_gds))
        {
            x_write_gdsn(m, out_nd, length(gds_colnm)+1L, scale.max, verbose)
            gds_colnm <- c(gds_colnm, split.cells[[x]])
        } else {
            if (verbose)
                pb <- txtProgressBar(min=0, max=length(ii), style=3, file=stderr())
            for (k in seq_along(ii))
            {
                v <- m[, k, drop=TRUE]
                v[v > scale.max] <- scale.max  # set a bound
                v[is.na(v)] <- 0
                scaled.data[, ii[k]] <- v
                if (verbose && (k %% 1000L==1L))
                    setTxtProgressBar(pb, k)
            }
            if (verbose)
            {
                setTxtProgressBar(pb, length(ii))
                close(pb)
            }
            CheckGC()
        }
    }
    if (is.character(use_gds))
    {
        on.exit()
        closefn.gds(outf)  # close the file first
        scaled.data <- scArray(use_gds, "scale.data")
        if (any(gds_colnm != colnames(object)))
        {
            i <- match(colnames(object), gds_colnm)
            scaled.data <- scaled.data[, i]
        }
        dimnames(scaled.data) <- object.names
        scaled.data <- scObj(scaled.data)
    }

    return(scaled.data)
}


####  Methods -- FindVariableFeatures()  ####

.row_var_std <- function(mat, mu, sd, vmax, verbose)
{
    stopifnot(is(mat, "SC_GDSMatrix"))
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

FindVariableFeatures.SC_GDSMatrix <- function(object,
    selection.method="vst", loess.span=0.3, clip.max="auto",
    mean.function=NULL, dispersion.function=NULL, num.bin=20,
    binning.method="equal_width", verbose=TRUE, ...)
{
    # check
    CheckDots(...)

    if (is.null(mean.function))
        mean.function <- Seurat:::FastExpMean
    if (is.null(dispersion.function))
        dispersion.function <- Seurat:::FastLogVMR

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


####  Methods -- RunPCA()  ####

RunPCA.SCArrayAssay <- function(object, assay=NULL, features=NULL, npcs=50,
    rev.pca=FALSE, weight.by.var=TRUE, verbose=TRUE, ndims.print=1:5,
    nfeatures.print=30, reduction.key="PC_", seed.use=42, ...)
{
    # check
    if (length(VariableFeatures(object))==0L && is.null(features))
    {
        stop("Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.")
    }
    data.use <- GetAssayData(object, "scale.data")
    if (NROW(data.use) == 0L)
        stop("Data has not been scaled. Please run ScaleData and retry.")

    # filter (need var > 0)
    features <- features %||% VariableFeatures(object)
    features.keep <- unique(features[features %in% rownames(data.use)])
    if (length(features.keep) < length(features))
    {
        features.exclude <- setdiff(features, features.keep)
        if (verbose)
        {
            warning(paste0("The following ", length(features.exclude),
                " features requested have not been scaled (running reduction without them): ",
                paste0(features.exclude, collapse = ", ")), immediate.=TRUE)
        }
    }
    features <- features.keep
    features.var <- rowVars(data.use[features, ])
    features.keep <- features[features.var > 0]
    if (length(features.keep) < length(features))
    {
        features.exclude <- setdiff(features, features.keep)
        if (verbose)
        {
            warning(paste0("The following ", length(features.exclude),
                " features requested have zero variance (running reduction without them): ",
                paste0(features.exclude, collapse = ", ")), immediate.=TRUE)
        }
    }
    features <- features.keep
    features <- features[!is.na(x = features)]
    data.use <- data.use[features, ]

    # run
    RunPCA(object=data.use,
        assay = assay,
        npcs = npcs,
        rev.pca = rev.pca,
        weight.by.var = weight.by.var,
        verbose = verbose,
        ndims.print = ndims.print,
        nfeatures.print = nfeatures.print,
        reduction.key = reduction.key,
        seed.use = seed.use,
        ...
    )
}

RunPCA.SC_GDSMatrix <- function(object, assay=NULL, npcs=50, rev.pca=FALSE,
    weight.by.var=TRUE, verbose=TRUE, ndims.print=1:5, nfeatures.print=30,
    reduction.key="PC_", seed.use=42, approx=TRUE, ...)
{
    x_check(object, "Calling RunPCA.SC_GDSMatrix() with %s ...")

    if (!is.null(seed.use)) set.seed(seed.use)
    pca_func <- if (isTRUE(approx)) runIrlbaSVD else runExactSVD
    if (rev.pca)
    {
        total.variance <- sum(colVars(object))
        npcs <- min(npcs, ncol(object)-1L)
        pca_rv <- pca_func(object, k=npcs, center=FALSE, scale=FALSE,
            deferred=FALSE, fold=1)
        sdev <- pca_rv$d / sqrt(max(1L, nrow(object)-1L))
        if (weight.by.var)
        {
            feature.loadings <- pca_rv$u %*% diag(pca_rv$d)
        } else {
            feature.loadings <- pca_rv$u
        }
        cell.embeddings <- pca_rv$v
    } else {
        total.variance <- sum(rowVars(object))
        npcs <- min(npcs, nrow(object)-1L)
        pca_rv <- pca_func(t(object), k=npcs, center=FALSE, scale=FALSE,
            deferred=FALSE, fold=1)
        feature.loadings <- pca_rv$v
        sdev <- pca_rv$d / sqrt(max(1L, ncol(object)-1L))
        if (weight.by.var)
        {
            cell.embeddings <- pca_rv$u %*% diag(pca_rv$d)
        } else {
            cell.embeddings <- pca_rv$u
        }
    }

    rownames(feature.loadings) <- rownames(object)
    colnames(feature.loadings) <- paste0(reduction.key, 1:npcs)
    rownames(cell.embeddings) <- colnames(object)
    colnames(cell.embeddings) <- colnames(feature.loadings)
    reduction.data <- CreateDimReducObject(
        embeddings = cell.embeddings,
        loadings = feature.loadings,
        assay = assay,
        stdev = sdev,
        key = reduction.key,
        misc = list(total.variance=total.variance)
    )
    if (verbose)
    {
        msg <- capture.output(print(
            x = reduction.data,
            dims = ndims.print,
            nfeatures = nfeatures.print
        ))
        message(paste(msg, collapse='\n'))
    }
    return(reduction.data)
}








# as.Seurat


