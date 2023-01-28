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


#######################################################################
# Internal functions

.find_filename <- function(fn)
{
    stopifnot(is.character(fn), length(fn)==1L, !is.na(fn))
    stopifnot(grepl("%d", fn, fixed=TRUE))
    i <- 1L
    s <- gsub("%d", "", fn, fixed=TRUE)
    while (file.exists(s))
    {
        i <- i + 1L
        s <- sprintf(fn, i)
    }
    s
}

# ensure row- and column-names are vectors, not arrays
.form_mat_names <- function(mat)
{
    if (!is.vector(rownames(mat)))
        rownames(mat) <- as.vector(rownames(mat))
    if (!is.vector(colnames(mat)))
        colnames(mat) <- as.vector(colnames(mat))
    mat
}

# check row- and column- names
.check_mat <- function(m)
{
    # check that dimnames of input counts are unique
    if (anyDuplicated(rownames(m)))
    {
        warning(
            "Non-unique features (rownames) present in the input matrix, making unique",
            call.=FALSE, immediate.=TRUE)
        rownames(m) <- make.unique(rownames(m))
    }
    if (anyDuplicated(colnames(m)))
    {
        warning(
            "Non-unique cell names (colnames) present in the input matrix, making unique",
            call.=FALSE, immediate.=TRUE)
        colnames(m) <- make.unique(colnames(m))
    }
    if (is.null(colnames(m)))
    {
        stop("No cell names (colnames) names present in the input matrix",
            call.=FALSE)
    }
    if (any(rownames(m) == ''))
    {
        stop("Feature names of counts matrix cannot be empty",
            call.=FALSE)
    }
    if (nrow(m) > 0L && is.null(rownames(m)))
    {
        stop("No feature names (rownames) names present in the input matrix",
            call.=FALSE)
    }
    # output
    m
}



#######################################################################

# Create an SCArrayAssay object from counts or data
CreateAssayObject2 <- function(counts, data, min.cells=0, min.features=0,
    check.matrix=FALSE, ...)
{
    if (missing(counts) && missing(data))
    {
        stop("Must provide either 'counts' or 'data'")
    } else if (!missing(counts) && !missing(data))
    {
        stop("Either 'counts' or 'data' must be missing; both cannot be provided")
    } else if (!missing(counts))
    {
        # if not DelayedArray
        # should call SeuratObject::CreateAssayObject() instead
        if (!is(counts, "DelayedArray"))
        {
            return(CreateAssayObject(counts,
                min.cells = min.cells, min.features = min.features,
                check.matrix = check.matrix, ...))
        }
        # check counts
        x_check(counts, "Calling CreateAssayObject2() with %s ...")
        counts <- .check_mat(counts)
        if (isTRUE(check.matrix)) CheckMatrix(counts)
        # filter based on min.features
        if (min.features > 0)
        {
            nfeatures <- colSums(counts > 0)
            counts <- counts[, which(nfeatures >= min.features)]
        }
        # filter genes on the number of cells expressing
        if (min.cells > 0)
        {
            num.cells <- rowSums(counts > 0)
            counts <- counts[which(num.cells >= min.cells), ]
        }
        # set data
        data <- counts <- scObj(counts)

    } else if (!missing(data))
    {
        # if not DelayedArray
        # should call SeuratObject::CreateAssayObject() instead
        if (!is(data, "DelayedArray"))
        {
            return(CreateAssayObject(data=data,
                min.cells = min.cells, min.features = min.features,
                check.matrix = check.matrix, ...))
        }
        # check data
        x_check(data, "Calling CreateAssayObject2() with %s ...")
        data <- .check_mat(data)
        if (min.cells != 0 || min.features != 0)
        {
            warning(
                "No filtering performed if passing to data rather than counts",
                call.=FALSE, immediate.=TRUE)
        }
        # set data
        data <- scObj(data)
        counts <- scObj(DelayedArray(new('matrix')))
    }

    # Ensure row- and column-names
    counts <- .form_mat_names(counts)
    data <- .form_mat_names(data)
    if (any(grepl('_', rownames(counts))) || any(grepl('_', rownames(data))))
    {
        warning(
            "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
            call.=FALSE, immediate.=TRUE)
        rownames(counts) <- gsub('_', '-', rownames(counts))
        rownames(data) <- gsub('_', '-', rownames(data))
    }
    if (any(grepl('|', rownames(counts), fixed=TRUE)) || any(grepl('|', rownames(data), fixed=TRUE)))
    {
        warning(
            "Feature names cannot have pipe characters ('|'), replacing with dashes ('-')",
            call.=FALSE, immediate.=TRUE)
        rownames(counts) <- gsub('|', '-', rownames(counts), fixed=TRUE)
        rownames(data) <- gsub('|', '-', rownames(data), fixed=TRUE)
    }
  
    # output
    new(Class = "SCArrayAssay",
        counts2 = counts, data2 = data, scale.data2 = NULL,
        meta.features = data.frame(row.names = rownames(data)),
        misc = list())
}



#######################################################################
# S3/S4 Methods for SCArrayAssay

# S3 method for CreateSeuratObject()
# Create a Seurat Object from a DelayedMatrix
CreateSeuratObject.DelayedMatrix <- function(counts, project='SeuratProject',
    assay='RNA', names.field=1, names.delim='_', meta.data=NULL, min.cells=0,
    min.features=0, row.names=NULL, ...)
{
    # check
    x_check(counts, "Calling CreateSeuratObject.DelayedMatrix() with %s ...")
    if (!is.null(meta.data))
    {
        if (!all(rownames(meta.data) %in% colnames(counts)))
        {
            warning("Some cells in meta.data not present in provided counts matrix",
                immediate.=TRUE)
        }
    }
    # new SCArrayAssay
    assay.data <- CreateAssayObject2(counts,
        min.cells=min.cells, min.features=min.features, row.names=row.names)
    # meta data & key
    if (!is.null(meta.data))
    {
        common.cells <- intersect(rownames(meta.data), colnames(assay.data))
        meta.data <- meta.data[common.cells, , drop = FALSE]
    }
    Key(assay.data) <- suppressWarnings(Seurat:::UpdateKey(tolower(assay)))
    # output
    CreateSeuratObject(assay.data, project, assay,
        names.field=names.field, names.delim=names.delim,
        meta.data=meta.data, ...)
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

x_row_scale <- function(x, center=TRUE, scale=TRUE, scale.max=10)
{
    if (center && scale)
    {
        v <- scRowMeanVar(x, na.rm=TRUE)
        m <- v[,1L]
        rsd <- 1 / sqrt(v[,2L])
        rsd[!is.finite(rsd)] <- 0
        x <- (x - m) * rsd
    } else {
        if (center)
        {
            m <- rowMeans(x)
            x <- x - m
        } else if (scale)
        {
            rsd <- 1 / rowSds(x, center=rep(0, nrow(x)))
            rsd[!is.finite(rsd)] <- 0
            x <- x * rsd
        }
    }
    scSetMax(x, scale.max)    # set a bound
}

x_append_gdsn <- function(mat, gdsn, verbose=TRUE)
{
    stopifnot(is(mat, "DelayedMatrix"))
    stopifnot(is(gdsn, "gdsn.class"))
    if (verbose)
        pb <- txtProgressBar(min=0, max=ncol(mat), style=3L, file=stderr())
    # block write
    blockReduce(function(bk, i, gdsn)
    {
        append.gdsn(gdsn, bk)
        if (verbose) setTxtProgressBar(pb, i+ncol(bk))
        i + ncol(bk)
    }, mat, init=0L, grid=colAutoGrid(mat), gdsn=gdsn)
    # finally
    if (verbose) close(pb)
    invisible()
}

ScaleData.SC_GDSMatrix <- function(object, features=NULL, vars.to.regress=NULL,
    latent.data=NULL, split.by=NULL, model.use='linear', use.umi=FALSE,
    do.scale=TRUE, do.center=TRUE, scale.max=10, use_gds=TRUE, verbose=TRUE,
    block.size=1000, min.cells.to.block=3000, ...)
{
    # check
    CheckDots(...)
    if (is.null(features)) features <- rownames(object)
    features <- intersect(features, rownames(object))
    object <- object[features, , drop=FALSE]
    object.names <- dimnames(object)
    # min.cells.to.block <- min(min.cells.to.block, ncol(object))
    # suppressWarnings(expr = Parenting(
    #     parent.find = "ScaleData.Assay",
    #     features = features,
    #     min.cells.to.block = min.cells.to.block
    # ))
    if (is.null(split.by)) split.by <- TRUE
    split.cells <- split(colnames(object), split.by)
    CheckGC()

    if (is.logical(use_gds))
    {
        if (isTRUE(use_gds))
            use_gds <- .find_filename("_scale_data%d.gds")
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
    outf <- NULL
    if (is.character(use_gds) && length(split.cells)>1L)
    {
        # use DelayedMatrix
        if (verbose) message("Writing to ", sQuote(use_gds))
        if (file.exists(use_gds))
        {
            outf <- openfn.gds(use_gds, readonly=FALSE)
        } else {
            outf <- createfn.gds(use_gds)
        }
        on.exit(closefn.gds(outf))  # in case fail
        out_nd <- add.gdsn(outf, "scale.data", storage="double",
            valdim=c(nrow(object), 0L), replace=TRUE)
    }

    # Scale each cell split
    lst <- lapply(names(split.cells), function(x)
    {
        m <- object
        # subsetting
        if (!identical(split.cells[[x]], colnames(m)))
        {
            ii <- match(split.cells[[x]], colnames(object))
            m <- object[, ii, drop=FALSE]
        }
        if (verbose)
        {
            if (length(split.cells)>1L && (do.scale || do.center))
            {
                message("Data split (", class(m)[1L], " [",
                    paste(dim(m), collapse=","), "]): ", x)
            }
        }
        # center & scale, m is DelayedMatrix
        m <- x_row_scale(m, do.center, do.scale, scale.max)
        # save
        if (!is.null(outf))
        {
            x_append_gdsn(m, out_nd, verbose)
            split.cells[[x]]  # output
        } else if (is.character(use_gds))
        {
            m  # output
        } else {
            # in-memory
            as.matrix(m)  # output
        }
    })

    # output
    if (!is.null(outf))
    {
        on.exit()
        closefn.gds(outf)  # close the file first
        scaled.data <- scArray(use_gds, "scale.data")
        gds_colnm <- unlist(lst)
        if (!identical(gds_colnm, colnames(object)))
        {
            i <- match(colnames(object), gds_colnm)
            scaled.data <- scaled.data[, i]
        }
        dimnames(scaled.data) <- object.names
    } else if (is.character(use_gds))
    {
        # length(split.cells) == 1, and a single DelayedMatrix
        scaled.data <- lst[[1L]]
    } else {
        # in-memory
        scaled.data <- do.call(cbind, lst)
        gds_colnm <- colnames(scaled.data)
        if (!identical(gds_colnm, colnames(object)))
        {
            i <- match(colnames(object), gds_colnm)
            scaled.data <- scaled.data[, i]
        }
    }
    return(scaled.data)
}


####  Methods -- FindVariableFeatures()  ####

.row_var_std <- function(x, mu, sd, vmax, verbose)
{
    # check
    stopifnot(is(x, "SC_GDSMatrix"))
    # initialize
    inv <- 1 / sd
    inv[!is.finite(inv)] <- 0
    if (verbose)
        pb <- txtProgressBar(min=0L, max=ncol(x), style=3L, file=stderr())
    # block read
    v <- blockReduce(function(bk, v, mu, inv, vmax, vb)
    {
        b <- pmin((bk - mu)*inv, vmax)^2L
        if (vb)
            setTxtProgressBar(pb, start(currentViewport())[2L])
        v + rowSums(b)
    }, x, 0, grid=colAutoGrid(x), mu=mu, inv=inv, vmax=vmax, vb=verbose)
    # finally
    if (verbose)
    {
        setTxtProgressBar(pb, ncol(x))
        close(pb)
    }
    # output
    v[!is.finite(v)] <- 0
    v / (ncol(x) - 1L)
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
        v <- scRowMeanVar(object)
        hvf.info <- data.frame(mean=v[,1L], variance=v[,2L])
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
        stop("selection.method != 'vst', not implemented yet.")
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
    if (is.null(features)) features <- VariableFeatures(object)
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

    # BiocSingular SVD functions
    pca_func <- if (isTRUE(approx)) runIrlbaSVD else runExactSVD

    if (!is.null(seed.use)) set.seed(seed.use)
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
            x=reduction.data, dims=ndims.print, nfeatures=nfeatures.print
        ))
        message(paste(msg, collapse='\n'))
    }
    return(reduction.data)
}








# as.Seurat


