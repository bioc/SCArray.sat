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
        warning("Non-unique features (rownames) present in the input matrix",
            ", making unique", call.=FALSE, immediate.=TRUE)
        rownames(m) <- make.unique(rownames(m))
    }
    if (anyDuplicated(colnames(m)))
    {
        warning("Non-unique cell names (colnames) present in the input matrix",
            ", making unique", call.=FALSE, immediate.=TRUE)
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
        stop("Must provide either 'counts' or 'data'.")
    } else if (!missing(counts) && !missing(data))
    {
        stop("Either 'counts' or 'data' must be missing; ",
            "both cannot be provided.")
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
        data <- counts <- scObj(counts)  # wrap SC_GDSMatrix if needed

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
        data <- scObj(data)  # wrap SC_GDSMatrix if needed
        counts <- scObj(DelayedArray(new('matrix')))
    }

    # Ensure row- and column-names
    counts <- .form_mat_names(counts)
    data <- .form_mat_names(data)
    if (any(grepl('_', c(rownames(counts), rownames(data)), fixed=TRUE)))
    {
        warning("Feature names cannot have underscores ('_'), ",
            "replacing with dashes ('-')", call.=FALSE, immediate.=TRUE)
        rownames(counts) <- gsub('_', '-', rownames(counts), fixed=TRUE)
        rownames(data) <- gsub('_', '-', rownames(data), fixed=TRUE)
    }
    if (any(grepl('|', c(rownames(counts), rownames(data)), fixed=TRUE)))
    {
        warning("Feature names cannot have pipe characters ('|'), ",
            "replacing with dashes ('-')", call.=FALSE, immediate.=TRUE)
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
    assay='rna_', names.field=1, names.delim='_', meta.data=NULL, min.cells=0,
    min.features=0, row.names=NULL, ...)
{
    # check
    x_check(counts, "Calling CreateSeuratObject.DelayedMatrix() with %s ...")
    if (!is.null(meta.data))
    {
        if (!all(rownames(meta.data) %in% colnames(counts)))
        {
            warning(
                "Some cells in meta.data not present in provided counts matrix",
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
    Key(assay.data) <- Seurat:::UpdateKey(tolower(assay))
    # output
    CreateSeuratObject(assay.data, project, assay,
        names.field=names.field, names.delim=names.delim,
        meta.data=meta.data, ...)
}



####  Methods -- NormalizeData()  ####

.log_norm <- function(x, scale.factor=1e4, verbose=TRUE)
{
    stopifnot(is(x, "DelayedArray"))
    if (verbose)
        .cat("Performing log-normalization")
    s <- scale.factor / colSums(x)
    log1p(sweep(x, 2L, s, `*`))
}

.clr_norm <- function(x, margin, verbose)
{
    if (margin == 1L)
    {
        if (verbose)
            .cat("Normalizing across features (CLR)")
        s <- exp(-rowSums(log1p(x)) / ncol(x))
        log1p(x * s)
    } else {
        if (verbose)
            .cat("Normalizing across cells (CLR)")
        s <- exp(-colSums(log1p(x)) / nrow(x))
        log1p(sweep(x, 2L, s, `*`))
    }
}

.rc_norm <- function(x, scale.factor=1, verbose=TRUE)
{
    stopifnot(is(x, "DelayedArray"))
    if (verbose)
        .cat("Performing relative-counts-normalization")
    s <- scale.factor / colSums(x)
    sweep(x, 2L, s, `*`)
}


NormalizeData.SC_GDSMatrix <- function(object,
    normalization.method="LogNormalize", scale.factor=1e4, margin=1,
    verbose=TRUE, ...)
{
    # check
    x_msg("Calling NormalizeData.SC_GDSMatrix() ...")
    CheckDots(...)
    if (!is.null(normalization.method))
    {
        stopifnot(is.character(normalization.method),
            length(normalization.method)==1L)
    }
    stopifnot(margin %in% c(1, 2))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (is.null(normalization.method)) return(object)
    # output
    switch(normalization.method,
        "LogNormalize" = .log_norm(object, scale.factor, verbose),
        "CLR" = .clr_norm(object, margin, verbose),
        "RC"  = .rc_norm(object, scale.factor, verbose),
        stop("Unknown or not implemented normalization method: ",
            normalization.method)
    )
}


####  Methods -- ScaleData()  ####

x_row_scale <- function(x, center=TRUE, scale=TRUE, scale.max=10)
{
    stopifnot(is.numeric(scale.max), length(scale.max)==1L)
    if (center && scale)
    {
        v <- scRowMeanVar(x, na.rm=TRUE)  # calculate mean and var together
        m <- v[,1L]              # mean
        rsd <- 1 / sqrt(v[,2L])  # var
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
    # set a bound and output
    if (is(x, "SC_GDSArray"))
    {
        scSetMax(x, scale.max)
    } else {
        x[x > scale.max] <- scale.max
        x
    }
}

x_regress_out <- function(x, latent.data=NULL, out_nd=NULL, 
    model.use=c("linear", "poisson", "negbinom"), verbose=TRUE)
{
    # check
    model.use <- match.arg(model.use)
    if (is.null(latent.data)) return(x)
    if (nrow(latent.data) != ncol(x))
        stop("Uneven number of cells between latent data and expression data.")
    # initialize
    vars.to.regress <- colnames(latent.data)
    fm <- as.formula(paste("GENE ~", paste(vars.to.regress, collapse="+")))
    if (model.use == "linear")
    {
        mat <- cbind(latent.data, x[1L, ])
        colnames(mat) <- c(colnames(latent.data), "GENE")
        qr <- lm(fm, mat, qr=TRUE)$qr
    }

    gd <- rowAutoGrid(x)
    if (verbose)
        pb <- txtProgressBar(0L, length(gd), style=3L, width=64L, file=stderr())
    lst <- blockApply(x, function(bk)
    {
        if (is(bk, "SparseArraySeed"))
            bk <- as(bk, "sparseMatrix")
        ans <- NULL
        if (is.null(out_nd))
            ans <- matrix(0, nrow=nrow(bk), ncol=ncol(bk))
        for (i in seq_len(nrow(bk)))
        {
            xx <- bk[i, ]
            mat <- cbind(latent.data, xx)
            colnames(mat) <- c(vars.to.regress, "GENE")
            res <- switch(model.use,
                linear = qr.resid(qr, xx),
                poisson = residuals(glm(fm, "poisson", mat), type="pearson"),
                negbinom = Seurat:::NBResiduals(fm, mat, rownames(x)[i])
            )
            # set
            if (is.null(out_nd))
            {
                ans[i, ] <- res
            } else {
                append.gdsn(out_nd, res)
            }
        }
        if (verbose)
            setTxtProgressBar(pb, currentBlockId())
        ans
    }, grid=gd, as.sparse=NA, BPPARAM=NULL)
    if (verbose) close(pb)

    # output
    if (is.null(out_nd))
    {
        ans <- do.call(rbind, lst)
        colnames(ans) <- colnames(x)
        ans
    } else {
        colnames(x)
    }
}

ScaleRegressOut <- function(object, features, vars.to.regress, latent.data,
    split.cells, use_gds, model.use='linear', use.umi=FALSE, verbose=TRUE)
{
    # check
    if (verbose)
        .cat("Regressing out: ", paste(vars.to.regress, collapse=", "))
    if (is.null(latent.data))
    {
        latent.data <- data.frame(row.names=colnames(object))
    } else {
        latent.data <- latent.data[colnames(object), , drop = FALSE]
        rownames(latent.data) <- colnames(object)
    }
    if (any(vars.to.regress %in% rownames(object)))
    {
        s <- intersect(vars.to.regress, rownames(object))
        if (verbose)
            .cat("Loading the matrix for ", paste(s, collapse=", "))
        x <- object[s, , drop=FALSE]
        latent.data <- cbind(latent.data, as.matrix(t(x)))
        remove(x)
    }
    notfound <- setdiff(vars.to.regress, colnames(latent.data))
    if (length(notfound) == length(vars.to.regress))
    {
        stop("None of the requested variables to regress ",
            "are present in the object.", call.=FALSE)
    } else if (length(notfound) > 0L)
    {
        warning("Requested variables to regress not in object: ",
            paste(notfound, collapse = ", "), call.=FALSE, immediate.=TRUE)
        vars.to.regress <- colnames(latent.data)
    }

    # use DelayedMatrix ?
    resid_gdsfn <- outf <- out_nd <- NULL
    if (is.character(use_gds))
    {
        resid_gdsfn <- paste0("_temp", use_gds)
        if (verbose)
            .cat("Writing to ", sQuote(resid_gdsfn))
        if (file.exists(use_gds))
        {
            warning("Overwriting the file '", resid_gdsfn, "'",
                call.=FALSE, immediate.=TRUE)
        }
        # create a GDS file for regression residuals
        outf <- createfn.gds(resid_gdsfn)
        on.exit(closefn.gds(outf))  # in case fail
    }
    if (verbose && x_warn_speed(object))
    {
        message("Regressing maybe faster with a larger block size ",
            "via setAutoBlockSize()")
    }

    # run regression
    lst <- lapply(seq_along(split.cells), FUN=function(i)
    {
        if (verbose && length(split.cells) > 1L)
            .cat("Regressing out variables from split ", names(split.cells)[i])
        ss <- split.cells[[i]]
        if (!is.null(outf))
        {
            # a new GDS node storing matrix
            out_nd <- add.gdsn(outf, sprintf("residuals.%d", i),
                storage="double", valdim=c(length(ss), 0L))
        }
        x_regress_out(
            object[, ss, drop=FALSE], latent.data[ss, , drop=FALSE],
            out_nd, model.use, verbose)
    })

    # finally
    if (is.null(outf))
    {
        # in-memory matrix
        x <- do.call(cbind, lst)
        gds_colnm <- colnames(x)
    } else {
        # on-disk backend
        on.exit()
        if (length(split.cells) == 1L)
        {
            closefn.gds(outf)  # close the GDS file
            x <- t(scArray(resid_gdsfn, "residuals.1"))
        } else {
            # merge
            closefn.gds(outf)  # close the GDS file
            outf <- openfn.gds(resid_gdsfn, readonly=FALSE)
            out_nd <- add.gdsn(outf, "residuals", storage="double",
                valdim=c(nrow(object), 0L))
            for (i in seq_along(split.cells))
            {
                # block write
                m <- scArray(outf, sprintf("residuals.%d", i))
                blockReduce(function(bk, v) append.gdsn(out_nd, t(bk)),
                    m, init=NULL, grid=rowAutoGrid(m))
            }
            remove(m)
            closefn.gds(outf)  # close the GDS file
            x <- scArray(resid_gdsfn, "residuals")
        }
        gds_colnm <- unlist(lst)
    }
    if (!identical(gds_colnm, colnames(object)))
    {
        i <- match(colnames(object), gds_colnm)
        object <- x[, i]
    } else
        object <- x

    use.umi <- ifelse(model.use!="linear", TRUE, use.umi)
    if (use.umi)
    {
        m <- rowMins(object, na.rm=TRUE)
        object <- log1p(object - m)
    }

    CheckGC()
    # output
    list(object=object, filename=resid_gdsfn)
}


ScaleData.SC_GDSMatrix <- function(object, features=NULL, vars.to.regress=NULL,
    latent.data=NULL, split.by=NULL, model.use='linear', use.umi=FALSE,
    do.scale=TRUE, do.center=TRUE, scale.max=10, block.size=1000,
    min.cells.to.block=3000, verbose=TRUE, use_gds=TRUE, rm_tmpfile=TRUE, ...)
{
    # check
    x_check(object, "Calling ScaleData.SC_GDSMatrix() with %s ...")
    CheckDots(...)
    if (!is.null(features))
    {
        features <- intersect(features, rownames(object))
        object <- object[features, , drop=FALSE]
    } else {
        features <- rownames(object)
    }
    if (anyNA(split.by))
        stop("'split.by' should not contain NA.")
    if (is.null(split.by)) split.by <- TRUE
    split.cells <- split(colnames(object), split.by)
    object.names <- dimnames(object)
    CheckGC()

    # use DelayedMatrix ?
    resid_gdsfn <- NULL    # temporary file name
    force_gds <- is.character(use_gds)
    if (is.logical(use_gds))
    {
        if (isTRUE(use_gds))
        {
            use_gds <- attr(use_gds, "gdsfn")
            if (!is.character(use_gds))
                use_gds <- .find_filename("_scale_data%d.gds")
        }
    } else if (!is.character(use_gds) || is.na(use_gds))
    {
        stop("'use_gds' should be FALSE, TRUE or a file name.")
    }

    if (!is.null(vars.to.regress))
    {
        v <- ScaleRegressOut(object, features, vars.to.regress,
            latent.data, split.cells, use_gds, model.use, use.umi, verbose)
        object <- v$object
        dimnames(object) <- object.names
        resid_gdsfn <- v$filename
    }

    if (verbose && (do.scale || do.center))
    {
        msg <- paste(na.omit(c(
            ifelse(do.center, 'centering', NA_character_),
            ifelse(do.scale, 'scaling', NA_character_))), collapse = ' and ')
        msg <- paste0(
            toupper(substr(msg, 1L, 1L)), substring(msg, 2L),
            ' data matrix (', class(object)[1L],
            ' [', paste(dim(object), collapse='x'), '])')
        .cat(msg)
    }

    # output variable
    out_nd <- NULL
    if (is.character(use_gds) &&
        (force_gds || length(split.cells)>1L || !is.null(vars.to.regress)))
    {
        # use DelayedMatrix
        if (verbose)
            .cat("Writing to ", sQuote(use_gds))
        if (is.null(vars.to.regress))
        {
            if (file.exists(use_gds))
            {
                warning("Overwriting the file '", use_gds, "'",
                    call.=FALSE, immediate.=TRUE)
            }
        }
        # create a GDS file
        outf <- createfn.gds(use_gds)
        on.exit(closefn.gds(outf))  # in case fail
        # a new GDS node storing matrix
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
                .cat("Data split (", class(m)[1L], " [",
                    paste(dim(m), collapse=","), "]): ", x)
            }
        }
        # center & scale, return a DelayedMatrix
        m <- x_row_scale(m, do.center, do.scale, scale.max)
        # save
        if (!is.null(out_nd))
        {
            x_append_gdsn(m, out_nd, verbose)
            split.cells[[x]]  # output colnames
        } else if (is.character(use_gds))
        {
            m  # output, no need a new GDS file
        } else {
            # in-memory
            as.matrix(m)  # output
        }
    })

    # output
    if (!is.null(out_nd))
    {
        on.exit()
        closefn.gds(outf)  # close the GDS file
        scaled.data <- scArray(use_gds, "scale.data")
        gds_colnm <- unlist(lst)
    } else if (is.character(use_gds))
    {
        # length(split.cells) == 1, and a single DelayedMatrix
        scaled.data <- lst[[1L]]
        gds_colnm <- colnames(scaled.data)
    } else {
        # in-memory
        scaled.data <- do.call(cbind, lst)
        gds_colnm <- colnames(scaled.data)
    }
    if (!identical(gds_colnm, colnames(object)))
    {
        i <- match(colnames(object), gds_colnm)
        scaled.data <- scaled.data[, i]
    }
    dimnames(scaled.data) <- object.names

    # clean temp data
    if (!is.null(resid_gdsfn) && isTRUE(rm_tmpfile))
    {
        if (verbose) .cat("Delete ", sQuote(resid_gdsfn))
        unlink(resid_gdsfn, force=TRUE)
    }

    CheckGC()
    return(scaled.data)
}


####  Methods -- FindVariableFeatures()  ####

.row_var_std <- function(x, mu, sd, vmax, verbose=TRUE)
{
    # check
    stopifnot(is(x, "SC_GDSMatrix"))
    # initialize
    inv <- 1 / sd
    inv[!is.finite(inv)] <- 0
    gd <- colAutoGrid(x)
    if (verbose)
        pb <- txtProgressBar(0L, length(gd), style=3L, width=64L, file=stderr())
    # block read
    v <- blockReduce(function(bk, v, mu, inv, vmax, vb)
    {
        b <- pmin((bk - mu)*inv, vmax)^2L
        if (vb) setTxtProgressBar(pb, currentBlockId())
        v + rowSums(b)
    }, x, 0, grid=gd, mu=mu, inv=inv, vmax=vmax, vb=verbose)
    # finally
    if (verbose) close(pb)
    # output
    v[!is.finite(v)] <- 0
    v / (ncol(x) - 1L)
}

.mean_disp_exp <- function(x, verbose=TRUE)
{
    # check
    stopifnot(is(x, "SC_GDSMatrix"))
    # initialize
    gd <- colAutoGrid(x)
    if (verbose)
        pb <- txtProgressBar(0L, length(gd), style=3L, width=64L, file=stderr())
    # block read
    v <- blockReduce(function(bk, v, vb)
    {
        bk <- expm1(bk)
        if (vb) setTxtProgressBar(pb, currentBlockId())
        v + c(rowSums(bk), rowSums(bk * bk))
    }, x, double(nrow(x)*2L), grid=gd, vb=verbose)
    # finally
    if (verbose) close(pb)
    # output
    v <- matrix(v, nrow(x), ncol=2L)
    s1 <- v[,1L]; s2 <- v[,2L]
    m <- s1 / ncol(x)
    vr <- (s2 - s1*s1/ncol(x)) / (ncol(x) - 1L)
    cbind(log1p(m), log(vr/m))
}

FindVariableFeatures.SC_GDSMatrix <- function(object,
    selection.method="vst", loess.span=0.3, clip.max="auto",
    mean.function=NULL, dispersion.function=NULL, num.bin=20,
    binning.method="equal_width", verbose=TRUE, ...)
{
    # check
    x_msg("Calling FindVariableFeatures.SC_GDSMatrix() ...")
    stopifnot(is.character(selection.method), length(selection.method)==1L)
    CheckDots(...)

    # check mean & dispersion functions
    if (!is.null(mean.function))
    {
        if (!identical(mean.function, Seurat:::FastExpMean) &&
            !identical(mean.function, Seurat::ExpMean))
            stop("User-defined 'mean.function' is not supported.")
    }
    if (!is.null(dispersion.function))
    {
        if (!identical(dispersion.function, Seurat:::FastLogVMR) &&
            !identical(dispersion.function, Seurat::LogVMR))
            stop("User-defined 'dispersion.function' is not supported.")
    }

    if (selection.method == "vst")
    {
        if (clip.max == "auto")
            clip.max <- sqrt(ncol(object))
        if (verbose)
            .cat("Calculating gene variances")
        v <- scRowMeanVar(object)
        hvf <- data.frame(mean=v[,1L], variance=v[,2L])
        hvf$variance[is.na(hvf$variance)] <- 0
        hvf$variance.expected <- 0
        hvf$variance.standardized <- 0
        not.const <- hvf$variance > 0
        fit <- loess(log10(variance) ~ log10(mean),
            hvf[not.const, ], span=loess.span)
        hvf$variance.expected[not.const] <- 10 ^ fit$fitted

        # get variance after feature standardization
        if (verbose)
        {
            .cat("Calculating feature variances ",
                "of standardized and clipped values")
        }
        hvf$variance.standardized <- .row_var_std(
            object, hvf$mean, sqrt(hvf$variance.expected),
            clip.max, verbose)
        colnames(hvf) <- paste0('vst.', colnames(hvf))

    } else {

        # calculate the mean and dispersion for each feature with
        #     normalized (log-transformed) data
        v <- .mean_disp_exp(object, verbose)
        f_mean <- v[, 1L]
        f_mean[is.na(f_mean)] <- 0
        f_dispersion <- v[, 2L]
        f_dispersion[is.na(f_dispersion)] <- 0
        names(f_mean) <- names(f_dispersion) <- rownames(object)
        data.x.breaks <- switch(binning.method,
            'equal_width' = num.bin,
            'equal_frequency' = quantile(f_mean[f_mean>0],
                probs=seq(0, 1, length.out=num.bin)),
            stop("Unknown binning method: ", binning.method)
        )
        data.x.bin <- cut(f_mean, breaks=data.x.breaks, include.lowest=TRUE)
        names(data.x.bin) <- names(f_mean)
        mean.y <- tapply(f_dispersion, data.x.bin, FUN=mean)
        sd.y <- tapply(f_dispersion, data.x.bin, FUN=sd)
        f_dispersion.scaled <-
            (f_dispersion - mean.y[as.numeric(data.x.bin)]) /
            sd.y[as.numeric(data.x.bin)]
        names(f_dispersion.scaled) <- names(f_mean)
        hvf <- data.frame(f_mean, f_dispersion, f_dispersion.scaled)
        colnames(hvf) <- paste0('mvp.',
            c('mean', 'dispersion', 'dispersion.scaled'))
    }

    # output
    rownames(hvf) <- rownames(object)
    hvf
}

