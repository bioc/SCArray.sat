#######################################################################
#
# Package name: SCArray.sat
#
# Description:
#     Large-scale single-cell RNA-seq data analysis using GDS files and Seurat
#
# Copyright (C) 2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


####  Methods -- RunPCA()  ####

RunPCA.SCArrayAssay <- function(object, assay=NULL, features=NULL, npcs=50,
    rev.pca=FALSE, weight.by.var=TRUE, verbose=TRUE, ndims.print=seq_len(5),
    nfeatures.print=30, reduction.key="PC_", seed.use=42, ...)
{
    # check
    x_msg("Calling RunPCA.SCArrayAssay() ...")
    if (length(VariableFeatures(object))==0L && is.null(features))
    {
        stop("Variable features haven't been set. ",
            "Run FindVariableFeatures() or provide a vector of feature names.")
    }

    if (verbose)
        .cat("Run PCA on the scaled data matrix ...")
    data.use <- GetAssayData(object, "scale.data")
    if (NROW(data.use) == 0L)
        stop("Data has not been scaled. Please run ScaleData and retry.")

    # check features
    if (is.null(features)) features <- VariableFeatures(object)
    f_keep <- unique(features[features %in% rownames(data.use)])
    if (length(f_keep) < length(features))
    {
        f_exclude <- setdiff(features, f_keep)
        if (verbose)
        {
            warning(paste0("The following ", length(f_exclude),
                " features requested have not been scaled ",
                "(running reduction without them): ",
                paste0(f_exclude, collapse = ", ")), immediate.=TRUE)
        }
    }
    features <- f_keep
    data.use <- data.use[features, ]

    # run
    RunPCA(object = data.use,
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
    weight.by.var=TRUE, verbose=TRUE, ndims.print=seq_len(5),
    nfeatures.print=30, reduction.key="PC_", seed.use=42, approx=TRUE, ...)
{
    x_check(object, "Calling RunPCA.SC_GDSMatrix() with %s ...")
    if (verbose)
    {
        n <- min(dim(object))
        .cat("Calculating the covariance matrix [", n, "x", n, "] ...")
        oldopt <- options(SCArray.progress.verbose=TRUE)
        on.exit(options(oldopt))
    }

    # BiocSingular SVD functions (Irlba or Exact algorithm)
    pca_func <- if (isTRUE(approx)) runIrlbaSVD else runExactSVD

    if (!is.null(seed.use)) set.seed(seed.use)
    if (rev.pca)
    {
        totvar <- sum(colVars(object))
        npcs <- min(npcs, ncol(object)-1L)
        # the input matrix has been centered
        # fold=1 to use covariance matrix
        pca_rv <- pca_func(object, k=npcs, center=FALSE, scale=FALSE,
            deferred=FALSE, fold=1)
        sdev <- pca_rv$d / sqrt(max(1L, nrow(object)-1L))
        if (weight.by.var)
        {
            f_loadings <- pca_rv$u %*% diag(pca_rv$d)
        } else {
            f_loadings <- pca_rv$u
        }
        c_embeddings <- pca_rv$v
    } else {
        totvar <- sum(rowVars(object))
        npcs <- min(npcs, nrow(object)-1L)
        # the input matrix has been centered
        # fold=1 to use covariance matrix
        pca_rv <- pca_func(t(object), k=npcs, center=FALSE, scale=FALSE,
            deferred=FALSE, fold=1)
        f_loadings <- pca_rv$v
        sdev <- pca_rv$d / sqrt(max(1L, ncol(object)-1L))
        if (weight.by.var)
        {
            c_embeddings <- pca_rv$u %*% diag(pca_rv$d)
        } else {
            c_embeddings <- pca_rv$u
        }
    }

    rownames(f_loadings) <- rownames(object)
    colnames(f_loadings) <- paste0(reduction.key, seq_len(npcs))
    rownames(c_embeddings) <- colnames(object)
    colnames(c_embeddings) <- colnames(f_loadings)
    reduction.data <- CreateDimReducObject(
        embeddings = c_embeddings,
        loadings = f_loadings,
        assay = assay,
        stdev = sdev,
        key = reduction.key,
        misc = list(total.variance = totvar)
    )
    if (verbose)
        print(reduction.data, dims=ndims.print, nfeatures=nfeatures.print)

    return(reduction.data)
}

