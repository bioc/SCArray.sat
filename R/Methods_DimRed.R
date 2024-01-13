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


################################

# Prepare data matrix for dimensionality reduction
xPrepDR <- function(object, features=NULL, slot="scale.data", verbose=TRUE)
{
    if (length(VariableFeatures(object))==0L && is.null(features))
    {
        stop("Variable features haven't been set. ",
            "Run FindVariableFeatures() or provide a vector of feature names.")
    }
    data.use <- GetAssayData(object, slot=slot)
    if (nrow(data.use)==0L && slot=="scale.data")
        stop("Data has not been scaled. Please run ScaleData and retry")
    if (is.null(features))
        features <- VariableFeatures(object)
    features.keep <- unique(features[features %in% rownames(data.use)])
    if (length(features.keep) < length(features))
    {
        features.exclude <- setdiff(features, features.keep)
        if (verbose)
        {
            warning("The following ", length(features.exclude),
                " features requested have not been scaled (running reduction without them): ",
                paste0(features.exclude, collapse=", "))
        }
    }
    features <- features.keep
    features.var <- rowVars(data.use)
    features.keep <- features[features.var > 0]
    if (length(features.keep) < length(features))
    {
        features.exclude <- setdiff(features, features.keep)
        if (verbose)
        {
            warning("The following ", length(features.exclude),
                " features requested have zero variance (running reduction without them): ",
                paste0(features.exclude, collapse=", "))
        }
    }
    features <- features.keep
    features <- features[!is.na(features)]
    data.use[features, ]
}



####  Methods -- RunPCA()  ####

RunPCA.SCArrayAssay <- function(object, assay=NULL, features=NULL, npcs=50,
    rev.pca=FALSE, weight.by.var=TRUE, verbose=TRUE, ndims.print=1:5,
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
    weight.by.var=TRUE, verbose=TRUE, ndims.print=1:5,
    nfeatures.print=30, reduction.key="PC_", seed.use=42, approx=TRUE,
    BPPARAM, ...)
{
    x_check(object, "Calling RunPCA.SC_GDSMatrix() with %s ...")
    if (verbose)
    {
        n <- min(dim(object))
        .cat("Calculating the covariance matrix [", n, "x", n, "] ...")
        oldopt <- options(SCArray.progress.verbose=TRUE)
        on.exit(options(oldopt))
    } else {
        oldopt <- options(SCArray.progress.verbose=FALSE)
        on.exit(options(oldopt))
    }
    if (missing(BPPARAM)) BPPARAM <- getAutoBPPARAM()
    if (is.null(BPPARAM)) BPPARAM <- SerialParam()
    if (!is(BPPARAM, "BiocParallelParam"))
        stop("'BPPARAM' should be NULL or a BiocParallelParam object.")

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
            deferred=FALSE, fold=1, BPPARAM=BPPARAM)
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
            deferred=FALSE, fold=1, BPPARAM=BPPARAM)
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


####  Methods -- RunUMAP()  ####

RunUMAP.Seurat_g <- function(...)
{
    x_msg("Calling RunUMAP.Seurat_g() ...")
    if (requireNamespace("future", quietly=TRUE))
    {
        bp <- getAutoBPPARAM()
        nb <- if (is.null(bp)) 1L else bpnworkers(bp)
        if (nb != future::nbrOfWorkers())
        {
            if (is.null(bp) || is(bp, "MulticoreParam"))
            {
                if (nb == 1L)
                {
                    old_plan <- future::plan(future::sequential)
                } else if (.Platform$OS.type == "windows")
                {
                    old_plan <- future::plan(future::multisession, workers=nb)
                } else {
                    old_plan <- future::plan(future::multicore, workers=nb)
                }
                on.exit({ future::plan(old_plan) })
            }
        }
    }
    NextMethod()
}


####  Methods -- RunICA()  ####

RunICA.SCArrayAssay <- function(object, assay=NULL, features=NULL, nics=50,
    rev.ica=FALSE, ica.function="icafast", verbose=TRUE, ndims.print=1:5,
    nfeatures.print=30, reduction.name="ica", reduction.key="ica_",
    seed.use=42, ...)
{
    x_msg("Calling RunICA.SCArrayAssay() ...")
    data.use <- xPrepDR(object, features=features, verbose=verbose)
    x_msg("  \\= Running ICA dimensionality reduction ...")
    RunICA(object=data.use, assay=assay,
        nics=nics, rev.ica=rev.ica, ica.function=ica.function,
        verbose=verbose, ndims.print=ndims.print,
        nfeatures.print=nfeatures.print, reduction.key=reduction.key,
        seed.use=seed.use, ...)
}


####  Methods -- RunSPCA()  ####

RunSPCA.SCArrayAssay <- function(object, assay=NULL, features=NULL, npcs=50,
    reduction.key="SPC_", graph=NULL, verbose=TRUE, seed.use=42, ...)
{
    x_msg("Calling RunSPCA.SCArrayAssay() ...")
    data.use <- xPrepDR(object, features=features, verbose=verbose)
    x_msg("  \\= Running supervised PCA ...")
    RunSPCA(object=data.use, assay=assay,
        npcs=npcs, reduction.key=reduction.key,
        graph=graph, verbose=verbose, seed.use=seed.use, ...)
}


####  Methods -- RunLDA()  ####

RunLDA.SCArrayAssay <- function(object, assay=NULL, labels, features=NULL,
    verbose=TRUE, ndims.print=1:5, nfeatures.print=30, reduction.key="LDA_",
    seed=42, ...)
{
    x_msg("Calling RunLDA.SCArrayAssay() ...")
    data.use <- xPrepDR(object, features=features, verbose = verbose)
    x_msg("  \\= Load the matrix")
    object <- as.matrix(t(data.use))
    x_msg("  \\= Running linear discriminant analysis ...")
    RunLDA(object = object, assay=assay,
        labels=labels, verbose=verbose, ndims.print=ndims.print,
        nfeatures.print=nfeatures.print, reduction.key=reduction.key,
        seed=seed, ...)
}


####  Methods -- RunSLSI()  ####

RunSLSI.SCArrayAssay <- function(object, assay=NULL, features=NULL, n=50,
    reduction.key="SLSI_", graph=NULL, verbose=TRUE, seed.use=42, ...)
{
    x_msg("Calling RunSPCA.SCArrayAssay() ...")
    data.use <- xPrepDR(object, features=features, slot="data", verbose=verbose)
    x_msg("  \\= Running supervised LSI dimensionality reduction ...")
    RunSLSI(object=data.use, assay=assay,
        n=n, reduction.key=reduction.key, graph=graph,
        verbose=verbose, seed.use=seed.use, ...)
}


