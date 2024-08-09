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


####  Methods -- FoldChange()  ####

FoldChange.SCArrayAssay <- function(object, cells.1, cells.2, features=NULL,
    slot="data", pseudocount.use=1, fc.name=NULL, mean.fxn=NULL, base=2,
    norm.method=NULL, ...)
{
    x_msg("Calling FoldChange.SCArrayAssay() ...")
    rv <- NextMethod()
    # make sure that the rownames are the gene names
    if (NROW(rv)) rownames(rv) <- rownames(object)
    rv
}



####  Methods -- FindMarkers()  ####

FindMarkers.SCArrayAssay <-
function(object, slot = "data", cells.1 = NULL, cells.2 = NULL,
    features = NULL, logfc.threshold = 0.25, test.use = "wilcox",
    min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE,
    max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
    min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1,
    mean.fxn = NULL, fc.name = NULL, base = 2, densify = FALSE,
    norm.method = NULL, ...)
{
    x_msg("Calling FindMarkers.SCArrayAssay() ...")
    if (is.null(pseudocount.use)) pseudocount.use <- 1
    data.slot <- ifelse(test.use %in% Seurat:::DEmethods_counts(),
        "counts", slot)
    data.use <- GetAssayData(object, data.slot)
    counts <- switch(data.slot,
        scale.data = GetAssayData(object, "counts"), numeric())

    # reset future plan
    if (requireNamespace("future", quietly=TRUE))
    {
        if (future::nbrOfWorkers() != 1L)
        {
            old_plan <- future::plan(future::sequential)
            on.exit({ future::plan(old_plan) })
        }
    }

    # get FC values
    fc.results <- FoldChange(object, slot=data.slot,
        cells.1=cells.1, cells.2=cells.2, features=features,
        pseudocount.use=pseudocount.use, mean.fxn=mean.fxn, fc.name=fc.name,
        base=base, norm.method=norm.method)

    # parameters
    pm <- list(
        object = NULL,
        slot = data.slot,
        counts = counts, cells.1 = cells.1, cells.2 = cells.2,
        features = features, logfc.threshold = logfc.threshold,
        test.use = test.use,
        min.pct = min.pct, min.diff.pct = min.diff.pct,
        verbose = FALSE,    # not use internal progression bar (e.g., pbsapply)
        only.pos = only.pos,
        max.cells.per.ident = max.cells.per.ident,
        random.seed = random.seed, latent.vars = latent.vars,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        pseudocount.use = pseudocount.use, fc.results = fc.results,
        densify = FALSE,    # ignored
        ...)

    # split data.use for faster data access
    x_msg("Calling FindMarkers.SCArrayAssay(), blocking data ...")
    gd <- scRowAutoGrid(data.use)
    rv <- blockApply(data.use, function(bk, pm)
    {
        if (is(bk, "COO_SparseArray")) bk <- as(bk, "SVT_SparseArray")
        pm$object <- bk
        pm$fc.results <- pm$fc.results[
            match(rownames(bk), rownames(pm$fc.results)), ]
        do.call(FindMarkers, pm)
    }, grid=gd, as.sparse=NA, verbose=verbose, pm=pm)
    # merge
    de.results <- do.call(rbind, rv)

    # output
    if (test.use %in% Seurat:::DEmethods_nocorrect())
    {
        de.results <- de.results[order(-de.results$power,
            -de.results[, 1]), ]
    } else {
        # need a p_val_adj
        de.results <- de.results[order(de.results$p_val, -de.results[, 1]), ]
        de.results$p_val_adj <- p.adjust(de.results$p_val,
            method="bonferroni", n=nrow(object))
    }
    de.results
}

