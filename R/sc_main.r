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


# Package-wide variable
.packageEnv <- new.env()


#######################################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")

.plural <- function(num) if (num > 1L) "s" else ""

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

# if getOption("SCArray.verbose")=TRUE, show message for debugging
x_check <- function(x, msg) SCArray:::x_check(x, msg)

# if getOption("SCArray.verbose")=TRUE, show message for debugging
x_msg <- function(msg) SCArray:::x_check(NULL, msg)

# write a matrix to a GDS node
x_append_gdsn <- function(mat, gdsn, verbose=TRUE)
{
    stopifnot(is(mat, "DelayedMatrix"))
    stopifnot(is(gdsn, "gdsn.class"))
    pb <- NULL
    if (verbose)
    {
        pb <- txtProgressBar(0L, ncol(mat), style=3L, width=64L, file=stderr())
        on.exit(close(pb))
    }
    # block write
    blockReduce(function(bk, i, gdsn, pb)
    {
        append.gdsn(gdsn, bk)
        if (!is.null(pb)) setTxtProgressBar(pb, i+ncol(bk))
        i + ncol(bk)
    }, mat, init=0L, grid=colAutoGrid(mat), gdsn=gdsn, pb=pb)
    # return
    invisible()
}

# return TRUE, if inefficient row enumeration is performed
x_warn_speed <- function(x, use_row=TRUE, max_bk_num=10L)
{
    stopifnot(is(x, "DelayedArray"))
    if (isTRUE(use_row))
        g <- rowAutoGrid(x)
    else
        g <- colAutoGrid(x)
    return(length(g) > max_bk_num)
}



#######################################################################

# GDS to SCArrayAssay
scNewAssayGDS <- function(gdsfile, name="counts", key="rna_", row_data=TRUE,
    check=TRUE, verbose=TRUE)
{
    # check
    x_msg("Calling scNewAssayGDS() ...")
    stopifnot(is.character(gdsfile) || inherits(gdsfile, "SCArrayFileClass"))
    stopifnot(is.character(name), length(name)==1L)
    stopifnot(is.character(key), length(key)==1L, !is.na(key))
    stopifnot(is.logical(row_data) || is.data.frame(row_data))
    stopifnot(is.logical(check), length(check)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # load gds data
    if (verbose)
    {
        s <- gdsfile
        if (inherits(s, "SCArrayFileClass")) s <- path(gdsfile)
        .cat("Input: ", s)
    }
    sce <- scExperiment(gdsfile, load.row=isTRUE(row_data))
    lst <- assays(sce)
    if (length(lst) == 0L) stop("No assay.")
    if (is.na(name)) name <- names(lst)[1L]
    m <- lst[[name]]
    if (is.null(m)) stop("No '", name, "' assay!")
    if (verbose)
        .cat("    ", name, ": ", nrow(m), " x ", ncol(m))
    # check feature IDs
    s <- rownames(m)
    if (isTRUE(check) && any(grepl('_', s)))
    {
        warning("Feature names cannot have underscores ('_'), ",
            "replacing with dashes ('-')", immediate.=TRUE)
        rownames(m) <- gsub('_', '-', s)
    }
    # meta data
    meta_data <- data.frame(row.names=rownames(m))
    if (isTRUE(row_data))
    {
        v <- rowData(sce)
        if (!is.null(v) && ncol(v)>0L)
        {
            v <- as.data.frame(v)
            if (!identical(rownames(m), rownames(v)))
            {
                stop("The rownames of 'rowData()' should be ",
                    "the same as 'count' matrix.")
            }
            meta_data <- v
        }
    }
    # key adjust if needed
    key <- Seurat:::UpdateKey(key)
    # output
    new(Class = "SCArrayAssay",
        counts2 = m, data2 = m, scale.data2 = NULL, key = key,
        meta.features = meta_data, misc = list())
}


# GDS to Seurat
scNewSeuratGDS <- function(gdsfile, assay.name=NULL, key=c(counts="rna_"),
    row_data=TRUE, col_data=TRUE, check=TRUE, verbose=TRUE)
{
    # check
    x_msg("Calling scNewSeuratGDS() ...")
    stopifnot(is.character(gdsfile) || inherits(gdsfile, "SCArrayFileClass"))
    stopifnot(is.null(assay.name) || is.character(assay.name))
    if (is.character(assay.name))
        stopifnot(length(assay.name) > 1L)
    stopifnot(is.character(key), length(key) > 0L)
    stopifnot(is.logical(row_data) || is.data.frame(row_data))
    stopifnot(is.logical(col_data) || is.data.frame(col_data))
    stopifnot(is.logical(check), length(check)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # load gds data
    if (verbose)
    {
        s <- gdsfile
        if (inherits(s, "SCArrayFileClass")) s <- path(gdsfile)
        .cat("Input: ", s)
        old_opt <- options(SCArray.progress.verbose=TRUE)
        on.exit(options(old_opt))
    }
    sce <- scExperiment(gdsfile, load.row=isTRUE(row_data),
        load.col=isTRUE(col_data))
    assay_lst <- assays(sce)
    if (length(assay_lst) == 0L) stop("No assay.")
    if (!("counts" %in% names(assay_lst)))
        stop("'counts' is not in the input GDS file!")
    if (verbose)
        .cat("    counts: ", nrow(sce), " x ", ncol(sce))

    # check assay names
    if (is.null(assay.name))
        assay.name <- names(assay_lst)
    assay.name <- unique(c("counts", assay.name[!is.na(assay.name)]))
    s <- setdiff(assay.name, names(assay_lst))
    if (length(s))
    {
        warning("No ", paste(s, collapse=", "), "!", call.=FALSE,
            immediate.=TRUE)
        assay.name <- intersect(assay.name, names(assay_lst))
    }
    if (is.na(key["counts"])) key <- c(counts="rna_", key)

    # check feature IDs
    s <- rownames(sce)
    if (isTRUE(check) && any(grepl('_', s)))
    {
        warning("Feature names cannot have underscores ('_'), ",
            "replacing with dashes ('-')", immediate.=TRUE)
        rownames(sce) <- gsub('_', '-', s)
        assay_lst <- assays(sce)
    }

    # counts
    if ("logcounts" %in% assay.name)
    {
        a <- CreateAssayObject2(counts(sce), assay_lst[["logcounts"]])
    } else {
        a <- CreateAssayObject2(counts(sce))
    }
    Key(a) <- unname(Seurat:::UpdateKey(tolower(key["counts"])))
    # meta information for features
    if (NCOL(rowData(sce)))
        a <- AddMetaData(a, as.data.frame(rowData(sce)))
    # meta information for cells
    meta.data <- NULL
    if (NCOL(colData(sce)))
        meta.data <- as.data.frame(colData(sce))
    object <- CreateSeuratObject(a, meta.data=meta.data)
    assay.name <- setdiff(assay.name, c("counts", "logcounts"))

    # Other assays
    for (nm in assay.name)
    {
        a <- CreateAssayObject2(counts=assay_lst[[nm]])
        k <- key[nm]
        if (is.na(k)) k <- paste0("rna_", nm)
        if (!grepl("_$", k)) k <- paste0(k, "_")
        Key(a) <- unname(Seurat:::UpdateKey(tolower(k)))
        object[[paste0("RNA_", nm)]] <- assay_lst[[nm]]
    }

    # output
    object
}

