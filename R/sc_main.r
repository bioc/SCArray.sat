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
    if (verbose)
        pb <- txtProgressBar(0L, ncol(mat), style=3L, width=64L, file=stderr())
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
scGetAssayGDS <- function(gdsfile, name="counts", key="rna_", row_data=TRUE,
    check=TRUE, verbose=TRUE)
{
    # check
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
    sce <- scExperiment(gdsfile)
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

