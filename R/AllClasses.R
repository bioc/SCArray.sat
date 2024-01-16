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
# Class definition

# Restrict to either dgCMatrix or SC_GDSMatrix
setClassUnion(name="UnionMatrix",
    members=c("dgCMatrix", "SC_GDSMatrix"))

# Restrict to NULL, dense matrix or SC_GDSMatrix
setClassUnion(name="UnionMatrix2",
    members=c("NULL", "matrix", "SC_GDSMatrix"))

# Extend Seurat Assay class, inherit from Seurat::Assay
setClass("SCArrayAssay", contains="Assay",
    slots = c(
        counts2 = "UnionMatrix",
        data2 = "UnionMatrix",
        scale.data2 = "UnionMatrix2"
    )
)

# Extend Seurat class
setClass("Seurat_g", contains="Seurat")

