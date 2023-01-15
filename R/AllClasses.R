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
# Class definition

setClassUnion(name="UnionMatrix",
    members=c("dgCMatrix", "SC_GDSMatrix"))

setClassUnion(name="UnionMatrix2",
    members=c("NULL", "matrix", "SC_GDSMatrix"))

setClass("SCArrayAssay", contains="Assay",
    slots = c(
        counts2 = "UnionMatrix",
        data2 = "UnionMatrix",
        scale.data2 = "UnionMatrix2"
    )
)

