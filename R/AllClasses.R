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
# Class definition

setClassUnion(name="UnionMatrix",
    members=c("dgCMatrix", "DelayedMatrix", "SC_GDSMatrix"))

setClass("SCArrayAssay", contains="Assay",
    slots = c(
        counts2 = "UnionMatrix",
        data2 = "UnionMatrix"
    )
)

