Large-scale single-cell RNA-seq data analysis using GDS files and Seurat
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

The package extends the [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) classes and functions to support GDS files as a DelayedArray backend for data representation. It introduces a new `SCArrayAssay` class (derived from the Seurat `Assay`), which wraps raw counts, normalized expressions and scaled data matrix based on DelayedMatrix. It is designed to integrate seamlessly with the SeuratObject and Seurat packages to provide common data analysis, with the optimized algorithms for GDS data files. Compared with Seurat, SCArray.sat significantly reduces the memory usage and can be applied to very large datasets.

![**Figure 1**: Overview of the SCArray framework.](vignettes/scarray_sat.svg)


## Bioconductor

v0.99.0

Package News: [NEWS](./NEWS)


## Package Maintainer

[Xiuwen Zheng](xiuwen.zheng@abbvie.com)


## Installation

* Bioconductor repository
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SCArray.sat")
```


## Examples

```R
suppressPackageStartupMessages({
    library(Seurat)
    library(SCArray.sat)
})

# an input GDS file with raw counts
fn <- system.file("extdata", "example.gds", package="SCArray")
fn

a <- scGetAssayGDS(fn)
class(a)          # new "SCArrayAssay"
is(a, "Assay")    # TRUE

# create a Seurat object with the SCArrayAssay object
d <- CreateSeuratObject(a)

d <- NormalizeData(d)
d <- FindVariableFeatures(d, nfeatures=500)
d <- ScaleData(d)

d <- RunPCA(d)
DimPlot(d, reduction="pca")

# check the internal data matrices
GetAssayData(d, "counts")        # SC_GDSMatrix
path(GetAssayData(d, "counts"))  # the file name of count data
GetAssayData(d, "data")          # SC_GDSMatrix
GetAssayData(d, "scale.data")    # SC_GDSMatrix
```


## See Also

* [SCArray](http://www.bioconductor.org/packages/SCArray): Large-scale single-cell RNA-seq data manipulation with GDS files
* [Seurat](https://cran.r-project.org/package=Seurat): A toolkit for quality control, analysis, and exploration of single cell RNA sequencing data.

