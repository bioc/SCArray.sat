Large-scale single-cell RNA-seq data analysis using GDS files and Seurat
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

The package extends [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) classes and functions to support GDS files as a DelayedArray backend for data representation. It provides the SCArrayAssay class (inherited from the Seurat Assay) to wrap the DelayedMatrix-based raw counts, normalized expression values and scaled data matrix. It is designed to integrate seamlessly with the SeuratObject and Seurat packages to provide the downstream data analysis, with the optimized algorithms for GDS data files.


## Bioconductor:

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
    library(SCArray.sat)
    library(Seurat)
})

```
