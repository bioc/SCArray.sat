suppressPackageStartupMessages({
    library(RUnit)
    library(Seurat)
    library(SCArray.sat)
})

# Enable debug information in SCArray & SCArray.sat
options(SCArray.verbose=TRUE)


test_sce_matrix <- function()
{
	# row count data in a GDS file
	fn <- system.file("extdata", "example.gds", package="SCArray")

	# load Assay
	x <- scGetAssayGDS(fn)
	d1 <- Seurat::CreateSeuratObject(x)  # new DelayedMatrix-based Assay
	m <- as(GetAssayData(d1, "counts"), "sparseMatrix")
	d0 <- Seurat::CreateSeuratObject(m)  # regular in-memory Assay

    # raw counts
	m0 <- GetAssayData(d0, "counts")
	m1 <- GetAssayData(d1, "counts")
	checkEquals(m0, as(m1, "sparseMatrix"), "row counts")

	# normalize, method: CLR (margin = 1 or 2)
	d1 <- NormalizeData(d1, normalization.method="CLR", margin=1)
	d0 <- NormalizeData(d0, normalization.method="CLR", margin=1)
	m0 <- GetAssayData(d0, "data")
	m1 <- GetAssayData(d1, "data")
	checkEquals(as.matrix(m0), as.matrix(m1), "normalized counts with CLR & margin=1")
	d1 <- NormalizeData(d1, normalization.method="CLR", margin=2)
	d0 <- NormalizeData(d0, normalization.method="CLR", margin=2)
	m0 <- GetAssayData(d0, "data")
	m1 <- GetAssayData(d1, "data")
	checkEquals(as.matrix(m0), as.matrix(m1), "normalized counts with CLR & margin=2")

	# normalize, method: RC
	d1 <- NormalizeData(d1, normalization.method="RC")
	d0 <- NormalizeData(d0, normalization.method="RC")
	m0 <- GetAssayData(d0, "data")
	m1 <- GetAssayData(d1, "data")
	checkEquals(m0, as(m1, "sparseMatrix"), "normalized counts with RC")

	# normalize, method: LogNormalize
	d1 <- NormalizeData(d1)
	d0 <- NormalizeData(d0)
	m0 <- GetAssayData(d0, "data")
	m1 <- GetAssayData(d1, "data")
	checkEquals(m0, as(m1, "sparseMatrix"), "normalized counts with LogNormalize")

    # feature subsets, method: mvp
	d1 <- FindVariableFeatures(d1, nfeatures=500, selection.method="mvp")
	d0 <- FindVariableFeatures(d0, nfeatures=500, selection.method="mvp")
	checkEquals(HVFInfo(d1), HVFInfo(d0), "HVFInfo")

    # feature subsets, method: disp
	d1 <- FindVariableFeatures(d1, nfeatures=500, selection.method="disp")
	d0 <- FindVariableFeatures(d0, nfeatures=500, selection.method="disp")
	checkEquals(HVFInfo(d1), HVFInfo(d0), "HVFInfo")

    # feature subsets, default method: vst
	d1 <- FindVariableFeatures(d1, nfeatures=500)
	d0 <- FindVariableFeatures(d0, nfeatures=500)
	checkEquals(HVFInfo(d1), HVFInfo(d0), "HVFInfo")

	# scale with regressing out
	set.seed(100)
	dd <- data.frame(x1=rnorm(ncol(m1)), x2=rnorm(ncol(m1)),
		row.names=colnames(m1))
	a1 <- ScaleData(GetAssay(d1), vars.to.regress=c("x1", "x2"), latent.data=dd)
	a0 <- ScaleData(GetAssay(d0), vars.to.regress=c("x1", "x2"), latent.data=dd)
	m0 <- GetAssayData(a0, "scale.data")
	m1 <- GetAssayData(a1, "scale.data")
	checkEquals(m0, as.matrix(m1), "scaled data with regressing out the variables")

	# scale
	d1 <- ScaleData(d1)
	d0 <- ScaleData(d0)
	m0 <- GetAssayData(d0, "scale.data")
	m1 <- GetAssayData(d1, "scale.data")
	checkEquals(m0, as.matrix(m1), "scaled data")

	# runPCA
	d1 <- RunPCA(d1, verbose=FALSE)
	d0 <- RunPCA(d0, verbose=FALSE)
	m0 <- Embeddings(d0[["pca"]])
	m1 <- Embeddings(d1[["pca"]])
	for (i in 1:32)
	{
		if (mean(abs(m0[,i]+m1[,i])) < mean(abs(m0[,i]-m1[,i])))
			checkEquals(m0[,i], -m1[,i], paste("PC", i))
		else
			checkEquals(m0[,i], m1[,i], paste("PC", i))
	}

	# finally
	rm(list=ls())
	gc(FALSE)
	unlink(dir(pattern="^(_temp)?_scale.*\\.gds$"), force=TRUE)
}
