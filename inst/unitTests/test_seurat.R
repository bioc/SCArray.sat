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

	# load Seurat object
	d1 <- scNewSeuratGDS(fn, verbose=FALSE)  # new DelayedMatrix-based Assay5
	m <- as(LayerData(d1, layer="counts"), "sparseMatrix")
	d0 <- CreateSeuratObject(m)  # regular in-memory Assay5

	# raw counts
	m0 <- LayerData(d0, layer="counts")
	m1 <- LayerData(d1, layer="counts")
	checkEquals(m0, as(m1, "sparseMatrix"), "row counts")

	# normalize, method: CLR (margin = 1 or 2)
	d1 <- NormalizeData(d1, normalization.method="CLR", margin=1)
	d0 <- NormalizeData(d0, normalization.method="CLR", margin=1)
	m0 <- LayerData(d0, layer="data")
	m1 <- LayerData(d1, layer="data")
	checkEquals(as.matrix(m0), as.matrix(m1), "normalized counts with CLR & margin=1")
	d1 <- NormalizeData(d1, normalization.method="CLR", margin=2)
	d0 <- NormalizeData(d0, normalization.method="CLR", margin=2)
	m0 <- LayerData(d0, layer="data")
	m1 <- LayerData(d1, layer="data")
	checkEquals(as.matrix(m0), as.matrix(m1), "normalized counts with CLR & margin=2")

	# normalize, method: RC
	d1 <- NormalizeData(d1, normalization.method="RC")
	d0 <- NormalizeData(d0, normalization.method="RC")
	m0 <- LayerData(d0, layer="data")
	m1 <- LayerData(d1, layer="data")
	checkEquals(m0, as(m1, "sparseMatrix"), "normalized counts with RC")

	# normalize, method: LogNormalize
	d1 <- NormalizeData(d1)
	d0 <- NormalizeData(d0)
	m0 <- LayerData(d0, layer="data")
	m1 <- LayerData(d1, layer="data")
	checkEquals(m0, as(m1, "sparseMatrix"), "normalized counts with LogNormalize")

	# feature subsets, default method: vst
	d1 <- FindVariableFeatures(d1, nfeatures=500)
	d0 <- FindVariableFeatures(d0, nfeatures=500)
	checkEquals(unname(HVFInfo(d1, method="vst")),
		unname(HVFInfo(d0, method="vst")), "HVFInfo")

	# VariableFeatures
	checkEquals(VariableFeatures(d1), VariableFeatures(d0), "VariableFeatures")

	# scale with regressing out
	set.seed(100)
	dd <- data.frame(x1=rnorm(ncol(m1)), x2=rnorm(ncol(m1)),
		row.names=colnames(m1))
	a1 <- ScaleData(GetAssay(d1), vars.to.regress=c("x1", "x2"), latent.data=dd)
	a0 <- ScaleData(GetAssay(d0), vars.to.regress=c("x1", "x2"), latent.data=dd)
	m0 <- LayerData(a0, layer="scale.data")
	m1 <- LayerData(a1, layer="scale.data")
	checkEquals(m0, as.matrix(m1), "scaled data with regressing out covariates")

	# scale with split.by
	ss <- rep(c(TRUE, FALSE), length.out=ncol(m1))
	a1 <- ScaleData(GetAssay(d1), split.by=ss)
	a0 <- ScaleData(GetAssay(d0), split.by=ss)
	m0 <- LayerData(a0, layer="scale.data")
	m1 <- LayerData(a1, layer="scale.data")
	checkEquals(m0, as.matrix(m1), "scaled data with split.by")

	# scale
	d1 <- ScaleData(d1)
	d0 <- ScaleData(d0)
	m0 <- LayerData(d0, layer="scale.data")
	m1 <- LayerData(d1, layer="scale.data")
	checkTrue(is(m1, "SC_GDSMatrix"), "scaled data is SC_GDSMatrix")
	checkEquals(m0, as.matrix(m1), "scaled data")

	# runPCA (Seurat >= 5.6.0 uses a different PCA backend for dense matrices,
	#   so the last few PCs are not compared)
	d1 <- RunPCA(d1, verbose=FALSE)
	d0 <- RunPCA(d0, verbose=FALSE)
	checkEquals(Stdev(d0[["pca"]]), Stdev(d1[["pca"]]), "PCA stdev")
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
