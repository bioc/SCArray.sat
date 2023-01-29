suppressPackageStartupMessages({
    library(RUnit)
    library(Seurat)
    library(SCArray.sat)
})


test_sce_matrix <- function()
{
	# row count data in a GDS file
	fn <- system.file("extdata", "example.gds", package="SCArray")

	# load Assay
	x <- scGetAssayGDS(fn)
	d1 <- Seurat::CreateSeuratObject(x)  # DelayedMatrix-based Assay
	m <- as(GetAssayData(d1, "counts"), "sparseMatrix")
	d0 <- Seurat::CreateSeuratObject(m)  # in-memory Assay

	m0 <- GetAssayData(d0, "counts")
	m1 <- GetAssayData(d1, "counts")
	checkEquals(m0, as(m1, "sparseMatrix"), "row counts")

	# normalize
	d1 <- NormalizeData(d1)
	d0 <- NormalizeData(d0)
	m0 <- GetAssayData(d0, "data")
	m1 <- GetAssayData(d1, "data")
	checkEquals(m0, as(m1, "sparseMatrix"), "normalized counts")

    # feature subsets
	d1 <- FindVariableFeatures(d1, nfeatures=500)
	d0 <- FindVariableFeatures(d0, nfeatures=500)
	checkEquals(HVFInfo(d1), HVFInfo(d0), "HVFInfo")

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
}
