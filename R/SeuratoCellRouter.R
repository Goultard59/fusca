as.CellRouter.Seurat <- function(x, assay = NULL, ...) {
  if (!PackageCheck('Seurat', error = FALSE)) {
    stop("Please install Seurat  before converting to a CellRouter object")
  }
  assay <- assay %||% Assays(object = x)
  if (!all(assay %in% Assays(object = x))) {
    stop("One or more of the assays you are trying to convert is not in the Seurat object")
  }
  if (DefaultAssay(object = x) %in% assay) {
    assay <- union(DefaultAssay(object = x), assay)
  }
  metadata <- x[[]]
  metadata$celltypes <- Idents(object = x)
  assay.list <- list()
  for (assayn in assay) {
    if (assayn == "Spatial") {
      logcounts = GetAssayData(object = x, assay = assayn, slot = "data")
      assay.data <- CreateAssay(logcounts, "ST", min.genes, min.cells, is.expr)
      assay.list <- c(assay.list, ST=assay.data)
    }
    if (assayn == "RNA") {
      logcounts = GetAssayData(object = x, assay = assayn, slot = "data")
      assay.data <- CreateAssay(logcounts, "RNA", min.genes, min.cells, is.expr);
      assay.list <- c(assay.list, RNA=assay.data)
    }
    cellrouterobject <- new(Class = "CellRouter", assays = assay.list)
  }
  metadata <- x[[]]
  metadata$celltypes <- Idents(object = x)
  for (assay in names(assay.list) {
    cellrouterobject@assays$[[assay]]@sampTab <- metadata
  }
  cellrouterobject@var.genes <- VariableFeatures(object = x)
  cellrouterobject@pca <- list(gene.loadings = Loadings(object = x, reduction = "pca"), cell.embeddings = Embeddings(object = x, reduction = "pca"), sdev = x@reductions$pca@stdev)
  cellrouterobject@umap <- list(cell.embeddings=Embeddings(object = immune.integrated, reduction = "umap"))
  return(cellrouterobject)
}
