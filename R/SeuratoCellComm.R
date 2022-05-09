CreateCellRouter <- function(rawdata, assay.type='RNA', min.genes, min.cells,
                             is.expr = 0){
  assay.data <- CreateAssay(rawdata, assay.type, min.genes, min.cells, is.expr);
  assay.list <- list(assay.data)
  names(assay.list) <- assay.type
  object <- new(Class = "CellRouter", assays = assay.list)
  # object@sampTab <- assay.data@sampTab
  return(object)
}


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
       
       
       
       
  cellrouterobject@assays$[[assay]]@sampTab
  SummarizedExperiment::colData(x = sce) <- S4Vectors::DataFrame(metadata)
  for (dr in FilterObjects(object = x, classes.keep = "DimReduc")) {
    assay.used <- DefaultAssay(object = x[[dr]])
    swap.exp <- assay.used %in% SingleCellExperiment::altExpNames(x = sce) & assay.used != orig.exp.name
    if (swap.exp) {
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = assay.used,
        saved = orig.exp.name
      )
    }
    SingleCellExperiment::reducedDim(x = sce, type = toupper(x = dr)) <- Embeddings(object = x[[dr]])
    if (swap.exp) {
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = orig.exp.name,
        saved = assay.used
      )
    }
  }
    object <- new(Class = "CellRouter", assays = assay.list)
  return(sce)
}
