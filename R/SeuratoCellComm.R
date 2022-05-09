as.SingleCellExperiment.Seurat <- function(x, assay = NULL, ...) {
  CheckDots(...)
  if (!PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
  }
  assay <- assay %||% Assays(object = x)
  if (!all(assay %in% Assays(object = x))) {
    stop("One or more of the assays you are trying to convert is not in the Seurat object")
  }
  if (DefaultAssay(object = x) %in% assay) {
    assay <- union(DefaultAssay(object = x), assay)
  }
  experiments <- list()
  for (assayn in assay) {
    assays <- list(
      counts = GetAssayData(object = x, assay = assayn, slot = "counts"),
      logcounts = GetAssayData(object = x, assay = assayn, slot = "data")
    )
    scaledata_a <- GetAssayData(object = x, assay = assayn, slot = "scale.data")
    if (isTRUE(x = all.equal(
      target = dim(x = assays[["counts"]]),
      current = dim(x = scaledata_a))
    )) {
      assays[["scaledata"]] <- scaledata_a
    }
    assays <- assays[sapply(X = assays, FUN = nrow) != 0]
    sume <- SummarizedExperiment::SummarizedExperiment(assays = assays)
    experiments[[assayn]] <- sume
  }
  # create one single cell experiment
  sce <- as(object = experiments[[1]], Class = "SingleCellExperiment")
  orig.exp.name <- names(x = experiments[1])
  if (packageVersion(pkg = "SingleCellExperiment") >= "1.14.0") {
    SingleCellExperiment::mainExpName(sce) <- names(x = experiments[1])
  }
  if (length(x = experiments) > 1) {
    sce <- SingleCellExperiment::SingleCellExperiment(sce, altExps = experiments)
    sce <- SingleCellExperiment::swapAltExp(
      x = sce,
      name = orig.exp.name,
      saved = NULL
    )
  }
  metadata <- x[[]]
  metadata$ident <- Idents(object = x)
  SummarizedExperiment::colData(x = sce) <- S4Vectors::DataFrame(metadata)
  for (assayn in assay) {
    if (assayn != orig.exp.name) {
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = assayn,
        saved = orig.exp.name
      )
      SummarizedExperiment::rowData(x = sce) <- S4Vectors::DataFrame(x[[assayn]][[]])
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = orig.exp.name,
        saved = assayn
      )
    }
  }
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
  return(sce)
}
