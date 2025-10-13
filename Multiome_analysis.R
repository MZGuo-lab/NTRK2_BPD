library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(SoupX)

set.seed(1234)

PREFIX="Multiome"
samples = c("Control_1","Control_2","BPD_1","BPD_3","BPD_4")

for (sample in samples) {
  
  cat(sample, "\n")
  
  counts <- Read10X_h5(paste0("../multiome_data/", sample, "_filtered_feature_bc_matrix.h5"))
  fragpath <- paste0("../multiome_data/", sample, "_atac_fragments.tsv.gz")
  
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  
  counts_RNA = counts$`Gene Expression`
  counts_ATAC = counts$Peaks
  
  obj <- CreateSeuratObject(
    counts = counts_RNA,
    project = sample,
    assay = "RNA"
  )
  
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts_ATAC,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  
  obj$barcode = rownames(obj@meta.data)
  obj$DonorID = sample
  
  print(obj)
  
  DefaultAssay(obj) <- "ATAC"
  
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  
  pMT = PercentageFeatureSet(obj, pattern="^MT-", assay="RNA")
  obj@meta.data$pMT_RNA = pMT[, 1]
  
  g1 = VlnPlot(
    object = obj,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
  )
  ggsave(file=paste0(sample, ".qc1.tiff"), plot=g1, width=16, height=6, dpi=300, units="in", compression="lzw")
  
  g2 = ggplot(obj@meta.data, aes(x=nFeature_RNA, y=pMT_RNA)) + geom_point(pt.size=0.001) + scale_x_log10() + ylim(0,100) + theme_bw()
  g2 = g2 + geom_vline(xintercept=500, col="red") + geom_hline(yintercept = 10, col="blue")
  ggsave(file=paste0(sample, ".qc2.tiff"), plot=g2, width=6, height=5, dpi=300, units="in", compression="lzw")
  
  g1 = VlnPlot(
    object = obj,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "pMT_RNA"),
    ncol = 5,
    pt.size = 0.01
  ) + theme(axis.title.x = element_blank())
  ggsave(file=paste0(sample, ".qc1.tiff"), plot=g1, width=12, height=4, dpi=300, units="in", compression="lzw")
  
  ##
  cells = NULL
  cells = subset(obj@meta.data, 
                 nCount_ATAC < 25000 &
                   nCount_RNA < 15000 &
                   nCount_ATAC > 1000 &
                   nFeature_RNA > 500 &
                   nucleosome_signal < 2 &
                   TSS.enrichment >=3 &
                   pMT_RNA<15)
  
  print(dim(cells))
  
  obj@meta.data$Selected=FALSE
  obj@meta.data[rownames(cells), "Selected"] = TRUE
  
  g2 = ggplot(obj@meta.data, aes(x=nFeature_RNA, y=pMT_RNA, col=Selected)) 
  g2 = g2 + geom_point(size=0.001) + ggtitle(paste0(dim(cells)[1], "/", dim(obj@meta.data)[1], " cells selected"))
  g2 = g2 + scale_x_log10() + ylim(0,100) + theme_bw() + theme(plot.title = element_text(color="orange"))
  g2 = g2 + geom_vline(xintercept=500, col="red")
  g2 = g2 + geom_hline(yintercept = 20, col="blue") + scale_color_manual(values=c("grey90","orange"))
  ggsave(file=paste0(sample, ".qc3.tiff"), plot=g2, width=4, height=3, dpi=300, units="in", compression="lzw")
  
  obj = subset(obj, Selected == TRUE)
  obj@meta.data = droplevels(obj@meta.data)
  saveRDS(obj, file=paste0(sample, ".rds"))
  
  
  # exprot for scrublet
  df <- data.frame(Gene=rownames(obj@assays$RNA@counts), as.matrix(obj@assays$RNA@counts), check.names = F)
  write.table(df, file=paste0("./scrublet/data/", sample,".counts.txt", sep=""), 
              sep = "\t", col.names = T, row.names = F, quote=FALSE)
  write.table(data.frame(barcode=rownames(obj@meta.data), obj@meta.data), 
              file=paste0("./scrublet/data/", sample,".meta.txt", sep=""), 
              sep = "\t", col.names = T, row.names = F, quote=FALSE)
}

# soupx
objlist = list()
grlist = list()

for (sample in samples) {
  
  raw = NULL
  i.counts = NULL
  
  i.obj = readRDS(file=paste0(sample, ".1.rds"))
  i.obj = subset(i.obj, scrublet.doublet == 0)
  
  i.cells = subset(cells, DonorID==sample)
  i.cells = droplevels(i.cells)
  i.suffix = unique(gsub("[[:print:]]+-1_","",rownames(i.cells)))
  
  i.file = paste0("../multiome_data/", sample, "_raw_feature_bc_matrix.h5")
  i.raw = Seurat::Read10X_h5(filename = i.file, use.names = T, unique.features = T)$`Gene Expression`
  i.raw = i.raw[, which(Matrix::colSums(i.raw)>0)]
  
  colnames(i.raw) = paste0(colnames(i.raw),"_",i.suffix)
  
  saveRDS(i.raw, file=paste0(sample, ".raw.gt0.rds"))
  
  i.counts = i.raw[, rownames(i.cells)]
  
  print(dim(i.counts)[2] == dim(i.cells)[1])
  
  i.obj = NormalizeData(i.obj)
  i.obj = CellCycleScoring(i.obj, s.features = Seurat::cc.genes.updated.2019$s.genes, 
                           g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
  pMT = PercentageFeatureSet(i.obj, pattern="^MT-")
  i.obj@meta.data$pMT = as.numeric(pMT[, 1])
  i.obj = SCTransform(i.obj, vars.to.regress = c("S.Score","G2M.Score","pMT"))
  i.obj = RunPCA(i.obj)
  i.obj = RunUMAP(i.obj, dims=1:30)
  i.obj = FindNeighbors(i.obj)
  i.obj = FindClusters(i.obj, method="igraph", algorithm=4)
  g = DimPlot(i.obj, reduction = "umap", group.by="seurat_clusters", label = T, repel = T, label.size=2)
  ggsave(file=paste0(sample, ".umap.seurat_clusters.tiff"), width=8, height=6, dpi=300, units="in", compression="lzw")
  
  
  sc = NULL
  out = NULL
  sc = SoupChannel(tod=i.raw, toc=i.obj@assays$RNA@counts)
  sc = setClusters(sc, setNames(i.obj@meta.data[, "seurat_clusters"], rownames(i.obj@meta.data)))
  sc = setDR(sc, DR=i.obj@reductions$umap@cell.embeddings[colnames(sc$toc), ], reductName="umap")
  tiff(file=paste0(sample, ".soupx.autoEstCont.tif"), width=6, height=5, res=300, units="in", compression="lzw")
  sc = autoEstCont(sc)
  dev.off()
  
  cat("\nCount correction\n")
  out = adjustCounts(sc)
  out = round(out)
  
  saveRDS(sc, file=paste0(sample, ".soupx.auto.sc.rds"))
  saveRDS(out, file=paste0(sample, ".soupx.auto.counts.rds"))
  
  sc = SoupChannel(tod=i.raw, toc=i.obj@assays$RNA@counts)
  sc = setClusters(sc, setNames(i.obj@meta.data[, "seurat_clusters"], rownames(i.obj@meta.data)))
  sc = setDR(sc, DR=i.obj@reductions$umap@cell.embeddings[colnames(sc$toc), ], reductName="umap")
  
  sc = setContaminationFraction(sc, 0.1)
  
  cat("\nCount correction\n")
  out = adjustCounts(sc)
  out = round(out)
  
  saveRDS(sc, file=paste0(sample, ".soupx.0.1.sc.rds"))
  saveRDS(out, file=paste0(sample, ".soupx.0.1.counts.rds"))
  
  i.obj[["SoupX"]] = CreateAssayObject(counts=out)
  
  objlist[[sample]] = i.obj
  
  grlist[[sample]] = i.obj@assays$ATAC@ranges
}


combined.peaks <- reduce(x = c(grlist$Control_1, grlist$BPD_1, grlist$BPD_3, grlist$Control_2, grlist$BPD_4))
print((combined.peaks))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


for (i in 1:length(objlist)) {
  i.sample = names(objlist)[i]
  i.obj = objlist[[i]]

  i.counts = FeatureMatrix(
    fragments = i.obj@assays$ATAC@fragments,
    features = combined.peaks,
    cells = as.character(i.obj$barcode)
  )
  i.counts = i.counts[, as.character(i.obj$barcode)]
  print(all(colnames(i.counts) == as.character(i.obj$barcode)))
  colnames(i.counts) = rownames(i.obj@meta.data)
  
  i.obj[["ATAC"]] = CreateChromatinAssay(i.counts, 
                                         fragments = i.obj@assays$ATAC@fragments,
                                         annotation = annotation)
}

combined <- merge(
  x = objlist$Control_1,
  y = list(objlist$BPD_1, objlist$BPD_3, objlist$Control_2, objlist$BPD_4)
)

saveRDS(objlist, file="objlist.rds")
saveRDS(combined, file="combined.obj.rds")


##################
## RPCA

DefaultAssay(obj) = "SoupX"
obj = NormalizeData(obj)
obj = CellCycleScoring(obj, s.features = Seurat::cc.genes.updated.2019$s.genes, 
                       g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)

objlist <- SplitObject(obj, split.by = "DonorID")

assay.use = "SoupX"
regress_cc_mt = TRUE
objlist <- lapply(X = objlist, FUN = function(x) {
  if (regress_cc_mt==TRUE) {
    x <- SCTransform(x, assay = assay.use, vars.to.regress = c("S.Score","G2M.Score","pMT_RNA"), verbose = FALSE)
  } else {
    x <- SCTransform(x, assay = assay.use, verbose = FALSE)
  }
  
})
features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000)
objlist <- PrepSCTIntegration(object.list = objlist, anchor.features = features)

objlist <- lapply(X = objlist, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = T)
})

obj.anchors <- FindIntegrationAnchors(object.list = objlist, 
                                      anchor.features = features, reduction = 'rpca', 
                                      normalization.method = "SCT",
                                      l2.norm = TRUE)

combined <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", verbose = T)

obj = combined

rm(combined); gc()

DefaultAssay(obj) <- "integrated"

npcs=200
obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
ElbowPlot(obj, ndims=50)

obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, return.model = T, min.dist = 0.5)
obj = FindNeighbors(obj, reduction = "pca", dims=1:npcs)
obj = FindClusters(obj, graph.name = "integrated_nn", algorithm=4, resolution = 0.8)

g1 = DimPlot(obj, reduction = "umap", group.by="DonorID", pt.size=0.001)
ggsave(filename = paste0(PREFIX, ".RNA.rpca.umap.DonorID.tiff"), width=7, height=6, dpi=300, units="in", compression="lzw")

g1 = DimPlot(obj, reduction = "umap", group.by="seurat_clusters", pt.size=0.001, label=T)
ggsave(filename = paste0(PREFIX, ".RNA.rpca.umap.clusters.tiff"), width=7, height=6, dpi=300, units="in", compression="lzw")


DefaultAssay(obj)= "SoupX"
obj = NormalizeData(obj)

saveRDS(obj, file=paste0(PREFIX, ".RNA.rpca.rds"))


##
DefaultAssay(obj) = "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = 10)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj, n=50)

objlist <- SplitObject(obj, split.by = "DonorID")

objlist <- lapply(X = objlist, FUN = function(x) {

  x <- FindTopFeatures(x, min.cutoff = 10)
  x <- RunTFIDF(x)
  x <- RunSVD(x, n=50)
  
})

features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000)

integration.anchors <- FindIntegrationAnchors(
  object.list = objlist,
  anchor.features = features,
  reduction = "rlsi",
  dims = 2:50
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = obj[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50)

integrated@meta.data$seurat_clusters = obj@meta.data$seurat_clusters
p2 <- DimPlot(integrated, group.by = "DonorID")

p2 <- DimPlot(integrated, group.by = "seurat_clusters", label=T)
saveRDS(integrated, file=paste0(PREFIX, ".ATAC.integraed.rds"))

##############################

DefaultAssay(obj) = "ATAC"

# process the combined dataset
obj <- FindTopFeatures(obj, min.cutoff = 10)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj, n=50)

obj@reductions$pcaUMAP = obj@reductions$umap

obj@reductions$ilsi = integrated@reductions$integrated_lsi
obj@reductions$ilsiUMAP = integrated@reductions$umap

obj <- FindMultiModalNeighbors(
  obj, reduction.list = list("pca", "ilsi"), 
  dims.list = list(1:200, 2:50), 
)

obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
               reduction.key = "wnnUMAP_", return.model = T, 
               min.dist = 0.5, seed.use = 42)

obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 4, resolution = 0.8, verbose = FALSE)

DefaultAssay(obj)= "SoupX"
obj = NormalizeData(obj)

saveRDS(obj, file=paste0(PREFIX, ".WNN.rds"))




