library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(reshape2)


samples = c("Control_ST_1", "Control_ST_2", "Control_ST_3", "Control_ST_4", 
            "Control_ST_5", "BPD_ST_1", "BPD_ST_2", "BPD_ST_3", "BPD_ST_4")

stats = data.frame(Sample= samples, nspot_all=0, nspot_sel=0, 
                   rna_BDNF=0, rna_NTRK2=0, sct_BDNF=0, sct_NTRK2=0)
rownames(stats)=samples

st.list = list()
for (i in 1:length(samples)) {
  i.sample = samples[i]
  obj = NULL
  obj =  Load10X_Spatial(data.dir="../visium/", i.sample, "/outs", 
                         slice=i.sample, filter.matrix = T, to.upper = FALSE)
  obj@meta.data$Sample = i.sample
  obj@meta.data$orig.ident = obj@meta.data$i.sample
  obj <- PercentageFeatureSet(obj, "^MT-", col.name = "percent_mito")
  

  tmp_mean = mean(obj$nFeature_Spatial)
  tmp_sd = sd(obj$nFeature_Spatial)
  
  if (SAMPLE %in% c("BPD_ST_2", "BPD_ST_3")) {
    nfeatures.min = 2200 
  } else {
    nfeatures.min = tmp_mean - 1.2*tmp_sd
  }
  
  cat(round(nfeatures.min, 2), "\t", sep="")
  
  nfeatures.min = max(750, nfeatures.min)
  
  
  pmt_mean = mean(obj$percent_mito)
  pmt_sd = sd(obj$percent_mito)
  
  if (SAMPLE %in% c("BPD_ST_2", "BPD_ST_3")) {
    pmt.max = 5 
  } else {
    pmt.max = pmt_mean + 1.2*pmt_sd
  }
  
  cat(round(pmt.max, 2), "\t", sep="")
  
  stats$nspot_all[i] = dim(obj)[2]
  
  obj = obj[, obj$nFeature_Spatial >= nfeatures.min & obj$percent_mito < pmt.max]
  
  st.list[[i.sample]] = obj
  
  cat(dim(obj)[2], "\n", sep = "")
  
  stats$nspot_sel[i] = dim(obj)[2]
  
}


st.list = lapply(st.list, function(x) {
  x = SCTransform(x, assay = "Spatial", verbose = F)
  DefaultAssay(x) = "SCT"
  x
})
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = st.list, 
                                      normalization.method = "SCT",
                                      verbose = FALSE, 
                                      anchor.features = st.features)
integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                            verbose = FALSE)

# export for cell2location spot filtering
cells = integrated@meta.data
cells$barcode = gsub("_[[:print:]]+","", rownames(cells))
cells$spotid = rownames(cells)
write.table(cells, file="spots.selected.txt", sep="\t", col.names = T, row.names = F, quote=F)

saveRDS(int.anchors, file="int.anchors.rds")
saveRDS(integrated, file="integrated.rds")


