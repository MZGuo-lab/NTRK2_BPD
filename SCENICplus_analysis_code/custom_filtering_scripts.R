library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
#library(zellkonverter)
library(Signac)
library(Seurat)
library(scDotPlot)
library(ggplot2)
library(RColorBrewer)
library(legendry)
library(grid)

source('custom_filtering.R')
# -------------------------
# Load Data
# -------------------------

data_path = 'outs/scplus_pipeline/Snakemake/outs/'
# TF–RE
cistromes = read_tsv(paste0(data_path,'cistromes_direct.TF_RE.tsv'))

# RE–TG
re_tg = read_tsv(paste0(data_path,"region_to_gene_adj.tsv"))
re_tg = re_tg %>%
  rename_with(.cols = -c(target, region), .fn = ~ paste0("r2g_", .x))

# TF–TG
tf_tg = read_tsv(paste0(data_path,"tf_to_gene_adj.tsv"))
tf_tg = tf_tg %>%
  rename_with(.cols = -c(TF, target), .fn = ~ paste0("tf2g_", .x))

# Seurat object
obj = readRDS('../1.convert_to_h5ad/obj.rds')
embs = obj@reductions$wnn.umap@cell.embeddings
write.table(embs,file = 'outs/scplus_pipeline/Snakemake/outs/embeddings_wnn.csv',sep=',')
# -------------------------
# Construct Triplets
# -------------------------
triplets_all = construct_triplets(cistromes, re_tg, tf_tg) %>% data.frame()
triplets_all$region = gsub(":", "-", triplets_all$region)
unique_TFs = unique(triplets_all$TF)
unique_TGs = unique(triplets_all$target)
unique_REs = unique(triplets_all$region)

# ------------------------
# Construct Feature Matrix
# -------------------------
fragments <- CreateFragmentObject('../1.convert_to_h5ad/data/fragment/subset.tsv.gz')
acc_mat = FeatureMatrix(
  fragments = fragments,
  features = StringToGRanges(unique_REs, sep = c("-", "-")),
  cells = Cells(obj),
  process_n = 5000
)
  
ATAC_new <- CreateChromatinAssay(
  counts    = acc_mat,
  fragments = fragments,
  genome    = genome(obj[["ATAC"]]),   
  annotation= Annotation(obj[["ATAC"]])  
)
  
obj[['ATAC_SCPlus']] <- ATAC_new
obj <- RunTFIDF(obj,assay = 'ATAC_SCPlus')
obj <- FindTopFeatures(obj, min.cutoff = 'q0',assay = 'ATAC_SCPlus')
obj <- RunSVD(obj,assay = 'ATAC_SCPlus')
obj$Celltype = droplevels(obj$Celltype)
saveRDS(obj, 'obj_subset_scplus.rds')

# -------------------------
# Cell-type Specific Filtering
# -------------------------
celltype_triplets = filter_triplets_by_celltype(triplets = triplets_all,
                                                obj = obj,
                                                group_use = 'Celltype_use',
                                                RNA_assay = 'SoupX',
                                                ATAC_assay = 'ATAC_SCPlus',
                                                
                                                compute_tf = TRUE, 
                                                compute_tg = TRUE, 
                                                compute_re = FALSE,
                                                
                                                tf_pct = -Inf,
                                                tf_fc = -Inf ,
                                                tf_features = NULL,
                                                default_TFs = unique_TFs,tf_pos_only = F,
                                                
                                                tg_pct = -Inf,
                                                tg_fc = -Inf,
                                                tg_features = NULL,
                                                default_TGs = unique_TGs,tg_pos_only = T,
                                                
                                                re_pct = -Inf,
                                                re_fc = -Inf,
                                                re_features = NULL,
                                                default_REs = unique_REs,re_pos_only = T,
                                                
                                                drop_zero_var = T)
# -------------------------
# Global Filtering
# -------------------------
global_triplets = apply_global_filters(celltype_triplets,
                                       tf2g_score_min = 0.01,
                                       tf2g_rho_min = 0.01,
                                       r2g_score_min = 0.01,
                                       r2g_rho_min = 0.01)

# Save final result
saveRDS(global_triplets, "outs/scplus_pipeline/Snakemake/outs/eRegulons_custom.rds")
write.table(global_triplets, "outs/scplus_pipeline/Snakemake/outs/eRegulons_custom.txt",sep='\t',quote = F,row.names = F)
write.table(global_triplets,
            "outs/scplus_pipeline/Snakemake/outs/eRegulons_custom.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)







######Filtering based on ScenicPlus eRegulons
library(stringr)
sc_reg = read_tsv('outs/scplus_pipeline/Snakemake/outs/eRegulon_direct.tsv')
# df is your data frame
m <- str_match(sc_reg$Region, "^(chr[^:]+):(\\d+)(?:-(\\d+))?$")
chr   <- m[,2]
start <- as.integer(m[,3])
end   <- ifelse(is.na(m[,4]), start + 1L, as.integer(m[,4]))

sc_reg$Region <- sprintf("%s-%d-%d", chr, start, end)
obj = readRDS('obj_subset_scplus.rds')
unique_TFs = unique(sc_reg$TF)
unique_TGs = unique(sc_reg$Gene)
unique_REs = unique(sc_reg$Region)


# tmp_sc_reg <- dplyr::slice_sample(sc_reg, n = 100)
# unique_TFs = unique(tmp_sc_reg$TF)
# unique_TGs = unique(tmp_sc_reg$Gene)
# unique_REs = unique(tmp_sc_reg$Region)
celltype_sc_reg = filter_triplets_by_celltype(triplets = sc_reg,
                                                obj = obj,
                                                group_use = 'Celltype_use',
                                                RNA_assay = 'SoupX',
                                                ATAC_assay = 'ATAC_SCPlus',
                                                
                                                compute_tf = TRUE, 
                                                compute_tg = FALSE, 
                                                compute_re = FALSE,
                                                
                                                tf_pct = -Inf,
                                                tf_fc = -Inf ,
                                                tf_features = NULL,
                                                default_TFs = unique_TFs,tf_pos_only = F, TF_col = 'TF',
                                                
                                                tg_pct = -Inf,
                                                tg_fc = -Inf,
                                                tg_features = NULL,
                                                default_TGs = unique_TGs,tg_pos_only = F, TG_col = 'Gene',
                                                
                                                re_pct = -Inf,
                                                re_fc = -Inf,
                                                re_features = NULL,
                                                default_REs = unique_REs,re_pos_only = F, RE_col = 'Region',
                                                
                                                drop_zero_var = T)


lfc_thr <- 0.5
tf_thr = 0.2
tf_pct_col <- "DE_TF_pct.1"
df2 <- celltype_sc_reg %>%
  mutate(
    # guards for NAs
    DE_TF_p_val_adj = replace_na(DE_TF_p_val_adj, 1),
    DE_TF_avg_log2FC = replace_na(DE_TF_avg_log2FC, 0),
    `DE_TF_pct.1` = replace_na(`DE_TF_pct.1`, 0),
    
    tf_pct_ok = .data[[tf_pct_col]] >= tf_thr,
    de_ok     = (DE_TF_p_val_adj < 0.05) & (abs(DE_TF_avg_log2FC) >= lfc_thr),
    
    # row-level support in THIS celltype+condition
    row_support = tf_pct_ok & de_ok
  ) %>%
  group_by(eRegulon_name) %>%
  mutate(
    # does this eRegulon pass in at least one of the 7 celltype+condition groups?
    eReg_support_any      = any(row_support, na.rm = TRUE),
    eReg_support_n_groups = n_distinct(celltype[row_support]),
    eReg_support_groups   = paste(sort(unique(celltype[row_support])), collapse = "; ")
  ) %>%
  ungroup() %>%
  mutate(
    # mark whether THIS triplet (TF–Region–Gene row) belongs to an eRegulon that passes
    triplet_pass = eReg_support_any
  )

df2 <- df2 %>%
  mutate(symbols = sub(".*_(\\S+)$", "\\1", eRegulon_name))

saveRDS(df2, "outs/scplus_pipeline/Snakemake/outs/eRegulon_direct_filtered.rds")
write.table(df2, "outs/scplus_pipeline/Snakemake/outs/eRegulon_direct_filtered.txt",sep='\t',quote = F,row.names = F)
write.table(df2,
            "outs/scplus_pipeline/Snakemake/outs/eRegulon_direct_filtered.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#####

df2 = readRDS("outs/scplus_pipeline/Snakemake/outs/eRegulon_direct_filtered.rds")
# 1) Inputs
features_in <- unique(df2[df2$triplet_pass,'TF'])$TF
p_fixed <- scDotPlot(
  object               = obj,
  features             = feature_order,     # your ordered genes
  group                = 'Celltype_use', # your ordered groups
  clusterRows = T,
  clusterColumns = T,
  scale = T,
  NumDetectedThreshold = 0.05,
  AverageThreshold     = 0,
  dotColors            = pal,
  flipPlot             = TRUE,fontSize = 15,fontFamily = 'Arial',
)
ggsave('outs/scplus_pipeline/Snakemake/outs/dotplot_selected_TF.tiff',plot = p_fixed,width = 25, height = 4,dpi = 600,compression='lzw')



# 1) Inputs
df3 = subset(df2, triplet_pass == T & symbols =='+/+')
write.table(df3, "outs/scplus_pipeline/Snakemake/outs/eRegulon_direct_filtered_positive.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
features_in <- unique(df2[df2$triplet_pass & df2$symbols == '+/+','TF'])$TF
col_order = read.table('outs/scplus_pipeline/Snakemake/outs/eregulon_labels.txt')
col_order = col_order$V1
col_order = gsub("_(?:direct|indirect)_\\+/\\+.*$", "", col_order, perl = TRUE)

obj$Celltype_use_dotplot =  sub("^(.+?)_(.+)$", "\\1 (\\2)", obj$Celltype_use)
row_order2 <- c(
  "gCap-C2 (BPD)", "gCap-C2 (Control)",
  "gCap-C1 (BPD)", "gCap-C1 (Control)",
  "aCap (BPD)", "aCap (Control)",
  "Venous (BPD)", "Venous (Control)",
  "SystemicVenous (BPD)", "SystemicVenous (Control)",
  "Lymphatics (BPD)", "Lymphatics (Control)",
  "Arterial (BPD)", "Arterial (Control)"
)

lab <- as.character(obj$Celltype_use)
lab <- trimws(lab)                      # just in case
DefaultAssay(obj) = 'SoupX'
obj$Celltype_use_dotplot = factor(obj$Celltype_use_dotplot, levels = rev(row_order2))

# build labels for the top axis
labs_top <- paste(rev(col_order), "+/+")

# make the plot
p_fixed <- scDotPlot(
  object               = obj,
  features             = rev(col_order),
  group                = "Celltype_use_dotplot",
  clusterRows          = FALSE,
  clusterColumns       = FALSE,
  scale                = TRUE,
  NumDetectedThreshold = 0.05,
  flipPlot             = TRUE,     # genes on x-axis
  fontSize             = 15,
  fontFamily           = "Arial"
)

# unwrap if scDotPlot returned a list
if (!inherits(p_fixed, "ggplot") && !is.null(p_fixed$plot)) {
  p_fixed <- p_fixed$plot
}

# left y-axis; typography and legends
p_fixed <- p_fixed +
  scale_y_discrete(position = "left") +
  theme(
    # bottom gene labels (italic)
    axis.text.x.bottom = element_text(
      family = "Arial", face = "italic", colour = "black", angle = 90, hjust = 1, vjust = 0.5
    ),
    axis.title.x  = element_text(family = "Arial", face = "plain", colour = "black"),
    axis.title.y  = element_text(family = "Arial", face = "plain", colour = "black"),
    axis.text.y   = element_text(family = "Arial", face = "plain", colour = "black", size = 17),
    legend.title  = element_text(family = "Arial", face = "plain", colour = "black"),
    legend.text   = element_text(family = "Arial", face = "plain", colour = "black"),
    strip.text    = element_text(family = "Arial", face = "plain", colour = "black")
  )

# add secondary top axis with eRegulon labels
p_fixed2 <- p_fixed +
  guides(
    x     = ggplot2::guide_axis(position = "bottom"),   # keep bottom genes
    x.sec = legendry::guide_axis_base(                  # add top eRegulons
      key = legendry::key_manual(
        aesthetic = rev(col_order),   # categories (= breaks)
        label     = labs_top          # what to show
      ),
      position = "top"
    )
  ) +
  theme(
    axis.text.x.top       = element_text(               # top labels: plain, tight, vertical
      family = "Arial", face = "plain", colour = "black",
      angle = 90, hjust = 0, vjust = 0.5
    ),
    axis.ticks.length.x.top = unit(1, "pt"),
    axis.title.x.top      = element_blank()
  )

# p_fixed2 is your final plot
p_fixed2
ggsave('outs/scplus_pipeline/Snakemake/outs/dotplot_selected_TF_positive.tiff',plot = p_fixed2,width = 20, height = 5.5,dpi = 600,compression='lzw')