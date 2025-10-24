library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(zellkonverter)
library(Signac)
library(Seurat)

# -------------------------
# Construct Triplets
# -------------------------
construct_triplets = function(tf_re, re_tg, tf_tg) {
  result = tf_re %>%
    inner_join(re_tg, by = "region", relationship = "many-to-many") %>%
    inner_join(tf_tg, by = c("TF", "target"), relationship = "many-to-many")
  return(result)
}

# -------------------------
# Cell-type Specific Filtering
# -------------------------
filter_triplets_by_celltype <- function(
    triplets, obj, group_use,
    RNA_assay = "RNA", ATAC_assay = "ATAC_SCPlus",
    # which blocks to run
    compute_tf = TRUE, compute_tg = TRUE, compute_re = TRUE,
    # TF params
    tf_pct = -Inf, tf_fc = -Inf, tf_pos_only = TRUE, tf_features = NULL, default_TFs = NULL, TF_col = 'TF',
    # TG params
    tg_pct = -Inf, tg_fc = -Inf, tg_pos_only = TRUE, tg_features = NULL, default_TGs = NULL, TG_col = 'TG',
    # RE params
    re_pct = -Inf, re_fc = -Inf, re_pos_only = TRUE, re_features = NULL, default_REs = NULL, RE_col = 'RE',
    # extras
    drop_zero_var = TRUE
) {
  stopifnot(group_use %in% colnames(obj@meta.data))
  
  # Helper: choose+sanitize features against an assay
  .features_present <- function(features_in, default_features, assay, drop_zero_var = TRUE) {
    feats <- if (is.null(features_in) || length(features_in) == 0 || all(is.na(features_in))) {
      default_features
    } else {
      features_in
    }
    feats <- unique(as.character(feats))
    if (is.null(feats) || length(feats) == 0) return(character(0))
    
    present <- intersect(feats, rownames(obj[[assay]]))
    if (length(present) == 0) return(character(0))
    
    if (drop_zero_var) {
      mtx <- GetAssayData(obj, assay = assay, slot = "data")
      nz  <- Matrix::rowSums(mtx[present, , drop = FALSE] != 0) > 0
      present <- present[nz]
    }
    present
  }
  
  # normalize RNA once if needed
  if (compute_tf || compute_tg) {
    DefaultAssay(obj) <- RNA_assay
    obj <- NormalizeData(obj, assay = RNA_assay)
  }
  
  # set identities by the grouping column
  Idents(obj) <- group_use
  cell_types <- unique(obj@meta.data[[group_use]])
  filtered_all <- list()
  
  for (ct in cell_types) {
    cat(paste0('---------------------------Calculating ', ct,'--------------------------------\n'))
    
    # --- DAR (RE) ---
    dar <- NULL
    if (compute_re) {
      cat('Calculating DAR\n')
      DefaultAssay(obj) <- ATAC_assay
      
      present_re <- .features_present(
        features_in     = re_features,
        default_features = default_REs,
        assay           = ATAC_assay,
        drop_zero_var   = drop_zero_var
      )
      
      if (length(present_re) > 0) {
        dar <- FindMarkers(
          object          = obj,
          ident.1         = ct,
          test.use        = "LR",
          latent.vars     = paste0("nCount_", ATAC_assay),
          min.pct         = re_pct,
          logfc.threshold = re_fc,
          group.by        = group_use,
          assay           = ATAC_assay,
          features        = present_re,
          only.pos        = re_pos_only
        )
        if (nrow(dar) > 0) {
          colnames(dar) <- paste0('DAR_', colnames(dar))
          dar$celltype  <- ct
          dar[,RE_col]    <- rownames(dar)
        }
      } else {
        message("RE: no valid peaks present in assay after filtering; skipping DAR for ", ct)
      }
    }
    
    # --- DEG for TF ---
    tf_deg <- NULL
    if (compute_tf) {
      cat('Calculating DEGs for TF\n')
      DefaultAssay(obj) <- RNA_assay
      
      present_tf <- .features_present(
        features_in     = tf_features,
        default_features = default_TFs,
        assay           = RNA_assay,
        drop_zero_var   = drop_zero_var
      )
      
      if (length(present_tf) > 0) {
        tf_deg <- FindMarkers(
          object          = obj,
          ident.1         = ct,
          min.pct         = tf_pct,
          logfc.threshold = tf_fc,
          group.by        = group_use,
          assay           = RNA_assay,
          features        = present_tf,
          only.pos        = tf_pos_only
        )
        if (nrow(tf_deg) > 0) {
          colnames(tf_deg) <- paste0('DE_TF_', colnames(tf_deg))
          tf_deg$celltype  <- ct
          tf_deg[,TF_col]        <- rownames(tf_deg)
        }
      } else {
        message("TF: no valid genes present in assay after filtering; skipping TF DE for ", ct)
      }
    }
    
    # --- DEG for TG ---
    tg_deg <- NULL
    if (compute_tg) {
      cat('Calculating DEGs for TG\n')
      DefaultAssay(obj) <- RNA_assay
      
      present_tg <- .features_present(
        features_in     = tg_features,
        default_features = default_TGs,
        assay           = RNA_assay,
        drop_zero_var   = drop_zero_var
      )
      
      if (length(present_tg) > 0) {
        tg_deg <- FindMarkers(
          object          = obj,
          ident.1         = ct,
          min.pct         = tg_pct,
          logfc.threshold = tg_fc,
          group.by        = group_use,
          assay           = RNA_assay,
          features        = present_tg,
          only.pos        = tg_pos_only
        )
        if (nrow(tg_deg) > 0) {
          colnames(tg_deg) <- paste0('DE_TG_', colnames(tg_deg))
          tg_deg$celltype  <- ct
          tg_deg[,TG_col]  <- rownames(tg_deg)
        }
      } else {
        message("TG: no valid genes present in assay after filtering; skipping TG DE for ", ct)
      }
    }
    
    # --- join what exists ---
    tmp_result <- triplets %>%
      dplyr::mutate(celltype = ct) %>%
      { if (!is.null(dar) && nrow(dar) > 0)
        dplyr::inner_join(., dar, by = c(RE_col,"celltype"), relationship = "many-to-many")
        else . } %>%
      { if (!is.null(tf_deg) && nrow(tf_deg) > 0)
        dplyr::inner_join(., tf_deg, by = c(TF_col,"celltype"), relationship = "many-to-many")
        else . } %>%
      { if (!is.null(tg_deg) && nrow(tg_deg) > 0)
        dplyr::inner_join(., tg_deg, by = c(TG_col,"celltype"), relationship = "many-to-many")
        else . }
    
    filtered_all[[ct]] <- tmp_result
  }
  
  dplyr::bind_rows(filtered_all)
}

# -------------------------
# Global Filtering
# -------------------------
apply_global_filters = function(triplets,
                                tf2g_score_min = 0.01,
                                tf2g_rho_min = 0.01,
                                r2g_score_min = 0.01,
                                r2g_rho_min = 0.01
) {
  triplets %>%
    filter(
      tf2g_importance >= tf2g_score_min,
      abs(tf2g_rho) >= tf2g_rho_min,
      r2g_importance >= r2g_score_min,
      abs(r2g_rho) >= r2g_rho_min
    )
}