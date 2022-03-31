
putLatentsSeurat <- function(seurat_object, latent){
  
  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)
  
  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)
  
  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)
  
  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))
  
  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}

getUnwantedGenes <- function(object){
  
  unwanted_variation <- c(grep("^LINC", rownames(object), value = T), grep("^AC", rownames(object), value = T),
                          grep("^AL", rownames(object), value = T),
                          grep("^MT-", rownames(object), value = T), grep("^RP", rownames(object), value = T))
  
}


getClonalityGenes <- function(object){
  
  clonality_genes <- c(grep("^TRAV", rownames(object), value = T), grep("^TRBV", rownames(object), value = T),
                       grep("^TRGV", rownames(object), value = T), grep("^TRDV", rownames(object), value = T),
                       grep("^IGLV", rownames(object), value = T), grep("^IGLC", rownames(object), value = T),
                       grep("^IGLL", rownames(object), value = T), grep("^IGKV", rownames(object), value = T),
                       grep("^IGHV", rownames(object), value = T), grep("^IGKC", rownames(object), value = T),
                       grep("^IGH", rownames(object), value = T),  grep("^IGK", rownames(object), value = T))
  
}




getQC <- function(seurat_object){
  
  ###################
  
  min_mito     <- 0
  max_mito     <- 10
  
  min_ribo     <- 0
  max_ribo     <- 100
  
  min_features <- 200
  max_features <- 6000
  
  min_counts   <- 1000
  max_counts   <- Inf
  
  
  ###################
  
  seurat_object@meta.data$barcode <- colnames(seurat_object)
  
  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  
  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()
  
  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)
  
  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))
  
  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))
  
  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)
  
}



plotQC <- function(seurat_object, folder){
  
  # min_mito     <- 0
  # max_mito     <- 15
  # 
  # min_ribo     <- 5
  # max_ribo     <- 50
  # 
  # min_features <- 300
  # max_features <- 5e3
  # 
  # min_counts   <- 1e3
  # max_counts   <- 30e3
  
  min_mito     <- 0
  max_mito     <- 10
  
  min_ribo     <- 0
  max_ribo     <- 100
  
  min_features <- 200
  max_features <- 6000
  
  min_counts   <- 1000
  max_counts   <- Inf
  
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  qc_df$cluster = Idents(seurat_object)
  
  plotQcViolin(qc_df, var_to_plot = "nFeature_RNA", grouping = "cluster", min = min_features, max = max_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_nFeature_RNA.png"), width = 6, height = 4)
  
  plotQcViolin(qc_df, var_to_plot = "nCount_RNA", grouping = "cluster", min = min_counts, max = max_counts) + scale_y_log10() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_nCount_RNA.png"), width = 6, height = 4)
  
  plotQcViolin(qc_df, var_to_plot = "percent.mt", grouping = "cluster", min = min_mito, max = max_mito) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_percent_mt.png"), width = 6, height = 4)
  
  plotQcViolin(qc_df, var_to_plot = "percent.ribo", grouping = "cluster", min = min_ribo, max = max_ribo) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "/violin_percent_ribo.png"), width = 6, height = 4)
  
  ## Scatter plots
  p <- qc_df %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + scale_x_log10() + scale_y_log10() +
    geom_vline(xintercept = min_counts, linetype = "dotted") +
    geom_vline(xintercept = max_counts, linetype = "dotted") +
    geom_hline(yintercept = min_features, linetype = "dotted") +
    geom_hline(yintercept = max_features, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder, "/scatter_counts_vs_genes.png"), width = 5, height = 4)
  
  p <- qc_df %>%
    ggplot(aes(percent.mt, nCount_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
    scale_y_log10() +
    geom_hline(yintercept = min_counts, linetype = "dotted") +
    geom_hline(yintercept = max_counts, linetype = "dotted") +
    geom_vline(xintercept = max_mito, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder,"/scatter_mito_vs_counts.png"), width = 5, height = 4)
  
  p <- qc_df %>%
    ggplot(aes(percent.mt, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
    scale_y_log10() +
    geom_hline(yintercept = min_features, linetype = "dotted") +
    geom_hline(yintercept = max_features, linetype = "dotted") +
    geom_vline(xintercept = max_mito, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder,"/scatter_mito_vs_genes.png"), width = 5, height = 4)
  
  p <- qc_df %>%
    ggplot(aes(percent.ribo, percent.mt, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
    geom_hline(yintercept = min_mito, linetype = "dotted") +
    geom_hline(yintercept = max_mito, linetype = "dotted") +
    geom_vline(xintercept = min_ribo, linetype = "dotted") +
    geom_vline(xintercept = max_ribo, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder, "/scatter_ribo_vs_mito.png"), width = 5, height = 4)
  
}
getLatentClustering <- function(seurat_object){
  
  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)
  
}



getDEGbyClusterImmune <- function(seurat_object, cluster, min_cells = 50){
  
  message(paste0("===== ", cluster, " ====="))
  
  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= min_cells) return(NULL)
  
  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$project_temp
  
  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()
  
  n1 <- n_df %>% filter(Var1 == "AML") %>% pull(Freq) >= min_cells
  n2 <- n_df %>% filter(Var1 == "other")       %>% pull(Freq) >= min_cells
  
  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  
  cluster_markers_2v1 <- NULL
  
  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "other",  ident.2 = "AML",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(testing = "AMLvOther") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  df <- rbind(cluster_markers_2v1)
  
  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)
  
}



getDEGbyClusterBM <- function(seurat_object, cluster, min_cells = 50){
  
  message(paste0("===== ", cluster, " ====="))
  
  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= min_cells) return(NULL)
  
  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$project_temp
  
  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()
  
  n1 <- n_df %>% filter(Var1 == "AML")   %>% pull(Freq) >= min_cells
  n2 <- n_df %>% filter(Var1 == "Normal BM Petti") %>% pull(Freq) >= min_cells
  n3 <- n_df %>% filter(Var1 == "Aplastic Anemia") %>% pull(Freq) >= min_cells
  
  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE
  
  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL
  
  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "Normal BM Petti",  ident.2 = "AML",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(testing = "AMLvBM") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n1 & n3) cluster_markers_3v1 <- FindMarkers(object = seurat_cluster, ident.1 = "Aplastic Anemia",  ident.2 = "AML",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(testing = "AMLvAA") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n2 & n3) cluster_markers_3v2 <- FindMarkers(object = seurat_cluster, ident.1 = "Normal BM Petti",  ident.2 = "Aplastic Anemia",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(testing = "AAvBM") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2)
  
  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)
  
}




preprocessSeurat <- function(orig_object, cells.to.use){
  
  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)
  
  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta
  
  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)
  
  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]
  
  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))
  
  ## Scale data
  object <- ScaleData(object, features = hvg)
  
  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))
  
  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)
  
  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df
  
  return(object)
  
}




plotLolliplot <- function(viz_df, timepoint_temp){
  
  df1 <- viz_df %>% group_by(cluster, timepoint, overall) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% filter(timepoint == timepoint_temp & overall == "R")
  df2 <- viz_df %>% group_by(cluster, timepoint, overall) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% filter(timepoint == timepoint_temp & overall == "N")
  
  df_tot <- left_join(df1, df2, by = "cluster") %>% mutate(log2fc = log2(freq.x/freq.y)) %>%
    mutate(dir = ifelse(log2fc > 1, "up", "unsigf")) %>%
    mutate(dir = ifelse(log2fc < -1, "down", dir))
  
  max_y <- abs(max(df_tot$log2fc))
  
  # ggplot(df_tot, aes(cluster,log2fc,fill=dir)) + geom_bar(stat = "identity") + coord_flip() + ylim(values = c(-max_y,max_y)) + scale_fill_manual(values = c("dodgerblue", "lightgrey", "salmon"))
  
  ggplot(df_tot, aes(log2fc, reorder(cluster, log2fc), fill=dir, size = n.x)) +
    geom_segment(aes(x = 0, xend = log2fc, y=reorder(cluster, log2fc), yend = reorder(cluster, log2fc)), color = "lightgrey", size = 0.5) +
    geom_point(shape = 21) + geom_vline(xintercept = 0)  + geom_vline(xintercept = -1, linetype = "dotted") + geom_vline(xintercept = 1, linetype = "dotted") +
    xlim(values = c(-max_y,max_y)) + scale_fill_manual(values = c("dodgerblue", "lightgrey", "salmon")) + labs(y = "", size = "nCells", fill = "") + add_guide
  
}

plotLatentUmap <- function(viz_df, cluster){
  
  viz_df_temp <- data.frame(viz_df, "seurat_cluster" = cluster)
  nClusters   <- unique(cluster) %>% length
  
  ## Visualise
  latent_umap_mean <- data.frame(aggregate(latent_umap_1 ~ seurat_cluster, viz_df_temp, median), latent_umap_2 = aggregate(latent_umap_2 ~ seurat_cluster, viz_df_temp, median)[,2])
  
  
  ## Plot UMAPs with TCRGP predictions highlighted
  ggplot() +
    geom_point(data = viz_df_temp, aes(x = latent_umap_1, y = latent_umap_2, color = seurat_cluster), size = 0.8) +
    
    # stat_ellipse(data = viz_df_temp, geom = "polygon", aes(x = latent_umap_1, y = latent_umap_2, color = seurat_cluster, fill = seurat_cluster), alpha = 0.1, lty = "dotted") +
    ggrepel::geom_label_repel(data = latent_umap_mean, aes(x = latent_umap_1, y = latent_umap_2, color = seurat_cluster, label = seurat_cluster), size = 5, color = "black") +
    
    theme_void() + theme(legend.position = "none") +
    scale_color_manual(values = getPalette3(nClusters)) +
    scale_fill_manual(values = getPalette3(nClusters)) + labs()
  
  
}

getReducedNames <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    x <- strsplit(str1, "[ ]")[[1]][c(1,3)]
    p[[i]] <- paste(x, collapse = " ")
    i <- i + 1
  }
  
  return(p)
  
  
}

getNewClusters <- function(clusters){
  
  clusters %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% extractClusterNumber() %>% getClusterPhenotypes()
  
}

reorderClusters <- function(cluster_vec){
  
  ## Get clusters in order
  clusters <- cluster_vec %>% unique()
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)
  
}

extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){
  
  str1 <- substr(str1, 1, nchar(str1) - 27)
  extractFileName(str1)
}

extractTimepoint <- function(strs){
  
  strs2 <-NULL
  i <- 1
  for(str1 in strs){
    strs2[[i]] <- strsplit(str1, "[_]")[[1]][2]
    i <- i + 1
  }
  
  # return(strs2)
  return(factor(strs2, levels = c("dg", "scr", "1", "2", "3", "4", "5")))
  
}


plotQcViolin <- function(viz_df, var_to_plot, grouping, min, max){
  
  ## Plot univariate violin plots with filter thresholds
  
  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable
  
  viz_df_temp <- viz_df %>% select(var_to_plot)
  
  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table
  
  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
    
    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +
    
    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +
    
    labs(x = "", title = var_to_plot) + theme(legend.position = "none")
  
}





getClusterPhenotypesImmune <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    # "0"  = "0 CD4",
    # "1"  = "1 CD8" ,
    # "2"  = "2 NK" ,
    # "3"  = "3 CD8" ,
    # "4"  = "4 CD4" ,
    # "5"  = "5 Monocytes" ,
    # "6"  = "6 B-cells" ,
    # "7"  = "7 T cell" ,
    # "8"  = "8 CD8" ,
    # "9"  = "9 Tregs" ,
    # 
    # "10" = "10 CD8" ,
    # "11" = "11 Monocytes macrophages" ,
    # "12" = "12 Monocytes" ,
    # "13" = "13 CD4 tregs" ,
    # "14" = "14 B-cell plasma" ,
    # "15" = "15 CD8" ,
    # "16" = "16 CD4 low quality" ,
    # "17" = "17 Monocytes" ,
    # "18" = "18 CD8 low quality" ,
    # "19" = "19 B-cells" ,
    # 
    # "20" = "20 CD8" ,
    # "21" = "21 Monocytes" ,
    # "22" = "22 Keratinocytes",
    # "23" = "23 CD8 low quality",
    # "24" = "24 CD8",
    # "25" = "25 pDC",
    # "26" = "26 CD8",
    # "27" = "27 Monocytes / DC",
    # "28" = "28 Monocytes macrophages low quality",
    # "29" = "29 CD8 low quality",
    # 
    # "30" = "30 fibroblasts" ,
    # "31" = "31 Monocytes" ,
    # "32" = "32 B-cells",
    # "33" = "33 Monocytes",
    # "34" = "34 Endothelial",
    # "35" = "35 Monocytes",
    # "36" = "36 MEP",
    # "37" = "37 CD8",
    # "38" = "38 CD8",
    # "39" = "39 CD8",
    
    "0" = "0 CD8 EM" ,
    "1" = "1 CD4 CM" ,
    "2" = "2 CD4 EMRA" ,
    "3" = "3 low quality" ,
    "4" = "4 NK CD56dim" ,
    "5" = "5 B-cell" ,
    "6" = "6 CD4 tregs tumor-like" ,
    "7" = "7 CD8 EMRA cytotoxic" ,
    "8" = "8 NK CD56dim" ,
    "9" = "9 CD8 skin RM" ,
    
    "10" = "10 CD4 fh" ,
    "11" = "11 CD8 IEL-like" ,
    "12" = "12 Monocytes CD16-" ,
    "13" = "13 CD8 IEL-like" ,
    "14" = "14 Monocytes CD16+" ,
    "15" = "15 CD8 cytotoxic" ,
    "16" = "16 Monocytes CD16-" ,
    "17" = "17 DC tissue" ,
    "18" = "18 CD4 Th17" ,
    "19" = "19 CD4 Th1-like" ,
    
    "20" = "20 Monocytes" ,
    "21" = "21 CD8 RCC" ,
    "22" = "22 B-cell plasma" ,
    "23" = "23 low quality" ,
    "24" = "24 CD4 tregs blood-like" ,
    "25" = "25 CD8 aplastic anemia" ,
    "26" = "26 CD8 late exhausted" ,
    "27" = "27 CD8 stem-like exhausted" ,
    "28" = "28 CD8 EMRA" ,
    "29" = "29 T-cell cycling" ,
    
    "30" = "30 Monocytes Macrophages" ,
    "31" = "31 CD8 low quality" ,
    "32" = "32 Keratinocytes low quality" ,
    "33" = "33 Monocytes CD16-" ,
    "34" = "34 DC plasma" ,
    "35" = "35 CD8 early exhausted" ,
    "36" = "36 CD8 IFNg" ,
    "37" = "37 DC LAMP3 TIM3+" ,
    "38" = "38 CD4 Trm" ,
    "39" = "39 Monocytes macrophages low quality" ,
    
    "40" = "40 CD4" ,
    "41" = "41 Monocytes Macrophages" ,
    "42" = "42 CD8" ,
    "43" = "43 Fibroblasts low quality" ,
    "44" = "44 ILC" ,
    "45" = "45 DC CD141" ,
    "46" = "46 low quality" ,
    "47" = "47 B-cell cycling" ,
    "48" = "48 B-cell plasma" ,
    "49" = "49 GMP" ,
    
    "50" = "50 CD8 GD" ,
    "51" = "51 B-cell memory" ,
    "52" = "52 B-cell plasma low quality" ,
    "53" = "53 Endothelial low quality" ,
    "54" = "54 Monocytes Macrophage MDSC-like" ,
    "55" = "55 Monocytes Macrophage MDSC-like" ,
    "56" = "56 Monocytes low quality" ,
    "57" = "57 NK CD56bright" ,
    "58" = "58 CD8" 
    
  ))
  
  return(clusters)
  
}




getClusterPhenotypes <- function(clusters){
  
  # input : vector of clusters for FHRB1641 latent umap
  
  # clusters <- as.character(cluster)
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 NK NK",
                                                    
                                                    "1"  = "1 CD8 CD8+:effector",
                                                    "2"  = "2 CD4 CD4+:CM/naive",
                                                    "3"  = "3 Hematopoiesis MPP",
                                                    "4"  = "4 Myelopoiesis Monocyte:CD16-",
                                                    "5"  = "5 NK NKT",
                                                    "6"  = "6 NK Adaptive",
                                                    "7"  = "7 CD4 CD4+:Treg",
                                                    "8"  = "8 CD8 CD8+:CM",
                                                    "9"  = "9 Unidentified Unidentified",
                                                    "10" = "10 Hematopoiesis GMP",
                                                    
                                                    "11" = "11 Hematopoiesis HSC",
                                                    "12" = "12 Unidentified Unidentified",
                                                    "13" = "13 Hematopoiesis CMP",
                                                    "14" = "14 Unidentified Unidentified",
                                                    "15" = "15 NK doublets",
                                                    "16" = "16 Myelopoiesis GMP/monocytes",
                                                    "17" = "17 B-cell B-cell:immature",
                                                    "18" = "18 NK NK",
                                                    "19" = "19 Hematopoiesis MEP",
                                                    "20" = "20 Myelopoiesis Monocyte:CD16+",
                                                    
                                                    "21" = "21 Unidentified Unidentified",
                                                    "22" = "22 Erythropoiesis Erythroblast",
                                                    "23" = "23 Myelopoiesis Monocyte:CD16-",
                                                    "24" = "24 Other MAIT",
                                                    "25" = "25 Erythropoiesis Erythroblast",
                                                    "26" = "26 Erythropoiesis Erythroblast",
                                                    "27" = "27 B-cell pre-CD34+",
                                                    "28" = "28 NK NK",
                                                    "29" = "29 T-cell Doublets",
                                                    "30" = "30 Myelopoiesis Monocyte:CD16-",
                                                    
                                                    "31" = "31 Other Neutrophils",
                                                    "32" = "32 B-cell Plasma",
                                                    "33" = "33 Hematopoiesis StemProg",
                                                    "34" = "34 Hematopoiesis StemProg",
                                                    "35" = "35 Unidentified Unidentified",
                                                    "36" = "36 Unidentified Unidentified",
                                                    "37" = "37 Erythropoiesis Erythropoiesis",
                                                    "38" = "38 Unidentified Unidentified",
                                                    "39" = "39 B-cell B-cell:immature",
                                                    "40" = "40 B-cell B-cell:immature",
                                                    
                                                    "41" = "41 CD8 CD8+:EM",
                                                    "42" = "42 Unidentified Unidentified",
                                                    "43" = "43 Hematopoiesis GMP",
                                                    "44" = "44 Unidentified Unidentified",
                                                    "45" = "45 Other Fibroblasts"
                                                    
  ))
  
  return(clusters)
  
}

getClusterPhenotypesBlueprint <- function(clusters){
  
  # input : vector of clusters for FHRB1641 latent umap
  
  # clusters <- as.character(cluster)
  
  clusters <- plyr::revalue(clusters, replace   = c("0 " = "CD4+ Tcm"                     ,
                                                    "1 " = "CD8+ Tem"                     ,
                                                    "2 " = "GMP"                          ,
                                                    "3 " = "NK cells"                     ,
                                                    "4 " = "CD8+ Tem"                     ,
                                                    "5 " = "Monocytes"                    ,
                                                    "6 " = "Monocytes"                    ,
                                                    "7 " = "Monocytes"                    ,
                                                    "8 " = "CD8+ Tem"                     ,
                                                    "9 " = "Monocytes"                    ,
                                                    "10" = "Class-switched memory B-cells",
                                                    "11" = "GMP"                          ,
                                                    "12" = "Erythrocytes"                 ,
                                                    "13" = "Monocytes"                    ,
                                                    "14" = "Monocytes"                    ,
                                                    "15" = "Monocytes"                    ,
                                                    "16" = "Monocytes"                    ,
                                                    "17" = "CD8+ Tem"                     ,
                                                    "18" = "CD8+ Tem"                     ,
                                                    "19" = "GMP"                          ,
                                                    "20" = "Plasma cells"                 ,
                                                    "21" = "CD8+ Tem"                     ,
                                                    "22" = "MEP"                          ,
                                                    "23" = "GMP"                          ,
                                                    "24" = "Neutrophils"                  ,
                                                    "25" = "Erythrocytes"                 ,
                                                    "26" = "CD8+ Tem"                     ,
                                                    "28" = "Monocytes"                    ,
                                                    "30" = "Adipocytes"
                                                    
  ))
  
  
  return(clusters)
  
}




fixSeurat <- function(seurat_object){
  
  ## Fix meta data if it brokes
  
  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)
  
  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]
  
  rownames(meta.data) <- meta.data$barcode
  
  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)
  
  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg
  
  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)
  
}





extractCoarsePhenotype <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }
  
  return(p)
  
}




getLatentUMAP <- function(seurat_object){
  
  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df
  
  return(seurat_object)
  
}


getTimepoints <- function(cluster){
  
  cluster <- plyr::revalue(as.factor(cluster), replace =
                             
                             c("001_scr"    = "001_scr",
                               "001_161017" = "001_1"  ,
                               "001_201117" = "001_2"  ,
                               "001_C22D1"  = "001_3"  ,
                               "001_211019" = "001_4"  ,
                               "001_181119" = "001_5"  ,
                               
                               "002_dg"   = "002_dg" ,
                               "002_scr"  = "002_scr",
                               "002_C3D1" = "002_1"  ,
                               "002_EOT"  = "002_2"  ,
                               
                               "008_scr"  = "008_scr",
                               "008_C3D1" = "008_1"  ,
                               
                               "009" = "009_dg" ,
                               
                               "012_SCR"  = "012_scr",
                               "012_C3D1" = "012_1"  ,
                               
                               "014_scr"  = "014_scr",
                               "014_C3D1" = "014_1"  ,
                               "014_C6D1" = "014_2"  ,
                               
                               "017_scr"  = "017_scr",
                               "017_C3D1" = "017_1"))
  
}




getClusterPhenotypesOld <- function(clusters){
  
  # input : vector of clusters for FHRB1641 latent umap
  
  # clusters <- as.character(cluster)
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 NK Adaptive",
                                                    
                                                    "1"  = "1 Hematopoiesis CMP/GMP",
                                                    "2"  = "2 Unidentified",
                                                    "3"  = "3 CD4+ CM/naive",
                                                    "4"  = "4 Myelopoiesis monocyte CD16-",
                                                    "5"  = "5 CD8+ effector",
                                                    "6"  = "6 CD4+ Treg",
                                                    "7"  = "7 Hematopoiesis GMP",
                                                    "8"  = "8 CD8+ exhausted/EM",
                                                    "9"  = "9 Hematopoiesis CMP/GMP",
                                                    "10" = "10 Unidentifieds",
                                                    
                                                    "11" = "11 Erythropoiesis",
                                                    "12" = "12 Unidentified",
                                                    "13" = "13 NK NKT",
                                                    "14" = "14 CD8+ exhausted/RM",
                                                    "15" = "15 Myelopoiesis",
                                                    "16" = "16 B cell immature",
                                                    "17" = "17 Unidentified",
                                                    "18" = "18 Hematopoiesis CMP",
                                                    "19" = "19 Myelopoiesis monocyte CD16+",
                                                    "20" = "20 Erythropoiesis erythroblast",
                                                    
                                                    "21" = "21 Myelopoiesis monocyte CD16-",
                                                    "22" = "22 immune MAIT",
                                                    "23" = "23 Hematopoiesis GMP",
                                                    "24" = "24 Erythropoiesis",
                                                    "25" = "25 NK CD56dim1",
                                                    "26" = "26 Erythropoiesis erythroblast",
                                                    "27" = "27 Unidentified",
                                                    "28" = "28 Myelopoiesis neutrophils",
                                                    "29" = "29 Myelopoiesis monocytes",
                                                    "30" = "30 NK CD56dim2",
                                                    
                                                    "31" = "31 Hematopoiesis",
                                                    "32" = "32 other tcell_doubles",
                                                    "33" = "33 B plasma",
                                                    "34" = "34 Erythropoiesis",
                                                    "35" = "35 CD8+ Tem",
                                                    "36" = "36 Myelopoiesis pDC",
                                                    "37" = "37 Erythropoiesis",
                                                    "38" = "38 B cell",
                                                    "39" = "39 B cell immature",
                                                    "40" = "40 Unidentified",
                                                    
                                                    "41" = "41 CMP/GMP",
                                                    "42" = "42 B pro CD34+",
                                                    "43" = "43 GMP",
                                                    "44" = "44 Fibroblasts"
  ))
  
  
  return(clusters)
  
}



removeTCRabData <- function(seurat_object){
  
  seurat_object@meta.data$cdr3s_nt <- NULL
  seurat_object@meta.data$tra_cdr3s_nt <- NULL
  seurat_object@meta.data$trb_cdr3s_nt <- NULL
  seurat_object@meta.data$cdr3s_aa <- NULL
  seurat_object@meta.data$tra_cdr3s_aa <- NULL
  seurat_object@meta.data$trb_cdr3s_aa <- NULL
  seurat_object@meta.data$clonotype_id <- NULL
  seurat_object@meta.data$cdr3s_nt               <- NULL
  seurat_object@meta.data$tra_cdr3s_nt <- NULL
  seurat_object@meta.data$trb_cdr3s_nt <- NULL
  seurat_object@meta.data$cdr3s_aa <- NULL
  seurat_object@meta.data$tra_cdr3s_aa <- NULL
  seurat_object@meta.data$trb_cdr3s_aa <- NULL
  seurat_object@meta.data$clonotype_id <- NULL
  seurat_object@meta.data$frequency <- NULL
  seurat_object@meta.data$proportion <- NULL
  seurat_object@meta.data$cdr3s_aa.1 <- NULL
  seurat_object@meta.data$cdr3s_nt.1 <- NULL
  seurat_object@meta.data$chain_tra <- NULL
  seurat_object@meta.data$v_tra <- NULL
  seurat_object@meta.data$d_tra <- NULL
  seurat_object@meta.data$j_tra <- NULL
  seurat_object@meta.data$cdr3s_nt_freq_tra <- NULL
  seurat_object@meta.data$cdr3s_aa_freq_tra <- NULL
  seurat_object@meta.data$v_freq_tra <- NULL
  seurat_object@meta.data$d_freq_tra <- NULL
  seurat_object@meta.data$j_freq_tra <- NULL
  seurat_object@meta.data$chain_trb             <- NULL
  seurat_object@meta.data$v_trb <- NULL
  seurat_object@meta.data$d_trb <- NULL
  seurat_object@meta.data$j_trb <- NULL
  seurat_object@meta.data$cdr3s_nt_freq_trb <- NULL
  seurat_object@meta.data$cdr3s_aa_freq_trb <- NULL
  seurat_object@meta.data$v_freq_trb <- NULL
  seurat_object@meta.data$d_freq_trb <- NULL
  seurat_object@meta.data$j_freq_trb <- NULL
  seurat_object@meta.data$tcr_type <- NULL
  seurat_object@meta.data$patient <- NULL
  seurat_object@meta.data$new_clonotypes_id  <- NULL
  
  return(seurat_object)
  
}
