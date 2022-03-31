## Create seurat-object from individual files
folders        <- list.files("data/liao_data/rnaseq/", recursive = T, full.names = T) # %>% grep(pattern = "filtered_feature_bc_matrix", value = T) # %>% grep(pattern = "petti", value = T)
scrnaseq_files <- lapply(folders, function(x){message(x); Read10X_h5(filename = x) %>% CreateSeuratObject(project = extractFileName(x), min.cells = 3, min.features = 200)})
liao_seurat    <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractFileName(folders))

## Basic QC
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                  "TUBB", "TYMS", "UBE2C")

liao_seurat  <- PercentageFeatureSet(liao_seurat, pattern = "^MT-", col.name = "percent.mt")
liao_seurat  <- PercentageFeatureSet(liao_seurat, pattern = "^RP", col.name = "percent.ribo")
liao_seurat  <- PercentageFeatureSet(liao_seurat, features = cycle.genes, col.name = "percent.cycle")
liao_seurat@meta.data$barcode   <- colnames(liao_seurat)

liao_seurat %>% plotQC(folder = "results/qc/")
liao_seurat_qc  <- liao_seurat %>% getQC()
clonality_genes <- getClonalityGenes(liao_seurat_qc)
unwanted_genes  <- getUnwantedGenes(liao_seurat_qc)

liao_seurat_qc <- liao_seurat_qc %>% preprocessSeurat(cells.to.use = colnames(liao_seurat_qc)) 


## Add meta data
library(GEOquery)
gds          <- GEOquery::getGEO("GSE145926")
df           <- gds$GSE145926_series_matrix.txt.gz@phenoData@data
colnames(df) <- make.names(colnames(df))
df           <- df %>% dplyr::select(title,characteristics_ch1:characteristics_ch1.3,description,disease.state.ch1:tissue.ch1)
df           <- df %>% filter(grepl("scRNA", title))

df$patient   <- c("C141", "C142", "C143", "C144", "C145", "C146", "C51", "C52", "C100", "C148", "C149", "C152")
df$type      <- ifelse(grepl("COVID", df$characteristics_ch1), "COVID19", "Healthy")
df$severity  <- ifelse(grepl("severe", df$characteristics_ch1.1), "severe", "mild")
df$severity  <- ifelse(grepl("severe", df$characteristics_ch1), "severe",  df$severity)

liao_seurat_qc$patient   <- gsub("\\.h5", "", liao_seurat_qc$orig.ident)
meta_df                  <- liao_seurat_qc@meta.data %>% left_join(df, by = "patient")
rownames(meta_df)        <- meta_df$barcode
liao_seurat_qc@meta.data <- meta_df

liao_seurat_qc$type     <- ifelse(grepl("COVID", liao_seurat_qc$characteristics_ch1), "COVID19", "Healthy")
liao_seurat_qc$severity <- ifelse(grepl("severe", liao_seurat_qc$characteristics_ch1.1), "severe", "mild")
liao_seurat_qc$severity <- ifelse(grepl("severe", liao_seurat_qc$characteristics_ch1), "severe",  liao_seurat_qc$severity)

fwrite(df, "data/liao_clinical.txt", sep = "\t", quote = F, row.names = F)
saveRDS(liao_seurat, "results/liao_seurat.rds")
saveRDS(liao_seurat_qc, "results/liao_seurat_qc.rds")



## Get scVI
idents.to.keep     <- table(liao_seurat_qc$orig.ident.x) %>% as.data.frame() %>% filter(Freq > 1) %>% pull(Var1)
cells.to.keep      <- liao_seurat_qc@meta.data %>% filter(orig.ident.x %in% idents.to.keep) %>% pull(barcode)
liao_seurat_qc     <- subset(liao_seurat_qc, cells = cells.to.keep)

liao_seurat_qc_diet <- DietSeurat(liao_seurat_qc)
liao_seurat_qc_diet@assays$RNA@data <- liao_seurat_qc_diet@assays$RNA@counts
SeuratDisk::SaveH5Seurat(liao_seurat_qc_diet, filename = "results/tcr/liao_seurat_qc_diet.h5Seurat")
SeuratDisk::Convert("results/tcr/liao_seurat_qc_diet.h5Seurat", dest = "h5ad")

latents <- fread("results/scvi/liao_seurat_qc_latent.csv")
liao_seurat_qc <- liao_seurat_qc %>% putLatentsSeurat(latent = latents)
saveRDS(liao_seurat_qc, "results/liao_seurat_qc.rds")


## Select cells with detected TCR
tcr_df <- lapply(list.files("data/liao_data/tcrseq/", full.names = T), FUN = function(x){
  fread(x) %>% mutate(barcode_uniq = barcode, barcode = paste0(gsub(".csv", ".h5", extractName(extractFileName(x))), "_", barcode_uniq))
}) %>% rbindlist()

liao_seurat_qc$tcr <- colnames(liao_seurat_qc) %in% tcr_df$barcode
liao_seurat_qc_tcr <- subset(liao_seurat_qc, tcr == T)
liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% preprocessSeurat(cells.to.use = colnames(liao_seurat_qc_tcr))

## Get scVI
idents.to.keep     <- table(liao_seurat_qc_tcr$orig.ident.x) %>% as.data.frame() %>% filter(Freq > 1) %>% pull(Var1)
cells.to.keep      <- liao_seurat_qc_tcr@meta.data %>% filter(orig.ident.x %in% idents.to.keep) %>% pull(barcode)
liao_seurat_qc_tcr <- subset(liao_seurat_qc_tcr, cells = cells.to.keep)

liao_seurat_qc_tcr_diet <- DietSeurat(liao_seurat_qc_tcr)
liao_seurat_qc_tcr_diet@assays$RNA@data <- liao_seurat_qc_tcr_diet@assays$RNA@counts
SeuratDisk::SaveH5Seurat(liao_seurat_qc_tcr_diet, filename = "results/tcr/liao_seurat_qc_tcr_diet.h5Seurat")
SeuratDisk::Convert("results/tcr/liao_seurat_qc_tcr_diet.h5Seurat", dest = "h5ad")

latents            <- fread("results/scvi/liao_seurat_qc_tcr_latent.csv")
liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% putLatentsSeurat(latent = latents)
liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% getLatentClustering()

## Remove cluster 4 (one outlier cell)
Idents(liao_seurat_qc_tcr) <- liao_seurat_qc_tcr$RNA_snn_res.0.2
DimPlot(liao_seurat_qc_tcr, label = T)

cells.to.keep      <- liao_seurat_qc_tcr@meta.data %>% filter(RNA_snn_res.0.2 %in% c(0:3)) %>% pull(barcode)
liao_seurat_qc_tcr <- subset(liao_seurat_qc_tcr, cells = cells.to.keep)
liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% getLatentClustering() %>% getLatentUMAP() %>% preprocessSeurat(cells.to.use = cells.to.keep)

Idents(liao_seurat_qc_tcr) <- liao_seurat_qc_tcr$RNA_snn_res.0.2
DimPlot(liao_seurat_qc_tcr, label = T)
liao_tcr_markers <- FindAllMarkers(liao_seurat_qc_tcr, test.use = "t", max.cells.per.ident = 1e3, only.pos = T) %>% filter(p_val_adj < 0.05)
liao_tcr_markers %>% filter(gene %in% van_galen_genes)

## Remove cluster 3 (monocyte doublets)
cells.to.keep      <- liao_seurat_qc_tcr@meta.data %>% filter(RNA_snn_res.0.2 %in% c(0:2)) %>% pull(barcode)
liao_seurat_qc_tcr <- subset(liao_seurat_qc_tcr, cells = cells.to.keep)
liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% getLatentUMAP() %>% getLatentClustering() %>% preprocessSeurat(cells.to.use = colnames(liao_seurat_qc_tcr))

Idents(liao_seurat_qc_tcr) <- liao_seurat_qc_tcr$RNA_snn_res.0.2
liao_tcr_markers <- FindAllMarkers(liao_seurat_qc_tcr, test.use = "t", max.cells.per.ident = 1e3, only.pos = T) %>% filter(p_val_adj < 0.05)
liao_tcr_markers %>% filter(gene %in% big_markers)

## Get new scVI
idents.to.keep     <- table(liao_seurat_qc_tcr$orig.ident.x) %>% as.data.frame() %>% filter(Freq > 1) %>% pull(Var1)
cells.to.keep      <- liao_seurat_qc_tcr@meta.data %>% filter(orig.ident.x %in% idents.to.keep) %>% pull(barcode)
liao_seurat_qc_tcr <- subset(liao_seurat_qc_tcr, cells = cells.to.keep)

liao_seurat_qc_tcr_diet <- DietSeurat(liao_seurat_qc_tcr)
liao_seurat_qc_tcr_diet@assays$RNA@data <- liao_seurat_qc_tcr_diet@assays$RNA@counts
SeuratDisk::SaveH5Seurat(liao_seurat_qc_tcr_diet, filename = "results/tcr/liao_seurat_qc_tcr_diet_qc.h5Seurat")
SeuratDisk::Convert("results/tcr/liao_seurat_qc_tcr_diet_qc.h5Seurat", dest = "h5ad")

latents            <- fread("results/scvi/liao_seurat_qc_tcr_qclatent.csv")
liao_seurat_qc_tcr@reductions$latent <- NULL
liao_seurat_qc_tcr@reductions$latent_umap <- NULL
liao_seurat_qc_tcr@reductions$umap <- NULL

liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% putLatentsSeurat(latent = latents) %>% getLatentClustering()

Idents(liao_seurat_qc_tcr) <- liao_seurat_qc_tcr$RNA_snn_res.0.2

## Remove obvius non-CD8+ cells (cl 0 and 4)
cells.to.keep      <- liao_seurat_qc_tcr@meta.data %>% filter(RNA_snn_res.0.2 %in% c(1:3)) %>% pull(barcode)
liao_seurat_qc_tcr <- subset(liao_seurat_qc_tcr, cells = cells.to.keep)
liao_seurat_qc_tcr <- liao_seurat_qc_tcr %>% getLatentUMAP() %>% getLatentClustering() %>% preprocessSeurat(cells.to.use = colnames(liao_seurat_qc_tcr))

Idents(liao_seurat_qc_tcr) <- liao_seurat_qc_tcr$RNA_snn_res.2
a <- DimPlot(liao_seurat_qc_tcr, label = T)
b <- FeaturePlot(liao_seurat_qc_tcr, features = c("CD8A", "CD8B"), label = T)
a+b + patchwork::plot_layout(design = "ABB")

## Take only CD8+ cells (remove cl 3 and 11)
cells.to.keep          <- liao_seurat_qc_tcr@meta.data %>% filter(!RNA_snn_res.2 %in% c(3,11)) %>% pull(barcode)
liao_seurat_qc_tcr_cd8 <- subset(liao_seurat_qc_tcr, cells = cells.to.keep)
liao_seurat_qc_tcr_cd8 <- liao_seurat_qc_tcr_cd8 %>% getLatentUMAP() %>% getLatentClustering() %>% preprocessSeurat(cells.to.use = colnames(liao_seurat_qc_tcr_cd8))

liao_seurat_qc_tcr_cd8_diet <- DietSeurat(liao_seurat_qc_tcr_cd8)
liao_seurat_qc_tcr_cd8_diet@assays$RNA@data <- liao_seurat_qc_tcr_cd8_diet@assays$RNA@counts
SeuratDisk::SaveH5Seurat(liao_seurat_qc_tcr_cd8_diet, filename = "results/tcr/liao_seurat_qc_tcr_diet_qc2.h5Seurat")
SeuratDisk::Convert("results/tcr/liao_seurat_qc_tcr_diet_qc2.h5Seurat", dest = "h5ad")

latents            <- fread("results/scvi/liao_seurat_qc_tcr_qclatent2.csv")
liao_seurat_qc_tcr_cd8@reductions$latent <- NULL
liao_seurat_qc_tcr_cd8@reductions$latent_umap <- NULL
liao_seurat_qc_tcr_cd8@reductions$umap <- NULL

liao_seurat_qc_tcr_cd8 <- liao_seurat_qc_tcr_cd8 %>% putLatentsSeurat(latent = latents) %>% getLatentClustering()

Idents(liao_seurat_qc_tcr_cd8) <- liao_seurat_qc_tcr_cd8$RNA_snn_res.0.2
saveRDS(liao_seurat_qc_tcr_cd8, "results/seurat_objects/liao_seurat_qc_tcr_cd8.rds")

liao_tcr_markers <- FindAllMarkers(liao_seurat_qc_tcr_cd8, test.use = "t", max.cells.per.ident = 1e3, only.pos = T) %>% filter(p_val_adj < 0.05)
fwrite(liao_tcr_markers, "results/tcr/liao_tcr_markers.txt", sep = "\t", quote = F, row.names = F)
