



savePath <- "../hcc_2023dataset_cello/"
dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)

test_text <- "hcc_2023dataset_cello_"

library(igraph)
library(VisCello)
#scriptPath <- "/kimdata/zhuqin/HSC/Preprocess.VisCello.eht/data-raw/preprocess/scripts/"
#scriptPath <- "~/Documents/CHOP/HSC/Preprocess.VisCello.eht/data-raw/preprocess/scripts/"
scriptPath <- "data-raw/scripts/"
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))
source(paste0(scriptPath, "Load_10x_latest.R"))
source(paste0(scriptPath, "iff_selection.R"))

source(paste0(scriptPath, "read_excel.R"))
set.seed(1)

sample_dirs <- c(
    "lib1"="../230214scRNA/cellrangeroutput/lib1/outs/multi/count/raw_feature_bc_matrix/",
    "lib2"="../230214scRNA/cellrangeroutput/lib2/outs/multi/count/raw_feature_bc_matrix/"
)


latest_cds <- combine_monocle_datasets(sample_dirs, minUMI = 1000, minGenes = 500)

# Move tag counts out of the expression matrix
tag_mapping_tbls<- read_excel_allsheets("../230214scRNA/CMO list.xlsx")
tag_mapping_tbls <- lapply(1:length(tag_mapping_tbls), function(i) {
    x = tag_mapping_tbls[[i]]
    x$Lib <-names(tag_mapping_tbls)[[i]]
    return(x)
})

names(tag_mapping_tbls) <- names(sample_dirs)

tag_mtx_list <- lapply(names(sample_dirs), function(dset) {
    cur_cds = latest_cds[,latest_cds$Dataset == dset]
    cur_mtx <- exprs(cur_cds)[grep("CMO", rownames(cur_cds), value = T), ]
    return(t(cur_mtx))
})
names(tag_mtx_list) <- names(sample_dirs)


saveRDS(tag_mapping_tbls, paste0(rdsPath, test_text, "tag_mapping_tbls.rds"))
saveRDS(tag_mtx_list, paste0(rdsPath, test_text, "tag_mtx_list.rds"))

# Remove tag "genes" from global transcriptome dataset
latest_cds <- latest_cds[!grepl("CMO", rownames(latest_cds)), ]

pData(latest_cds)$Total_mRNAs <- Matrix::colSums(exprs(latest_cds))

latest_cds <- estimateSizeFactors(latest_cds)
latest_cds <- detectGenes(latest_cds, min_expr = 0.1)
FM <- exprs(latest_cds)
FM <- Matrix::t(Matrix::t(FM)/sizeFactors(latest_cds))
FM <- log1p(FM)
latest_cds@assayData$norm_exprs <- Matrix(FM, sparse=T)

# deMULTIplex2
library(deMULTIplex2)

lib_name = "lib1"
tag_mtx = tag_mtx_list[[lib_name]]
lib1_res <- demultiplexTags(tag_mtx,
                       plot.path = plotPath,
                       plot.name = paste0("demux2_", lib_name),
                       plot.diagnostics = T)
table(lib1_res$final_assign)


lib_name = "lib2"
tag_mtx = tag_mtx_list[[lib_name]]
lib2_res <- demultiplexTags(tag_mtx,
                            plot.path = plotPath,
                            plot.name = paste0("demux2_", lib_name),
                            plot.diagnostics = T)
table(lib2_res$final_assign)

assign_res <- list(
    "lib1" = lib1_res, 
    "lib2" = lib2_res
)

saveRDS(assign_res, paste0(rdsPath, test_text, "assign_res.rds"))


assign_cbn <- c(assign_res$lib1$final_assign, assign_res$lib2$final_assign)
latest_cds$tag_assign <- assign_cbn[colnames(latest_cds)]
latest_cds$tag_assign[is.na(latest_cds$tag_assign)] = "negative"
tag_mapping_cbn <- do.call(rbind, tag_mapping_tbls)

tag_mapping_cbn$Tag_Lib <- paste0(tag_mapping_cbn$Tag, "_", tag_mapping_cbn$Lib)
latest_cds$Lib <- latest_cds$Dataset; latest_cds$Dataset <- NULL
latest_cds$Tag_Lib <- paste0(latest_cds$tag_assign, "_", latest_cds$Lib)

latest_cds$Dataset <- tag_mapping_cbn$Sample[match(latest_cds$Tag_Lib, tag_mapping_cbn$Tag_Lib)]
saveRDS(latest_cds, paste0(rdsPath, "latest_cds_20230404.rds"))


eset <- new("ExpressionSet",
            assayData = assayDataNew( "environment", exprs=exprs(latest_cds), norm_exprs = latest_cds@assayData$norm_exprs),
            phenoData =  new("AnnotatedDataFrame", data = pData(latest_cds)),
            featureData = new("AnnotatedDataFrame", data = fData(latest_cds)))

mito_genes <- grep("^MT-", fData(eset)$gene_short_name, value = T)
mito_cds <- eset[which(fData(eset)$gene_short_name %in% mito_genes),]
mito_values <- as.matrix(exprs(mito_cds))
mito_expr_prop <- sweep(mito_values, 2, Matrix::colSums(exprs(eset)), "/")
mito_expr_prop_sum <- colSums(mito_expr_prop)
pdf(paste0(plotPath, "all_data_mitofrac.pdf"))
hist(mito_expr_prop_sum, breaks=100)
dev.off()
pData(eset)$MitoFraction <- mito_expr_prop_sum
saveRDS(eset, paste0(savePath, "eset.rds"))

# Rename the sample names to standard
eset$Sample <-eset$Dataset 

sample_mapping <- c(
'CAF+Mac+Organoid 2' = 'F+M+O_2',
'CAF+Mac+Organoid 6' = 'F+M+O_6',
'CAF+Mac+Tumoroid 2' = 'F+M+T_2',
'CAF+Mac+Tumoroid 6' = 'F+M+T_6',
'CAF12' = "F_12",
'CAF5' = "F_5",
'M + CAF' = "M+F",
'M with CM from Organoid 2' = "M+CMO_2",
'M with CM from Organoid 6' = "M+CMO_6",
'M with CM from Tumoroid 2' = "M+CMT_2",
'M with CM from Tumoroid 6' = "M+CMT_6",
'Macrophage' = "M",
'O2 + M' = "M+O_2",
'O6 + M'= "M+O_6",
'Organoid 2' = "O_2",
'Organoid 6' = "O_6",
'T2 + M' = "M+T_2",
'T6 + M' = "M+T_6",
'Tumoroid 2' = "T_2",
'Tumoroid 6' = "T_6"
)

eset$Dataset  <- sample_mapping[eset$Sample]

eset$SampleType <- sapply(strsplit(eset$Dataset, "_"), function(x)x[1])
eset$Patient <- sapply(strsplit(eset$Dataset, "_"), function(x)x[2])
eset$Patient[eset$Patient == "2"] = "8"
eset$Patient[eset$Patient == "6"] = "24"

eset$Dataset <- ifelse(!is.na(eset$Patient), paste0(eset$SampleType, "_", eset$Patient), eset$SampleType)
eset$Patient[eset$Patient %in% c("5","12")] = NA
saveRDS(eset, paste0(savePath, "eset.rds"))



library(Seurat)
expr_all <- exprs(eset)
rownames(expr_all) <- fData(eset)$gene_short_name
seurat_global <- CreateSeuratObject(counts = expr_all, project = "screen1", min.cells = 0, min.features = 0)
identical(rownames(seurat_global@meta.data), rownames(pData(eset))) # must be TRUE
seurat_global@assays$RNA@meta.features <- fData(eset)
rownames(seurat_global@assays$RNA@meta.features) <- rownames(seurat_global)
seurat_global@meta.data <- cbind(seurat_global@meta.data, pData(eset))
seurat_global <- NormalizeData(seurat_global, normalization.method = "LogNormalize", scale.factor = 10000)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_global<- CellCycleScoring(seurat_global, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

saveRDS(seurat_global, paste0(rdsPath, "seurat_global.rds"))


eset$Cell_cycle_phase <- NA
eset$Cell_cycle_phase[match(colnames(seurat_global), colnames(eset))] <- as.character(seurat_global$Phase)
saveRDS(eset, paste0(savePath, "eset.rds"))


################# All data #################
test_text <- "all_cells_"
test_name <- "All cells"
oidx <- 1:ncol(eset)
cds_oidx <- eset[, oidx]
newvis <- new("Cello", name = test_name, idx=oidx)
expressed_g <- rowMeans(exprs(eset) > 0) > .1
cds_oidx <- cds_oidx[expressed_g,]
max_numpc <- 30
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=30, res = 5e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist <- list()
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))




################# All data, mito filtered #################
test_text <- "coculture23_singlets_mito_"
test_name <- "Coculture23 singlets [low mito]"
filter_std <- !is.na(eset$Dataset) & eset$MitoFraction < .2
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .01
cds_oidx <- cds_oidx[expressed_gene,]
newvis <- new("Cello", name = test_name, idx=oidx)
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = cds_oidx$Dataset, cluster_min_cell_num = 500, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]
max_numpc <- 50
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
set.seed(1)
for(n_component in c(2)){
    for(n_pc in c(10, 15, 30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC, IFG]`,k=30, res = 5e-4, louvain_path = "data-raw/scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))

saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))





test_text <- "coculture23_singlets_mito_seurat_"
test_name <- "Coculture23 singlets [low mito, seurat]"
dep_vis <- clist$`Coculture23 singlets [low mito]`
oidx <- dep_vis@idx
newvis <- new("Cello", name = test_name, idx=oidx)
seu_oidx <- seurat_global[,oidx]
seu_oidx <- FindVariableFeatures(seu_oidx, selection.method = "vst", nfeatures = 3000)
seu_oidx <- ScaleData(seu_oidx)
max_numpc <- 100
seu_oidx <- RunPCA(seu_oidx, features = VariableFeatures(object = seu_oidx), npcs=max_numpc)
#saveRDS(seu_oidx, paste0(rdsPath, test_text, "seu_oidx.rds"))

pca_proj <- as.data.frame(seu_oidx@reductions$pca@cell.embeddings)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30,100)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
#newvis@pmeta <- data.frame(Leiden_Cluster = factor(leiden_clus(newvis@proj$PCA, k = 30)))
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC]`,k=30, res = 5e-4, louvain_path = "./scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))




# cell type assignment

dep_vis <- clist$`Coculture23 singlets [low mito]`
dep_clus <- as.character(dep_vis@pmeta$Louvain_Cluster)
dep_st <- eset$SampleType[dep_vis@idx]
dep_clus[dep_clus %in% c(9,13) & grepl("F", dep_st)] <- "CAF"
dep_clus[dep_clus %in% c(11)  & grepl("M", dep_st) ] <- "NK cell"
dep_clus[dep_clus %in% c(1,10) & grepl("M", dep_st) ] <- "Macrophage"
dep_clus[dep_clus %in% c(6,4,5,12,8,7,2,3) & grepl("O", dep_st)] <- "Epithelial(Organoid)"
dep_clus[dep_clus %in% c(6,4,5,12,8,7,2,3) & grepl("T", dep_st)] <- "Epithelial(Tumoroid)"
dep_clus[!dep_clus %in% c("CAF", "NK cell", "Macrophage", "Epithelial(Organoid)", "Epithelial(Tumoroid)")] <- NA
pData(eset)$Cell_type <- NA
pData(eset)$Cell_type[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(savePath, "eset.rds"))






test_text <- "coculture23_Macrophage_seurat_"
test_name <- "Coculture23 Macrophage [low mito, seurat]"
oidx <- which(eset$Cell_type == "Macrophage")
newvis <- new("Cello", name = test_name, idx=oidx)
seu_oidx <- seurat_global[,oidx]
seu_oidx <- FindVariableFeatures(seu_oidx, selection.method = "vst", nfeatures = 3000)
seu_oidx <- ScaleData(seu_oidx)
max_numpc <- 100
seu_oidx <- RunPCA(seu_oidx, features = VariableFeatures(object = seu_oidx), npcs=max_numpc)
#saveRDS(seu_oidx, paste0(rdsPath, test_text, "seu_oidx.rds"))

pca_proj <- as.data.frame(seu_oidx@reductions$pca@cell.embeddings)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30,100)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
#newvis@pmeta <- data.frame(Leiden_Cluster = factor(leiden_clus(newvis@proj$PCA, k = 30)))
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC]`,k=30, res = 5e-4, louvain_path = "data-raw/scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))






test_text <- "coculture23_CAF_seurat_"
test_name <- "Coculture23 CAF [low mito, seurat]"
oidx <- which(eset$Cell_type == "CAF")
newvis <- new("Cello", name = test_name, idx=oidx)
seu_oidx <- seurat_global[,oidx]
seu_oidx <- FindVariableFeatures(seu_oidx, selection.method = "vst", nfeatures = 3000)
seu_oidx <- ScaleData(seu_oidx)
max_numpc <- 100
seu_oidx <- RunPCA(seu_oidx, features = VariableFeatures(object = seu_oidx), npcs=max_numpc)
#saveRDS(seu_oidx, paste0(rdsPath, test_text, "seu_oidx.rds"))

pca_proj <- as.data.frame(seu_oidx@reductions$pca@cell.embeddings)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30,100)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
#newvis@pmeta <- data.frame(Leiden_Cluster = factor(leiden_clus(newvis@proj$PCA, k = 30)))
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC]`,k=30, res = 5e-4, louvain_path = "data-raw/scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))




test_text <- "coculture23_tumoroid_seurat_"
test_name <- "Coculture23 Tumoroid [low mito, seurat]"
oidx <- which(eset$Cell_type == "Epithelial(Tumoroid)")
newvis <- new("Cello", name = test_name, idx=oidx)
seu_oidx <- seurat_global[,oidx]
seu_oidx <- FindVariableFeatures(seu_oidx, selection.method = "vst", nfeatures = 3000)
seu_oidx <- ScaleData(seu_oidx)
max_numpc <- 100
seu_oidx <- RunPCA(seu_oidx, features = VariableFeatures(object = seu_oidx), npcs=max_numpc)
#saveRDS(seu_oidx, paste0(rdsPath, test_text, "seu_oidx.rds"))

pca_proj <- as.data.frame(seu_oidx@reductions$pca@cell.embeddings)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30,100)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
#newvis@pmeta <- data.frame(Leiden_Cluster = factor(leiden_clus(newvis@proj$PCA, k = 30)))
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC]`,k=30, res = 5e-4, louvain_path = "data-raw/scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))





test_text <- "coculture23_organoid_seurat_"
test_name <- "Coculture23 Organoid [low mito, seurat]"
oidx <- which(eset$Cell_type == "Epithelial(Organoid)")
newvis <- new("Cello", name = test_name, idx=oidx)
seu_oidx <- seurat_global[,oidx]
seu_oidx <- FindVariableFeatures(seu_oidx, selection.method = "vst", nfeatures = 3000)
seu_oidx <- ScaleData(seu_oidx)
max_numpc <- 100
seu_oidx <- RunPCA(seu_oidx, features = VariableFeatures(object = seu_oidx), npcs=max_numpc)
#saveRDS(seu_oidx, paste0(rdsPath, test_text, "seu_oidx.rds"))

pca_proj <- as.data.frame(seu_oidx@reductions$pca@cell.embeddings)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30,100)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
#newvis@pmeta <- data.frame(Leiden_Cluster = factor(leiden_clus(newvis@proj$PCA, k = 30)))
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC]`,k=30, res = 5e-4, louvain_path = "data-raw/scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))


