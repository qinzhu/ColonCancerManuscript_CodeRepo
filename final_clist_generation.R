


library(VisCello)
savePath <- "../hcc_final_cello/"
dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)
mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

library(igraph)
#scriptPath <- "/kimdata/zhuqin/HSC/Preprocess.VisCello.eht/data-raw/preprocess/scripts/"
#scriptPath <- "~/Documents/CHOP/HSC/Preprocess.VisCello.eht/data-raw/preprocess/scripts/"
scriptPath <- "../hcc_cello/scripts/"
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))
source(paste0(scriptPath, "Load_10x_latest.R"))
source(paste0(scriptPath, "iff_selection.R"))

source(paste0(scriptPath, "read_excel.R"))
set.seed(1)

eset <- readRDS(paste0(savePath, "eset.rds"))
clist <- readRDS(paste0(savePath, "clist.rds"))


# add from old clist 
# clist_old <- readRDS("../hcc_cello/clist.rds")
# use_vis <- c("All in vivo cells [latest, IFF]", "Epithelial cells [Final, wt OGN, IFF]", "Colon cells [IFF]",  "Coculture experiment [IFF]")
# clist_idx_updated <- lapply(clist_old[use_vis], function(x) {
#     missing_total <- sum(!rownames(x@proj$PCA) %in% colnames(eset))
#     message(paste0("missing ", missing_total, " cells"))
#     use_cell <- which(rownames(x@proj$PCA) %in% colnames(eset))
#     x@proj <- lapply(x@proj, function(y) {
#         y[use_cell,]
#     })
#     x@pmeta <- x@pmeta[use_cell,, drop=F]
#     x@idx <- match(rownames(x@proj$PCA), colnames(eset))
#     return(x)
# })

clist_old <- readRDS("../hcc_final_cello/archive/20210411/clist.rds")
clist_idx_updated <- lapply(clist_old, function(x) {
    missing_total <- sum(!rownames(x@proj$PCA) %in% colnames(eset))
    message(paste0("missing ", missing_total, " cells"))
    use_cell <- which(rownames(x@proj$PCA) %in% colnames(eset))
    x@proj <- lapply(x@proj, function(y) {
        y[use_cell,]
    })
    x@pmeta <- x@pmeta[use_cell,, drop=F]
    x@idx <- match(rownames(x@proj$PCA), colnames(eset))
    return(x)
})


clist <- clist_idx_updated
saveRDS(clist, paste0(savePath,"clist.rds"))


# Cell cycle assignment
library(Seurat)
hcc_expr <- exprs(eset)
rownames(hcc_expr) <- fData(eset)$gene_short_name[match(rownames(hcc_expr), rownames(eset))]
seurat.hcc <- CreateSeuratObject(counts = hcc_expr, min.cells = 5)
seurat.hcc <- AddMetaData(seurat.hcc, metadata = pData(eset[,dep_vis@idx]))
seurat.hcc <- NormalizeData(seurat.hcc , verbose = T)
#seurat.hcc <- ScaleData(seurat.hcc, vars.to.regress="nCount_RNA")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat.hcc<- CellCycleScoring(seurat.hcc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

eset$Cell_cycle_phase <- NA
eset$Cell_cycle_phase[match(colnames(seurat.hcc), colnames(eset))] <- as.character(seurat.hcc$Phase)
saveRDS(eset, "../hcc_final_cello/eset.rds")




################# All data #################
test_text <- "all_cells"
test_name <- "All cells"
oidx <- 1:ncol(eset)
cds_oidx <- eset[, oidx]
newvis <- new("Cello", name = test_name, idx=oidx)
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .05
cds_oidx <- cds_oidx[expressed_gene,]
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = cds_oidx$Dataset, cluster_min_cell_num = 500, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]
max_numpc <- 50
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=30, res = 5e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
# clist <- list()
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))


################# All data, mito filtered #################
test_text <- "all_cells_mito_"
test_name <- "All cells [low mito]"
filter_std <- pData(eset)$MitoFraction < .2
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
newvis <- new("Cello", name = test_name, idx=oidx)
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .01
cds_oidx <- cds_oidx[expressed_gene,]
dep_clus <- clist$`All cells`@pmeta$Cluster[match(oidx, clist$`All cells`@idx)]
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = dep_clus, cluster_min_cell_num = 500, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]
max_numpc <- 50
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=30, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))


################# All in vivo, mito filtered #################
test_text <- "all_invivo_cells_20210411_"
test_name <- "All in vivo cells [20210411, IFF]"
filter_std <- eset$SampleType %in% c("Colon", "Liver", "LiverMet", "Tumor") & eset$MitoFraction < .2
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .01
cds_oidx <- cds_oidx[expressed_gene,]
newvis <- new("Cello", name = test_name, idx=oidx)
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = cds_oidx$Cell_type, cluster_min_cell_num = 500, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
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
#newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=20, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[4]],k=20, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "/Users/qinzhu/opt/anaconda3/envs/r-reticulate/bin/python", random_seed = NULL))

saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))


dep_vis <- clist$`All in vivo cells [20210220, IFF]`
dep_clus <- as.character(dep_vis@pmeta$Cluster)
dep_clus[dep_clus %in% c(1,15)] <- "Epithelial(Normal)" #
dep_clus[dep_clus %in% c(2,3,9,14,13,16)] <- "Epithelial(Tumor)" #
dep_clus[dep_clus %in% c(10)] <- "Endothelial" #
dep_clus[dep_clus %in% c(8)] <- "Fibroblast"#
dep_clus[dep_clus %in% c(12)] <- "Myofibroblast" #
dep_clus[dep_clus %in% c(5)] <- "T cell" #
dep_clus[dep_clus %in% c(19)] <- "B cell" #
dep_clus[dep_clus %in% c(4)] <- "Plasma cell" #
dep_clus[dep_clus %in% c(11,6,7)] <- "Myeloid" #
dep_clus[dep_clus %in% c(18)] <- "Mast cell" # 
dep_clus[dep_clus %in% c(17)] <- "TROP2+ Liver Progenitor" #
pData(eset)$Cell_type <- NA
pData(eset)$Cell_type[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(savePath, "eset.rds"))

# Correct for tumor epithelial
dep_vis <- clist$`All in vivo cells [20210411, IFF]`
dep_clus <- as.character(dep_vis@pmeta$Cluster)
dep_idx <- dep_vis@idx
eset$new_celltype <- eset$Cell_type
eset$new_celltype[dep_idx[dep_clus %in% c(21,6,4,10)]] <- "Epithelial(Tumor)" #
eset$Cell_type <- eset$new_celltype
eset$new_celltype <- NULL
saveRDS(eset, paste0(savePath, "eset.rds"))


# Compute co-culture experiment
test_text <- "coculture_experiment_all_"
test_name <- "Coculture experiment all [IFF]"
filter_std <- grepl("hr", eset$Dataset) & eset$MitoFraction < .2
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[2]],k=30, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
newvis <- readRDS(paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))



# Compute co-culture experiment
test_text <- "coculture_experiment_finalset_"
test_name <- "Coculture experiment final set [IFF]"
filter_std <- grepl("hr", eset$Dataset) & eset$MitoFraction < .2 & (!eset$Dataset %in% c("Tumoroid_96_24hr", "Tumoroid+Macrophage_96_24hr", "Macrophage_b2_24hr"))
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[2]],k=30, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
newvis <- readRDS(paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))

dep_vis <- clist$`Coculture experiment final set [IFF]`
dep_clus <- as.character(dep_vis@pmeta$Cluster)
dep_clus[dep_clus %in% c(1,2,5,7,8,9) & grepl("Organoid", eset$SampleType[dep_vis@idx])] <- "Epithelial(Organoid, co-culture)" 
dep_clus[dep_clus %in% c(1,2,5,7,8,9) & grepl("Tumoroid", eset$SampleType[dep_vis@idx])] <- "Epithelial(Tumoroid, co-culture)" 
dep_clus[dep_clus %in% c(3,4,6)] <- "Macrophage (co-culture)" #
dep_clus[dep_clus %in% c(1,2,5,7,8,9) & eset$SampleType[dep_vis@idx] == "Macrophage"] <- "Ambiguous"
pData(eset)$Cell_type[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(savePath, "eset.rds"))


# Compute co-culture experiment 20201216
test_text <- "coculture_experiment_20201216_"
test_name <- "Coculture experiment 20201216 [IFF]"
filter_std <- eset$Dataset %in% c("Organoid+Macrophage_86_36hr", "Organoid+Macrophage_86_36hr(2)", "Tumoroid+Macrophage_86_36hr", "Macrophage_36hr", "Organoid_86_36hr", "Tumoroid_86_36hr") & eset$MitoFraction < .2
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=20, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))

saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))



# Compute co-culture experiment 2
test_text <- "coculture_experiment_20201216_final_"
test_name <- "Coculture experiment 20201216 [IFF, final]"
dep_vis <- clist$`Coculture experiment 20201216 [IFF]`
dep_clus <- dep_vis@pmeta$Cluster
oidx <- dep_vis@idx[dep_clus != 5]
cds_oidx <- eset[, oidx]
cds_oidx <- cds_oidx[,colSums(exprs(cds_oidx)) >= 2000] # Require umi total > 2000 per cell
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .01
cds_oidx <- cds_oidx[expressed_gene,]
newvis <- new("Cello", name = test_name, idx=match(colnames(cds_oidx), colnames(eset)))
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = dep_clus[match(colnames(cds_oidx), rownames(dep_vis@proj$PCA))], cluster_min_cell_num = 500, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=20, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))




test_text <- "coculture_experiment_20201201_"
test_name <- "Coculture experiment 20201201 [IFF]"
filter_std <- grepl("92|96", eset$Dataset) | eset$Dataset %in% c("Macrophage_b1_24hr", "Macrophage_b2_24hr")  & eset$MitoFraction < .2
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=20, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))



# Compute myeloid in vivo
################# Myeloid #################
test_text <- "myeloid_cell_invivo_"
test_name <- "Myeloid cells in vivo"
dep_vis <- clist$`All in vivo cells [20210220, IFF]`
filter_std <- (1:ncol(eset)) %in% dep_vis@idx & eset$Cell_type %in% c("Myeloid")
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
dep_clus <- dep_vis@pmeta$Cluster[match(oidx, dep_vis@idx)]
newvis <- new("Cello", name = test_name, idx=oidx)
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster =dep_clus, cluster_min_cell_num = 10, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]
max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[2]],k=30, res = 1e-2, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))



# Highlight macrophages

test_text <- "macrophages_invivo_"
test_name <- "Monocytes/macrophages in vivo"
dep_vis <- clist$`Myeloid cells in vivo`
filter_std <- (1:ncol(eset)) %in% dep_vis@idx[!dep_vis@pmeta$Cluster %in% c(7,9)] & eset$Cell_type_myeloid %in% c("Monocytes/Macrophage") & as.numeric(eset$Patient) < 40
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
newvis <- new("Cello", name = test_name, idx=oidx)
cur_clus <- dep_vis@pmeta$Cluster[match(oidx, dep_vis@idx)]
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster =cur_clus, cluster_min_cell_num = 10, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]

max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
#newvis@proj[["PCA"]] <- dep_vis@proj$PCA[match_idx,]
match_idx <- match(colnames(cds_oidx), rownames(dep_vis@proj$PCA))
newvis@proj[["UMAP-2D [zoom]"]] <- dep_vis@proj$`UMAP-2D [10PC, IFG]`[match_idx,]
newvis@pmeta <- dep_vis@pmeta[match_idx,,drop=F]
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))







# Update myeloid
test_text <- "myeloid_cell_invivo_20210411_"
test_name <- "Myeloid cells in vivo [20210411]"
dep_vis <- clist$`All in vivo cells [20210411, IFF]`
filter_std <- (1:ncol(eset)) %in% dep_vis@idx & eset$Cell_type %in% c("Myeloid")
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
dep_clus <- dep_vis@pmeta$Cluster[match(oidx, dep_vis@idx)]
newvis <- new("Cello", name = test_name, idx=oidx)
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster =dep_clus, cluster_min_cell_num = 10, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]
max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
#newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[4]],k=20, res = 1e-3, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))

dep_vis <- clist$`Myeloid cells in vivo [20210411]`
dep_clus <- as.character(dep_vis@pmeta$Cluster)
dep_clus[dep_clus %in% c(4,8)] <- "cDC2 (CD1C+)"
dep_clus[dep_clus %in% "cDC2 (CD1C+)" & exprs(eset)[fData(eset)$gene_short_name == "LILRA4", dep_vis@idx] > 0] <- "pDC (LILRA4+)"
dep_clus[dep_clus %in% c(9)] <- "cDC1 (BATF3+)"
dep_clus[dep_clus %in% c(11)] <- "LAMP3+ DC (LAMP3+)"
dep_clus[dep_clus %in% c(12)] <- "Macrophage-T cell doublet"
dep_clus[dep_clus %in% c(13)] <- "Unannotated"
dep_clus[dep_clus %in% c(1,2,3,5,6,7,10)] <- "Macrophage"

pData(eset)$Cell_type_myeloid <- NA
pData(eset)$Cell_type_myeloid[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(savePath, "eset.rds"))


test_text <- "macrophages_invivo_20210411_"
test_name <- "Monocytes/macrophages in vivo [20210411]"
dep_vis <- clist$`Myeloid cells in vivo [20210411]`
filter_std <- eset$Cell_type_myeloid %in% c("Macrophage") 
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
newvis <- new("Cello", name = test_name, idx=oidx)
cur_clus <- dep_vis@pmeta$Cluster[match(oidx, dep_vis@idx)]
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster =cur_clus, cluster_min_cell_num = 10, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]

match_idx <- match(colnames(cds_oidx), rownames(dep_vis@proj$PCA))
newvis@proj[["UMAP-2D [zoom]"]] <- dep_vis@proj$`UMAP-2D [10PC, batch regressed]`[match_idx,]
newvis@pmeta <- dep_vis@pmeta[match_idx,,drop=F]
max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}

#newvis@proj[["PCA"]] <- dep_vis@proj$PCA[match_idx,]
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))








test_text <- "coculture_experiment_P92_"
test_name <- "Coculture experiment P92 [IFF]"
filter_std <- grepl("92", eset$Dataset) | eset$Dataset %in% c("Macrophage_b1_24hr", "Macrophage_b2_24hr")  & eset$MitoFraction < .2
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
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=20, res = 1e-4, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))





# Update fibroblast/myofibroblast
test_text <- "mesenchymal_cell_invivo_20221224_"
test_name <- "Mesenchymal cells in vivo [20221224]"
dep_vis <- clist$`All in vivo cells [20210411, IFF]`
eset$temp_clus <- NA; eset$temp_clus[dep_vis@idx] = dep_vis@pmeta$Cluster
filter_std <- (1:ncol(eset)) %in% dep_vis@idx & eset$Cell_type %in% c("Fibroblast","Myofibroblast") & (!eset$Patient %in% c("40") & eset$temp_clus %in% c(15,17,8))
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .05
newvis <- new("Cello", name = test_name, idx=oidx)
cds_oidx <- cds_oidx[expressed_gene,]
max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
#newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}

for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, global subset]")]]<-dep_vis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]][match(oidx, dep_vis@idx),]
    }
}

#newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[4]],k=20, res = 1e-3, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))








# Update T cells
test_text <- "T_cell_invivo_20221224_"
test_name <- "T cells in vivo [20221224]"
dep_vis <- clist$`All in vivo cells [20210411, IFF]`
eset$temp_clus <- NA; eset$temp_clus[dep_vis@idx] = dep_vis@pmeta$Cluster
filter_std <- (1:ncol(eset)) %in% dep_vis@idx & eset$Cell_type == "T cell" & (!eset$Patient %in% c("40") )
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .05
newvis <- new("Cello", name = test_name, idx=oidx)
cds_oidx <- cds_oidx[expressed_gene,]
max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
#newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}

for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, global subset]")]]<-dep_vis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]][match(oidx, dep_vis@idx),]
    }
}
# py_install("louvain", pip=TRUE, pip_ignore_installed=TRUE)
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [10PC, batch regressed]`,k=20, res = 1e-3, louvain_path = "../hcc_cello/scripts/python/louvain.py", python_home = NULL, random_seed = NULL))

saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))






test_text <- "T_cell_invivo_20221224_IFF_"
test_name <- "T cells in vivo [20221224, IFF]"
dep_vis <- clist$`T cells in vivo [20221224]`
oidx <- dep_vis@idx
cds_oidx <- eset[, oidx]
dep_clus <- dep_vis@pmeta$Louvain_Cluster[match(oidx, dep_vis@idx)]
newvis <- new("Cello", name = test_name, idx=oidx)
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster =dep_clus, cluster_min_cell_num = 10, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]

max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
#newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}

# py_install("louvain", pip=TRUE, pip_ignore_installed=TRUE)
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [30PC, batch regressed]`,k=15, res = 5e-3, louvain_path = "../hcc_cello/scripts/python/louvain.py", python_home = NULL, random_seed = NULL))

saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))


dep_vis <- clist$`T cells in vivo [20221224, IFF]`
dep_clus <- as.character(dep_vis@pmeta$Louvain_Cluster)
dep_clus[dep_clus %in% c(11)] <- "NK cell"
dep_clus[dep_clus %in% c(16)] <- "NKT cell"
dep_clus[dep_clus %in% c(17,2,14,5)] <- "FOXP3+ Treg"
dep_clus[dep_clus %in% c(4,1,7)] <- "CCR7+ Tn/Tm"
dep_clus[dep_clus %in% c(3,10,15,12)] <- "Cycling T cell"
dep_clus[dep_clus %in% c(9,18,6)] <- "GZMA/GZMB+ Teff/Tem"
dep_clus[dep_clus %in% c(8)] <- "KLRG1+ SLEC"
dep_clus[dep_clus %in% c(13)] <- "T-Carcinoma cell doublet"
pData(eset)$Cell_type_T <- "unannotated"
pData(eset)$Cell_type_T[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(savePath, "eset.rds"))








# Update B cells
test_text <- "B_cell_invivo_20221224_"
test_name <- "B cells in vivo [20221224]"
dep_vis <- clist$`All in vivo cells [20210411, IFF]`
eset$temp_clus <- NA; eset$temp_clus[dep_vis@idx] = dep_vis@pmeta$Cluster
filter_std <- (1:ncol(eset)) %in% dep_vis@idx & eset$Cell_type %in% c("B cell", "Plasma cell") & (!eset$Patient %in% c("40") )
oidx <- which(filter_std)
cds_oidx <- eset[, oidx]
expressed_gene <- rowMeans(exprs(cds_oidx) > 0) > .05
newvis <- new("Cello", name = test_name, idx=oidx)
cds_oidx <- cds_oidx[expressed_gene,]
max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
#newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}

for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, global subset]")]]<-dep_vis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]][match(oidx, dep_vis@idx),]
    }
}
# py_install("louvain", pip=TRUE, pip_ignore_installed=TRUE)
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [10PC, batch regressed]`,k=20, res = 1e-3, louvain_path = "../hcc_cello/scripts/python/louvain.py", python_home = NULL, random_seed = NULL))

saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))





test_text <- "B_cell_invivo_20221224_IFF_"
test_name <- "B cells in vivo [20221224, IFF]"
dep_vis <- clist$`B cells in vivo [20221224]`
oidx <- dep_vis@idx
cds_oidx <- eset[, oidx]
dep_clus <- dep_vis@pmeta$Louvain_Cluster[match(oidx, dep_vis@idx)]
newvis <- new("Cello", name = test_name, idx=oidx)
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster =dep_clus, cluster_min_cell_num = 10, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]

max_numpc <- 100
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text,"irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}
# Regress batch
residualModelFormulaStr <- "~Batch"
batch_code <- sapply(strsplit(cds_oidx$Dataset, "_"), function(x)x[2])
cds_oidx$Batch <- ifelse(as.numeric(batch_code) >=40, "1", "0")
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
#newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, batch regressed]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component, n_neighbors = 30L)
    }
}

# py_install("louvain", pip=TRUE, pip_ignore_installed=TRUE)
newvis@pmeta <- data.frame(Louvain_Cluster = louvain_clus(data = newvis@proj$`UMAP-2D [10PC, IFG]`,k=15, res = 5e-3, louvain_path = "../hcc_cello/scripts/python/louvain.py", python_home = NULL, random_seed = NULL))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(savePath,"clist.rds"))



dep_vis <- clist$`B cells in vivo [20221224, IFF]`
dep_clus <- as.character(dep_vis@pmeta$Louvain_Cluster)
dep_clus[dep_clus %in% c(11,16)] <- "B cell"
dep_clus[dep_clus %in% c(23,24,5,21,17)] <- "Cycling B/Plasma cell"
dep_clus[dep_clus %in% c(2,13)] <- "IGHA+ Plasma cell (Normal)"
dep_clus[dep_clus %in% c(3,8)] <- "IGHA+ Plasma cell (Tumor)"
dep_clus[dep_clus %in% c(10,18,15,7,4,12,14,9,1,20,19)] <- "IGHG+ Plasma cell"
dep_clus[dep_clus %in% c(22,6)] <- "IGHM+ Plasma cell"
pData(eset)$Cell_type_B <- "unannotated"
pData(eset)$Cell_type_B[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(savePath, "eset.rds"))


