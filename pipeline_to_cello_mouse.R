



library(cellrangerRkit)
library(stringi)
library(monocle)
library(dplyr)
library(openxlsx)
library(reshape2)
library(R.utils)
library(VisCello)
library(igraph)

mousePath <- "../mouse_colon_cello/";dir.create(mousePath)
scriptPath <- paste0("../hcc_cello//scripts/")
#sourceDirectory(scriptPath)
source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))
source(paste0(scriptPath, "compute_go.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "Load_10x_latest.R"))
set.seed(2020)
exprs <- Biobase::exprs

rdsPath <- paste0(mousePath, "rds/"); dir.create(rdsPath)
plotPath <- paste0(mousePath,"plots/"); dir.create(plotPath)
sheetPath <- paste0(mousePath,"sheets/"); dir.create(sheetPath)

#dataPath <- "/kimdata/zhuqin/Ningproj/"
#dataPath <- "/Volumes/MPalace/LAB/KLBackup/zhuqin/Ningproj/"
sample_dirs <- c(
                 "Colon1" = paste0("~/Documents/CHOP/NingProj/Ningcolon/Control_td/outs/raw_gene_bc_matrices/mm10_tomato/"),
                 "Colon2" = paste0("~/Documents/CHOP/NingProj/CLNL10xSIColAOM/Colon2/outs/raw_gene_bc_matrices/mm10_tomato/"),
                 "Tumor ApcMin" =paste0("~/Documents/CHOP/NingProj/Ningcolon/ApcMin_td/outs/raw_gene_bc_matrices/mm10_tomato/"),
                 "Tumor AOM" = paste0("~/Documents/CHOP/NingProj/CLNL10xSIColAOM/TumorAOM/outs/raw_gene_bc_matrices/mm10_tomato/"),
                 "Tumor APKS" = paste0("~/Documents/CHOP/NingProj/APKS_outs/raw_feature_bc_matrix/"),
                 "ApcFlox" = paste0("~/Documents/CHOP/NingProj/merge190522_0601/ApcMerge/outs/raw_feature_bc_matrix/"),
                 "BcatFlox" = paste0("~/Documents/CHOP/NingProj/merge190522_0601/BcatMerge/outs/raw_feature_bc_matrix/")
)

cbn_cds <- combine_monocle_datasets(sample_dirs,minUMI = 1000, minGenes = 1000)


dim(cbn_cds)
cbn_cds <- estimateSizeFactors(cbn_cds)
cbn_cds <- detectGenes(cbn_cds, min_expr = 0.1)

pData(cbn_cds)$Total_mRNAs <- Matrix::colSums(exprs(cbn_cds))

sample_order <- names(sample_dirs)
sorder <- stringi::stri_extract_last(colnames(cbn_cds), regex = "\\d+")
sdset <- factor(sample_order[as.numeric(sorder)], levels= names(sample_dirs))
pData(cbn_cds)$Dataset <- sdset

cbn_cds@assayData$exprs<-as(cbn_cds@assayData$exprs, 'sparseMatrix')
FM <- normalize_expr_data2(cbn_cds, "log", 1, use_order_gene = F)
cbn_cds@assayData$norm_exprs <- Matrix(FM, sparse=T)

mito_genes <- grep("^mt-", fData(cbn_cds)$gene_short_name, value = F)
mito_values <- as.matrix(exprs(cbn_cds[mito_genes,]))
mito_expr_prop <- sweep(mito_values, 2, Matrix::colSums(exprs(cbn_cds)), "/")
mito_expr_prop_sum <- colSums(mito_expr_prop)
pdf(paste0(plotPath, "all_data_mitofrac.pdf"))
hist(mito_expr_prop_sum, breaks=100)
dev.off()
pData(cbn_cds)$MitoFraction <- mito_expr_prop_sum
#saveRDS(cbn_cds, paste0(rdsPath, "cbn_cds_200416.rds"))
saveRDS(cbn_cds, paste0(rdsPath, "cbn_cds_200430.rds"))

#cbn_cds <- readRDS(paste0(rdsPath, "cbn_cds_200416.rds"))
cbn_cds <- readRDS(paste0(rdsPath, "cbn_cds_200430.rds"))
fmeta <- fData(cbn_cds)
eset <- new("ExpressionSet",
            assayData = assayDataNew( "environment", exprs=exprs(cbn_cds), norm_exprs = cbn_cds@assayData$norm_exprs),
            phenoData =  new("AnnotatedDataFrame", data = pData(cbn_cds)),
            featureData = new("AnnotatedDataFrame", data = fmeta))
saveRDS(eset, paste0(mousePath, "eset.rds"))


eset <- readRDS(paste0(mousePath, "eset.rds"))
clist <- readRDS(paste0(mousePath, "clist.rds"))
#clist <- list()

################# Mouse all dataset combine #################
test_text <- "mouse_cbn_"
test_name <- "Mouse all data"
oidx <- 1:ncol(cbn_cds)
newvis <- new("Cello", name = test_name, idx=oidx)
cds_oidx <- filter_cds(cds=cbn_cds[, oidx], min_detect=1, min_numc_expressed = 10, min_disp_ratio=1)
max_numpc <- 100
irlba_res <- compute_pca_cds(cds_oidx, num_dim = max_numpc, scvis=NULL, use_order_gene = T, residualModelFormulaStr = NULL, return_type="irlba")
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2,3)){
    for(n_pc in c(20,50)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, VEG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(newvis@proj[[2]],k=30, res = 1e-3, louvain_path = paste0(scriptPath, "python/louvain.py")))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))


################# Mouse all dataset combine, IFF #################
test_text <- "mouse_cbn_iff_"
test_name <- "Mouse all data [IFF, low mito]"
oidx <- which(eset$MitoFraction <= .2)
dep_vis <- clist$`Mouse all data`
dep_clus <- dep_vis@pmeta$Cluster
use_ifg <- ifg_select(data = exprs(eset), cluster = dep_clus, cluster_min_cell_num = 100, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 

newvis <- new("Cello", name = test_name, idx=oidx)
cds_oidx <- eset[rownames(fData(eset)) %in% use_ifg, oidx]
max_numpc <- 50
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,20,50)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
#source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[3]],k=30, res = 1e-5, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))

# Assign cell types
dep_vis <- clist$`Mouse all data [IFF, low mito]`
dep_clus <- as.character(dep_vis@pmeta$Cluster)
dep_clus[dep_clus %in% c(2)] <- "Goblet cells"
dep_clus[dep_clus %in% c(1)] <- "Epithelial cells [EC,TA,SC]"
dep_clus[dep_clus %in% c(7,8)] <- "Myeloid"
dep_clus[dep_clus %in% c(16)] <- "Endothelial"
dep_clus[dep_clus %in% c(15)] <- "Fibroblast"
dep_clus[dep_clus %in% c(14)] <- "B cell"
dep_clus[dep_clus %in% c(12,10)] <- "Plasma cells"
dep_clus[dep_clus %in% c(9)] <- "T cells"
dep_clus[dep_clus %in% c(13)] <- "EEC cells"
dep_clus[dep_clus %in% c(5,6)] <- "Epithelial cells (tumor)"
dep_clus[dep_clus %in% c(11)] <- "EEC cells (organoid)"
dep_clus[dep_clus %in% c(3,4)] <- "Epithelial cells (organoid)"

pData(eset)$Cell_type <- NA
pData(eset)$Cell_type[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(mousePath, "eset.rds"))


################# Mouse all in vivo, IFF #################
test_text <- "mouse_vivo_iff_"
test_name <- "Mouse in vivo [IFF, low mito]"
dep_vis <- clist$`Mouse all data [IFF, low mito]`
filter_std <- !eset$Dataset %in% c("ApcFlox", "BcatFlox") & 1:ncol(eset) %in% dep_vis@idx
oidx <- which(filter_std)
cds_oidx <- eset[,oidx]

use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = cds_oidx$Cell_type, cluster_min_cell_num = 30, min_cluster_expr_fraction = .1, gini_cut_qt = .75, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 

newvis <- new("Cello", name = test_name, idx=oidx)
cds_oidx <- eset[rownames(fData(eset)) %in% use_ifg, oidx]
max_numpc <- 50
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
#source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(data = newvis@proj[[2]],k=30, res = 1e-5, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))


# Refine cell types
dep_vis <- clist$`Mouse in vivo [IFF, low mito]`
dep_clus <- as.character(dep_vis@pmeta$Cluster)
dep_clus[dep_clus %in% c(2)] <- "Goblet cells"
dep_clus[dep_clus %in% c(1)] <- "Epithelial cells [EC,TA,SC]"
dep_clus[dep_clus %in% c(3)] <- "Myeloid"
dep_clus[dep_clus %in% c(10)] <- "Endothelial"
dep_clus[dep_clus %in% c(9)] <- "Fibroblast"
dep_clus[dep_clus %in% c(6)] <- "B & Plasma cell"
dep_clus[dep_clus %in% c(7)] <- "T cells"
dep_clus[dep_clus %in% c(8)] <- "EEC cells"
dep_clus[dep_clus %in% c(4,5)] <- "Epithelial cells (tumor)"

pData(eset)$Cell_type_refined <- NA
pData(eset)$Cell_type_refined[dep_vis@idx] <- dep_clus

saveRDS(eset, paste0(mousePath, "eset.rds"))



################# Mouse all dataset combine, IFF #################
test_text <- "mouse_myeloid_"
test_name <- "Mouse myeloid"
oidx <- which(eset$Cell_type_refined == "Myeloid")
dep_vis <- clist$`Mouse in vivo [IFF, low mito]`
dep_clus <- dep_vis@pmeta$Cluster
cds_oidx <- eset[, oidx]
expressed_g <- rowMeans(exprs(cds_oidx) > 0) > .1
cds_oidx <- cds_oidx[expressed_g,]

newvis <- new("Cello", name = test_name, idx=oidx)
max_numpc <- 50
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
newvis@proj[["UMAP-2D Zoom"]] <- dep_vis@proj$`UMAP-2D [10PC, IFG]`[match(rownames(pca_proj), rownames(dep_vis@proj$`UMAP-2D [10PC, IFG]`)), ]
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
library(phateR)
proj_PHATE <- phate(t(cds_oidx@assayData$norm_exprs))
newvis@proj[["PHATE-2D"]] <- proj_PHATE$embedding

#source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(newvis@proj$`UMAP-2D Zoom`,k=30, res = 5e-3, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))

newvis <- readRDS(paste0(rdsPath, test_text, "cello.rds"))
de_list <- read_excel_allsheets(paste0("~/Dropbox/ColonManuscript/sheets/", "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))
max_gnum = 100
sig_deg <- lapply(de_list, function(x) {
    x$gene_short_name[1:min(max_gnum, nrow(x))]
})
sig_deg_mouse <- lapply(sig_deg, function(x) {
    y <- mouse_to_human_symbol(gene_symbol = x, in.type = "hs",HMD_HumanPhenotype = HMD_HumanPhenotype)
    return(y[!is.na(y)])
})
cds_oidx <- eset[,oidx]
expressed_g <- rowMeans(exprs(cds_oidx) > 0) > .01
cds_oidx <- cds_oidx[expressed_g,]
cds_oidx <- cds_oidx[fData(cds_oidx)$gene_short_name %in% unlist(sig_deg_mouse), ]
max_numpc <- 10
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
for(n_component in c(2)){
    for(n_pc in c(5,10)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, DEG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}

clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))




################# Normal colon #################
test_text <- "mouse_colon_"
test_name <- "Mouse normal colon"
oidx <- which(eset$Dataset %in% c("Colon1", "Colon2") & eset$MitoFraction <= .1)
cds_oidx <- eset[, oidx]
expressed_g <- rowMeans(exprs(cds_oidx) > 0) > .1
cds_oidx <- cds_oidx[expressed_g,]

newvis <- new("Cello", name = test_name, idx=oidx)
max_numpc <- 50
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
library(phateR)
proj_PHATE <- phate(t(cds_oidx@assayData$norm_exprs))
newvis@proj[["PHATE-2D"]] <- proj_PHATE$embedding
#source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(newvis@proj$`UMAP-2D [10PC, IFG]`,k=30, res = 5e-3, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))



################# Normal colon, IFF #################
test_text <- "mouse_colon_epithelial_iff_"
test_name <- "Mouse normal colon epithelial [IFF]"
dep_vis <- clist$`Mouse normal colon`
oidx <- which(1:ncol(eset) %in% dep_vis@idx & eset$Cell_type %in% c("EEC cells", "EEC cells (organoid)", "Epithelial cells (tumor)", "Epithelial cells [EC,TA,SC]", "Goblet cells"))
dep_clus <- dep_vis@pmeta$Cluster
cds_oidx <- eset[, oidx]
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = dep_clus[match(oidx, dep_vis@idx)], cluster_min_cell_num = 30, min_cluster_expr_fraction = .1, gini_cut_qt = .8, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]

residualModelFormulaStr <- "~Dataset"
FM <- cds_oidx@assayData$norm_exprs
if (!is.null(residualModelFormulaStr)) {
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds_oidx), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
}

newvis <- new("Cello", name = test_name, idx=oidx)
max_numpc <- 50
irlba_res <- prcomp_irlba(t(FM), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
# library(phateR)
# proj_PHATE <- phate(t(cds_oidx@assayData$norm_exprs))
# newvis@proj[["PHATE-2D"]] <- proj_PHATE$embedding
#source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(newvis@proj$`UMAP-2D [10PC, IFG]`,k=15, res = 1e-3, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))



################# Epithelial, IFF #################
test_text <- "mouse_epithelial_iff_"
test_name <- "Mouse epithelial [IFF]"
dep_vis <- clist$`Mouse in vivo [IFF, low mito]`
oidx <- which(eset$Cell_type %in% c("EEC cells", "EEC cells (organoid)", "Epithelial cells (tumor)", "Epithelial cells [EC,TA,SC]", "Goblet cells") & 1:ncol(eset) %in% dep_vis@idx)
#dep_clus <- dep_vis@pmeta$Cluster[match(oidx, dep_vis@idx)]
cds_oidx <- eset[, oidx]
use_ifg <- ifg_select(data = exprs(cds_oidx), cluster = cds_oidx$Dataset, cluster_min_cell_num = 30, min_cluster_expr_fraction = .1, gini_cut_qt = .8, compute_go = F, filePath = plotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
cds_oidx <- cds_oidx[use_ifg,]
newvis <- new("Cello", name = test_name, idx=oidx)
max_numpc <- 50
irlba_res <- prcomp_irlba(t(cds_oidx@assayData$norm_exprs), n = max_numpc, center = TRUE, scale. = TRUE)
saveRDS(irlba_res, paste0(rdsPath, test_text, "irlba_res.rds"))
pca_proj <- as.data.frame(irlba_res$x)
rownames(pca_proj) <- colnames(cds_oidx)
newvis@proj[["PCA"]] <- pca_proj
for(n_component in c(2)){
    for(n_pc in c(10,30)) {
        newvis@proj[[paste0("UMAP-", n_component, "D [", n_pc, "PC, IFG]")]]<-compute_umap_pca(pca_proj,use_dim = n_pc, n_component=n_component)
    }
}
# library(phateR)
# proj_PHATE <- phate(t(cds_oidx@assayData$norm_exprs))
# newvis@proj[["PHATE-2D"]] <- proj_PHATE$embedding
#source(paste0(scriptPath, "compute_dimR.R"))
newvis@pmeta <- data.frame(Cluster = louvain_clus(newvis@proj$`UMAP-2D [10PC, IFG]`,k=15, res = 1e-3, louvain_path = "~/Documents/CHOP/NingProj/hcc_cello/scripts/python/louvain.py", python_home = "~/anaconda3/envs/r-reticulate/bin/python"))
saveRDS(newvis, paste0(rdsPath, test_text, "cello.rds"))
clist[[test_name]] <- newvis
saveRDS(clist, paste0(mousePath, "clist.rds"))

