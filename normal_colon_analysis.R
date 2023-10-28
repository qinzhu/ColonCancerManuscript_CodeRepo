


# Figure to show cell type composition in normal colon



library(VisCello)
library(igraph)
savePath <- "../hcc_final_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)
scriptPath <- paste0("../hcc_cello/", "scripts/")

mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

eset <- readRDS("../hcc_final_cello/eset.rds")
clist <- readRDS("../hcc_final_cello/clist.rds")


test_text <- "normal_human_colon_"
dep_vis <- clist$`Colon cells [IFF]`
eset_normal <- eset[, dep_vis@idx]

cur_proj <- cbind(dep_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset_normal))
cur_proj$Cluster <- dep_vis@pmeta$Cluster

cur_proj$Colon_cell_type <- NA
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(13,8,22,2,14,15,18)] = "TA"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(4)] = "SC"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(6,20,9)] = "PRO"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(1,11,3,19,5,16,12,10)] = "CC"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(21)] = "EEC"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(17)] = "G"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(7)] = "PLC"

cc_level <- c("SC", "TA", "PRO", "CC", "G", "EEC", "PLC")
cur_proj$Colon_cell_type <- factor(cur_proj$Colon_cell_type, levels = cc_level)
ctype_color <- get_factor_color(cc_level, "Set2")
#dataset_color <- c(dataset_color[5:length(dataset_color)], dataset_color[1:4])
names(ctype_color) <- cc_level

label_data <- cur_proj %>% group_by_at("Colon_cell_type") %>% summarize_at(c("UMAP_1", "UMAP_2"), median)

library(ggrastr)
#library(ggrepel)
g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Colon_cell_type),size=.5, stroke=0) +
    #geom_text(data = label_data, aes(UMAP_1, UMAP_2, label=Colon_cell_type)) + 
    scale_color_manual(values = ctype_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "colon_ctype_umap_rast.pdf"), g1, width = 3, height=1.5, units = "in", device = "pdf")


cc_color <- c("G1" = "#66bd63", "S" = "#fee08b", "G2M" = "#d73027")
cur_proj$Cell_cycle_phase <- factor(cur_proj$Cell_cycle_phase, levels = names(cc_color))
g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Cell_cycle_phase),size=.5, stroke=0) +
    #geom_text(data = label_data, aes(UMAP_1, UMAP_2, label=Dataset)) + 
    scale_color_manual(values = cc_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "colon_cellcycle_umap_rast.pdf"), g1, width = 3, height=1.5, units = "in", device = "pdf")

dataset_color <- get_factor_color(unique(cur_proj$Dataset), "Accent")
names(dataset_color) <- unique(cur_proj$Dataset)

g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Dataset),size=.5, stroke=0) +
    #geom_text(data = label_data, aes(UMAP_1, UMAP_2, label=Dataset)) + 
    scale_color_manual(values = dataset_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "colon_dataset_umap_rast.pdf"), g1, width = 3.3, height=1.5, units = "in", device = "pdf")


# Heatmap of top DEGs
set.seed(2020)
bg_count <- 1000
max_count <- 500
eset_normal$Colon_cell_type <- cur_proj$Colon_cell_type
saveRDS(eset_normal, paste0(mstatePath, test_text, "eset_normal.rds"))

eset_de <- eset_normal
eset_de$de_group <- eset_normal$Colon_cell_type
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

pdeg_list<- list()
for(ctype in cc_level){
    test_pair <- c(ctype, "Background")
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <-  which(pData(eset_de)$de_group != test_pair[1])
    group2_idx <- sample(group2_idx, min(length(group2_idx), bg_count))
    test_idx <- c(group1_idx, group2_idx)
    test_clus <- c(rep(test_pair[1], length(group1_idx)), rep(test_pair[2], length(group2_idx)))
    cur_cds <- eset_de[, test_idx]
    cur_cds$test_cluster <- factor(test_clus, levels = test_pair)
    feature_data <- fData(cur_cds)
    prioritized_genes <- runsSeq(dat=as.matrix(exprs(cur_cds)), group=cur_cds$test_cluster, fdata = feature_data, order_by="pvalue", p_cutoff= .05, min_mean = 0, min_log2fc = 0, id_col = id_col, name_col = name_col)
    prioritized_genes <- lapply(prioritized_genes, function(x) {
        x$proportion1 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_pair[1]]) > 0)
        x$proportion2 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_pair[2]]) > 0)
        return(x)
    })
    pdeg_list[[ctype]] <- prioritized_genes[[1]]
}
saveRDS(pdeg_list, paste0(mstatePath, test_text, "colon_ctype_deg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})

saveRDS(deg_list, paste0(mstatePath,test_text,"colon_ctype_deg_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"colon_ctype_deg_sig_genes.xlsx"))

# Pheatmap
eset_plot <- eset_normal
max_cell = 100
max_gene <- 5
plot_sample <- unlist(lapply(names(table(eset_plot$Colon_cell_type)), function(x) {
    cur_idx <- which(eset_normal$Colon_cell_type == x)
    sample(cur_idx, min(max_cell, length(cur_idx)), replace = F)
}))
plot_gene <- as.character(unlist(lapply(deg_list[names(deg_list) != "PRO"], function(x) {
    x$gene_id[1:min(max_gene, length(x$gene_id))]
})))
#plot_gene_name <- fData(eset_plot)$gene_short_name[match(plot_gene, fData(eset_plot)$gene_id)]
eset_plot <- eset_plot[match(plot_gene, rownames(eset_plot)), plot_sample]
#eset_plot <- eset_plot[rowMeans(exprs(eset_plot) > 0) > .2, ]
plot_meta <- pData(eset_plot)
plot_meta <- plot_meta[, c("Colon_cell_type", "Dataset"), drop=F]
plot_meta <- plot_meta[order(plot_meta[["Colon_cell_type"]], plot_meta$Dataset),]
non_na_cells_ordered <- rownames(plot_meta)
value <- eset_plot@assayData$norm_exprs[,non_na_cells_ordered]
rownames(value) <- fData(eset_plot)$gene_short_name[match(rownames(value), rownames(eset_plot))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

cc_level <- c("SC", "TA", "PRO", "EC", "G", "EEC", "PLC")
ctype_color <- get_factor_color(cc_level, "Set2")
#dataset_color <- c(dataset_color[5:length(dataset_color)], dataset_color[1:4])
names(ctype_color) <- cc_level

dataset_color <- get_factor_color(unique(eset_plot$Dataset), "Accent")
names(dataset_color) <- unique(eset_plot$Dataset)

pdf(paste0(mplotPath, test_text, "colon_ctype_deg", ".pdf"), width = 5.2, height=4.1)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = F, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_colors = list(Colon_cell_type = ctype_color, Dataset = dataset_color),
                   annotation_names_col = FALSE, 
                   fontsize = 8)
dev.off()








test_text <- "normal_mouse_colon_"
eset_mouse <- readRDS("../mouse_colon_cello/eset.rds")
clist_mouse <- readRDS("../mouse_colon_cello/clist.rds")

# Assign cell cycle
library(Seurat)
HMD_HumanPhenotype <- HMD_HumanPhenotype <- read.delim("~/Documents/CHOP/HSC/Preprocess.VisCello.eht/data-raw/public_resource/HMD_HumanPhenotype.rpt", header=FALSE)
s.genes <- as.character(HMD_HumanPhenotype$V5[match(cc.genes$s.genes, HMD_HumanPhenotype$V1)])
g2m.genes <- as.character(HMD_HumanPhenotype$V5[match(cc.genes$g2m.genes, HMD_HumanPhenotype$V1)])

mouse_expr <- exprs(eset_mouse)
rownames(mouse_expr) <- fData(eset_mouse)$gene_short_name[match(rownames(mouse_expr), rownames(eset_mouse))]
seurat.mm <- CreateSeuratObject(counts = mouse_expr, min.cells = 5)
seurat.mm <- AddMetaData(seurat.mm, metadata = pData(eset_mouse[,dep_vis@idx]))
seurat.mm <- NormalizeData(seurat.mm , verbose = T)
seurat.mm<- CellCycleScoring(seurat.mm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

eset_mouse$Cell_cycle_phase <- NA
eset_mouse$Cell_cycle_phase[match(colnames(seurat.mm), colnames(eset_mouse))] <- as.character(seurat.mm$Phase)
saveRDS(eset_mouse, "../mouse_colon_cello/eset.rds")


dep_vis <- clist_mouse$`Mouse normal colon epithelial [IFF]`

cur_proj <- cbind(dep_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset_mouse[, dep_vis@idx]))
cur_proj$Cluster <- dep_vis@pmeta$Cluster

cur_proj$Colon_cell_type <- NA
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(25,23,5,15,20,2,14,18,24)] = "TA"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(3,22,10)] = "SC"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(12,7,9)] = "PRO"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(8,21,19,6,26)] = "CC"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(28,27)] = "EEC"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(13,16,17,1,4,11)] = "G"
cur_proj$Colon_cell_type[cur_proj$Cluster %in% c(29)] = "unannotated"


cc_level <- c("SC", "TA", "PRO", "CC", "G", "EEC", "unannotated")
cur_proj$Colon_cell_type <- factor(cur_proj$Colon_cell_type, levels = cc_level)
ctype_color <- get_factor_color(cc_level, "Set2")
#dataset_color <- c(dataset_color[5:length(dataset_color)], dataset_color[1:4])
names(ctype_color) <- cc_level
cur_proj <- cur_proj[cur_proj$Colon_cell_type != "unannotated",]
label_data <- cur_proj %>% group_by_at("Colon_cell_type") %>% summarize_at(c("UMAP_1", "UMAP_2"), median)
#library(ggrepel)
g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Colon_cell_type),size=.5, stroke=0) +
    #geom_text(data = label_data, aes(UMAP_1, UMAP_2, label=Colon_cell_type)) + 
    scale_color_manual(values = ctype_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "colon_ctype_umap_rast.pdf"), g1, width = 3.5, height=1.4, units = "in", device = "pdf")

g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Cell_cycle_phase),size=.5, stroke=0) +
    #geom_text(data = label_data, aes(UMAP_1, UMAP_2, label=Colon_cell_type)) + 
    scale_color_manual(values = cc_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "colon_cellcycle_umap_rast.pdf"), g1, width = 3.1, height=1.5, units = "in", device = "pdf")


dataset_color <- get_factor_color(unique(cur_proj$Dataset), "Accent")
names(dataset_color) <- unique(cur_proj$Dataset)

g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Dataset),size=.5, stroke=0) +
    #geom_text(data = label_data, aes(UMAP_1, UMAP_2, label=Dataset)) + 
    scale_color_manual(values = dataset_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "colon_dataset_umap_rast.pdf"), g1, width = 3.3, height=1.5, units = "in", device = "pdf")

# Heatmap of top DEGs
set.seed(2020)
bg_count <- 1000
max_count <- 500
eset_normal <- eset_mouse[, rownames(cur_proj)]
eset_normal$Colon_cell_type <- cur_proj$Colon_cell_type
saveRDS(eset_normal, paste0(mstatePath, test_text, "eset_normal.rds"))

eset_de <- eset_normal
eset_de$de_group <- eset_normal$Colon_cell_type
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

pdeg_list<- list()
for(ctype in cc_level){
    test_pair <- c(ctype, "Background")
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <-  which(pData(eset_de)$de_group != test_pair[1])
    group2_idx <- sample(group2_idx, min(length(group2_idx), bg_count))
    test_idx <- c(group1_idx, group2_idx)
    test_clus <- c(rep(test_pair[1], length(group1_idx)), rep(test_pair[2], length(group2_idx)))
    cur_cds <- eset_de[, test_idx]
    cur_cds$test_cluster <- factor(test_clus, levels = test_pair)
    feature_data <- fData(cur_cds)
    prioritized_genes <- runsSeq(dat=as.matrix(exprs(cur_cds)), group=cur_cds$test_cluster, fdata = feature_data, order_by="pvalue", p_cutoff= .05, min_mean = 0, min_log2fc = 0, id_col = id_col, name_col = name_col)
    prioritized_genes <- lapply(prioritized_genes, function(x) {
        x$proportion1 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_pair[1]]) > 0)
        x$proportion2 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_pair[2]]) > 0)
        return(x)
    })
    pdeg_list[[ctype]] <- prioritized_genes[[1]]
}
saveRDS(pdeg_list, paste0(mstatePath, test_text, "colon_ctype_deg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})

saveRDS(deg_list, paste0(mstatePath,test_text,"colon_ctype_deg_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"colon_ctype_deg_sig_genes.xlsx"))

# Pheatmap
eset_plot <- eset_normal
max_cell = 100
max_gene <- 5
plot_sample <- unlist(lapply(names(table(eset_plot$Colon_cell_type)), function(x) {
    cur_idx <- which(eset_normal$Colon_cell_type == x)
    sample(cur_idx, min(max_cell, length(cur_idx)), replace = F)
}))
low_expressed_gene <- which(rowMeans(exprs(eset_plot)>0)<0.025)
plot_gene <- as.character(unlist(lapply(deg_list[names(deg_list) != "PRO"], function(x) {
    use_id <- x$gene_id[!x$gene_id %in% names(low_expressed_gene)]
    use_id[1:min(max_gene, length(use_id))]
})))
plot_gene <- plot_gene[!is.na(plot_gene)]
#plot_gene_name <- fData(eset_plot)$gene_short_name[match(plot_gene, fData(eset_plot)$gene_id)]
eset_plot <- eset_plot[match(plot_gene, rownames(eset_plot)), plot_sample]
#eset_plot <- eset_plot[rowMeans(exprs(eset_plot) > 0) > .2, ]
plot_meta <- pData(eset_plot)
plot_meta <- plot_meta[, c("Colon_cell_type", "Dataset"), drop=F]
plot_meta <- plot_meta[order(plot_meta[["Colon_cell_type"]], plot_meta$Dataset),]
non_na_cells_ordered <- rownames(plot_meta)
value <- eset_plot@assayData$norm_exprs[,non_na_cells_ordered]
rownames(value) <- fData(eset_plot)$gene_short_name[match(rownames(value), rownames(eset_plot))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

cc_level <- c("SC", "TA", "PRO", "CC", "G", "EEC", "unannotated")
ctype_color <- get_factor_color(cc_level, "Set2")
names(ctype_color) <- cc_level
ctype_color <- ctype_color[names(ctype_color) != "unannotated"]
#dataset_color <- c(dataset_color[5:length(dataset_color)], dataset_color[1:4])

dataset_color <- get_factor_color(unique(eset_plot$Dataset), "Accent")
names(dataset_color) <- unique(eset_plot$Dataset)
plot_meta$Colon_cell_type <- factor(as.character(plot_meta$Colon_cell_type), levels = c("SC", "TA", "PRO", "CC", "G", "EEC"))
pdf(paste0(mplotPath, test_text, "colon_ctype_deg", ".pdf"), width = 5.2, height=4.1)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = F, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_colors = list(Colon_cell_type = ctype_color, Dataset = dataset_color),
                   annotation_names_col = FALSE, 
                   fontsize = 8)
dev.off()



