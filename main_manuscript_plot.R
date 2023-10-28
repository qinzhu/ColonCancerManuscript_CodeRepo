

library(VisCello)

is_windows = FALSE
if(is_windows) {
    home_path = "C:/Users/qinzh/"
} else {
    home_path = "~"
}
mplotPath <- paste0(home_path, "/Dropbox/ColonManuscript/subplots/")
msheetPath <- paste0(home_path, "/Dropbox/ColonManuscript/sheets/")
mstatePath<- paste0(home_path, "/Dropbox/ColonManuscript/states/")


scriptPath <- paste0("../hcc_cello/scripts/")
source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")




# Figure 1 stats
test_text <- "fig_data_stats_"
cur_vis <- clist$`All cells [low mito]`
cur_proj <- pData(eset)[cur_vis@idx,]

cur_proj <- cur_proj[!cur_proj$Patient %in% c("40", "96") & !is.na(cur_proj$Dataset),]
use_levels <- names(table(cur_proj$Dataset))
length(use_levels)
# use_levels <- c(grep("Colon_", use_levels, value = T), grep("Tumor_", use_levels, value = T), 
#                 grep("Liver_", use_levels, value = T), grep("LiverMet_", use_levels, value = T),
#                 grep("Organoid_", use_levels, value = T), grep("Tumoroid_", use_levels, value = T),
#                 grep("Macrophage_", use_levels, value = T))
cur_proj$Dset2 <- cur_proj$Dataset
cur_proj$Dset2[cur_proj$Dset2 == 'Organoid+Macrophage_86_36hr(2)'] = 'Organoid+Macrophage_86_36hr'
cur_proj$Dset2[cur_proj$Dset2 == 'Macrophage_24hr'] = 'Macrophage_24hr(1)'
cur_proj$Dset2[cur_proj$Dset2 == 'Macrophage_b1_24hr'] = 'Macrophage_24hr(2)'
use_levels <- c('Colon_07', 'Colon_08', 'Colon_09', 
                'Tumor_07', 'Tumor_08', 'Tumor_09', 'Tumor_24', 'Tumor_26', 'Tumor_27', 'Tumor_28', 
                'Tumor_44', 'Tumor_49', 'Tumor_80', 'Tumor_81', 'Tumor_83', 'Tumor_84', 'Tumor_86', 'Tumor_87', 
                'Liver_09', 'Liver_44', 'LiverMet_09', 'LiverMet_27', 'LiverMet_44', 'LiverMet_49', 
                'Organoid_08', 'Organoid_24', 'Organoid_26', 'Organoid_28', 
                'Tumoroid_08', 'Tumoroid_24', 'Tumoroid_28', 
                'Tumoroid_28_24hr', 'Tumoroid_28_48hr', 
                'Tumoroid_86_36hr', 
                'Organoid_86_36hr', 'Organoid_92_24hr', 'Tumoroid_92_24hr', 
                'Macrophage_24hr(1)', 'Macrophage_36hr', 'Macrophage_48hr', 'Macrophage_24hr(2)', 
                'Organoid+Macrophage_86_36hr', 'Organoid+Macrophage_92_24hr', 
                'Tumoroid+Macrophage_28_24hr', 'Tumoroid+Macrophage_28_48hr', 'Tumoroid+Macrophage_86_36hr', 'Tumoroid+Macrophage_92_24hr')
any(!names(table(cur_proj$Dset2)) %in% use_levels); any(!use_levels %in% names(table(cur_proj$Dset2)))

cur_proj$Dset2 <- factor(cur_proj$Dset2, levels = use_levels)
df_med <- cur_proj %>% group_by(Dset2) %>% summarize(med_g = median(num_genes_expressed), med_umi = median(Total_mRNAs))
g1<-ggplot(cur_proj, aes(x=Dset2, y=num_genes_expressed,fill=Dset2)) + 
    geom_boxplot(lwd=0.3, outlier.size = 0.1, show.legend = F) + 
    geom_point(data = df_med, aes(x=Dset2, y=med_g), size = 0, stroke=0, shape = 22) +
    #scale_fill_manual(values = combine_dataset_color) +
    scale_y_continuous(limits = c(0, 8000))+
    labs(x="", y = "#Genes/cell")+
    guides(fill=F)+
    theme_bw() + 
    theme(text = element_text(size=6), 
          legend.text=element_text(size=6),
          axis.text = element_text(size=6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "right", 
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "dataset_numg.eps"), g1, width = 6, height=3, units = "in", device = "eps")

write.csv(cur_proj[,c("Dset2", "num_genes_expressed")], paste0(msheetPath, test_text, "dataset_numg.csv"))

g1<-ggplot(cur_proj, aes(x=Dset2, y=Total_mRNAs,fill=Dset2)) + 
    geom_boxplot(lwd=0.3, outlier.size = 0.1, show.legend = F) + 
    #geom_point(data = df_med, aes(x=Dset2, y=med_umi), size = 0, stroke=0, shape = 22) +
    #scale_fill_manual(values = combine_dataset_color) +
    scale_y_log10(limits = c(1e3,100000))+
    labs(x="", y = "#UMIs/cell")+
    guides(fill=F)+
    theme_bw() + 
    theme(text = element_text(size=6), 
          legend.text=element_text(size=6),
          axis.text = element_text(size=6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "right", 
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "dataset_numumi.eps"), g1, width = 6, height=3, units = "in", device = "eps")

write.csv(cur_proj[,c("Dset2", "Total_mRNAs")], paste0(msheetPath, test_text, "dataset_numumi.csv"))

plot_proj <- as.data.frame(table(cur_proj$Dset2))
colnames(plot_proj) <- c("Dataset", "Cell_count")
g1 <- ggplot(plot_proj, aes(x = Dataset, y = Cell_count, fill = Dataset)) +
    geom_bar(position="dodge", stat="identity") +
    guides(fill=F)+
    theme_bw() + 
    theme(text = element_text(size=6), 
          legend.text=element_text(size=6),
          axis.text = element_text(size=6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "right", 
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "dataset_cell_count.eps"), g1, width = 6, height=3, units = "in", device = "eps")


cur_proj$Cell_type_broad <- ifelse(grepl("Epithelial", cur_proj$Cell_type), "Epithelial", "Non-epithelial")
plot_proj <- cur_proj
plot_proj$Dataset <- as.character(plot_proj$Dataset)
plot_proj$Dataset[plot_proj$Dataset == "Tumoroid_86_36hr"] = "Tumoroid_86"
plot_proj$Dataset[plot_proj$Dataset == "Organoid_86_36hr"] = "Organoid_86"
plot_proj$Dataset[plot_proj$Dataset == "Tumoroid_92_24hr"] = "Tumoroid_92"
plot_proj$Dataset[plot_proj$Dataset == "Organoid_92_24hr"] = "Organoid_92"
plot_proj <- plot_proj[!grepl("hr", plot_proj$Dataset),]
use_levels <- names(table(plot_proj$Dataset))
use_levels <- use_levels[c(grep("Colon_", use_levels), 
                grep("Tumor_", use_levels),
                grep("Liver_", use_levels), 
                grep("LiverMet_", use_levels),
                grep("Organoid_", use_levels),
                grep("Tumoroid_", use_levels))]
                
plot_proj$Dataset <- factor(plot_proj$Dataset, levels = use_levels)
plot_proj$Cell_type_broad <- factor(plot_proj$Cell_type_broad, levels = c("Non-epithelial", "Epithelial"))
#plot_proj <- as.data.frame(table(plot_proj$Dataset, plot_proj$Cell_type_broad))
#colnames(plot_proj) <- c("Dataset", "Cell_type", "Cell_count")
g1 <- ggplot(plot_proj, aes(x = Dataset, fill = Cell_type_broad)) +
    geom_bar() +
    scale_fill_manual(values = c("Epithelial" = "#1f78b4", "Non-epithelial" = "#fb9a99"))+
    theme_bw() + 
    xlab("") + 
    ylab("Cell count") + 
    theme(text = element_text(size=6), 
          legend.text=element_text(size=7),
          axis.text = element_text(size=7),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "top", 
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "dataset_cell_compos.pdf"), g1, width = 4, height=2.5, units = "in", device = "pdf")

write.csv(plot_proj[,c("Dataset", "Cell_type_broad")], paste0(msheetPath, test_text, "dataset_cell_compos.csv"))



# Figure for invivo global UMAP
# test_text <- "fig_invivo_0606_"
# test_text <- "fig_invivo_0816_"
test_text <- "fig_invivo_221225_"
cur_vis <- clist$`All in vivo cells [20210411, IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
cur_proj$Louvain_Cluster = cur_vis@pmeta$Cluster
filter_std = !(cur_proj$Patient %in% c("40") | cur_proj$Louvain_Cluster %in% c(13, 28) | 
    (cur_proj$Louvain_Cluster == "23" & cur_proj$Cell_type != "Mast cell"))
cur_proj <- cur_proj[filter_std,]

ctype_color <- c(
    "Epithelial(Normal)" = "#1f78b4",
    "Epithelial(Tumor)" = "#e31a1c",
    "Fibroblast" = "#cab2d6", 
    "Myofibroblast" = "#fb9a99",
    "Endothelial"= "#33a02c",
    "Myeloid" = "#6a3d9a",
    "T cell" ="#a6cee3",
    "B cell" = "#fdbf6f",
    "Plasma cell" = "#ff7f00",
    "Mast cell" = "#ae017e",
    "TROP2+ Liver Progenitor" = "#b2df8a",
    "Macrophage-T cell doublet" = "#b15928"
)
library(ggrastr)
cur_proj$Cell_type <- factor(cur_proj$Cell_type, levels = names(ctype_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_type", pal=ctype_color, size = .3, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "cell_type.pdf"), g1, width = 5, height=2.8, units = "in", device = "pdf")


write.csv(cur_proj[,c("UMAP_1", "UMAP_2", "Cell_type")], paste0(msheetPath, test_text, "cell_type.csv"))


stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)
cur_proj$Sample_Type <- factor(cur_proj$SampleType, levels = names(stype_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Sample_Type", size = .3, na.col = "lightgrey", pal = stype_color, legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "sampletype.pdf"), g1, width =4, height=2.8, units = "in", device = "pdf")


write.csv(cur_proj[,c("UMAP_1", "UMAP_2", "SampleType")], paste0(msheetPath, test_text, "sampletype.csv"))


# Final revision comment, small clusters in the umap, sample composition
cur_proj$Cluster_clean <- factor(as.numeric(factor(as.character(cur_proj$Louvain_Cluster))))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cluster_clean", size = .3, na.col = "lightgrey", pal = "Set1", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "louvain_cluster.pdf"), g1, width =4, height=2.8, units = "in", device = "pdf")

write.csv(cur_proj[,c("UMAP_1", "UMAP_2", "Cluster_clean")], paste0(msheetPath, test_text, "louvain_cluster.csv"))

g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cluster_clean", size = .3, na.col = "lightgrey", pal = "Set1", legend=F,onplotAnnot = T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 

small_clusters <- names(which(table(cur_proj$Cluster_clean) < 1e3))
cluster_count <- as.data.frame.matrix(table(cur_proj$Patient, cur_proj$Cluster_clean))
cluster_count_sm <- cluster_count[,small_clusters]
cluster_frac_sm <- t(t(cluster_count_sm) / colSums(cluster_count_sm))
pdf(paste0(mplotPath, test_text, "cluster_count_sm", ".pdf"), width = 3.32, height = 3)
pheatmap(cluster_frac_sm, cluster_rows = F, cluster_cols = F, display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, color = colorRampPalette(c("white", "#c73647"))(10))
dev.off()

write.csv(cluster_count, paste0(msheetPath, test_text, "cluster_count.csv"))


# DE to show cancer epithelial cell state specific gene expression

set.seed(2020)
test_text <- "rare_epi_cluster_de_"
test_clusters <- c(10,19,20)
use_cell <- rownames(cur_proj)[cur_proj$Cluster_clean %in% test_clusters]
eset_de <- eset[, use_cell]
eset_de$Cluster_clean <- cur_proj$Cluster_clean[cur_proj$Cluster_clean %in% test_clusters]

bg_count <- 2000
max_count <- 500

eset_de$de_group <- eset_de$Cluster_clean
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

pdeg_list<- list()
for(test_c in test_clusters){
    print(test_c)
    group1_idx <- which(pData(eset_de)$de_group == test_c)
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <- which(pData(eset_de)$de_group != test_c)
    group2_idx <- sample(group2_idx, min(length(group2_idx), bg_count))
    test_idx <- c(group1_idx, group2_idx)
    test_clus <- c(rep(test_c, length(group1_idx)), rep("background", length(group2_idx)))
    cur_cds <- eset_de[, test_idx]
    cur_cds$test_cluster <- factor(test_clus, levels = c(test_c, "background"))
    feature_data <- fData(cur_cds)
    prioritized_genes <- runsSeq(dat=as.matrix(exprs(cur_cds)), group=cur_cds$test_cluster, fdata = feature_data, order_by="pvalue", p_cutoff= .05, min_mean = 0, min_log2fc = 0, id_col = id_col, name_col = name_col)
    
    prioritized_genes <- lapply(prioritized_genes, function(x) {
        x$proportion1 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_c]) > 0)
        x$proportion2 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == "background"]) > 0)
        return(x)
    })
    pdeg_list[[as.character(test_c)]] <- prioritized_genes[[1]]
}

saveRDS(pdeg_list, paste0(mstatePath, test_text, "pdeg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))

# Heatmap of expression of DEGs for rare clusters
eset_plot <- eset_de
plot_meta <- pData(eset_de)[, c("Cluster_clean", "Patient", "SampleType"), drop=F]
plot_meta <- plot_meta
plot_meta <- plot_meta[order(plot_meta[["Cluster_clean"]], plot_meta$Patient),]
non_na_cells_ordered <- rownames(plot_meta)
max_g = 10
plotg <- lapply(deg_list, function(x) x[1:max_g,c("gene_id", "gene_name")])
plotg <- do.call(rbind, plotg)
value <- eset_plot@assayData$norm_exprs[plotg$gene_id,non_na_cells_ordered]
rownames(value) <- fData(eset_plot)$gene_short_name[match(rownames(value), rownames(eset_plot))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

pdf(paste0(mplotPath, test_text, "deg_expr", ".pdf"), width = 4, height=3)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = F, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 6)
dev.off()


png(paste0(mplotPath, test_text, "deg_expr", ".png"), width = 4*300, height=3*300, res = 300)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = F, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 6)
dev.off()


















show_ctype <- c(
    "Epithelial(Normal)",
    "Epithelial(Tumor)" ,
    "Fibroblast" , 
    "Myofibroblast" ,
    "Endothelial",
    "Myeloid",
    "T cell" ,
    "B cell",
    "Plasma cell",
    "Mast cell"
)

# Wrong heatmap
# compos_frac <- compos / rowSums(compos)
# pdf(paste0(mplotPath, test_text, "compos_frac_heatmap", ".pdf"), width = 4, height = 1.88)
# pheatmap(compos_frac[c(1,2),show_ctype]*100, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.0f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "red"))(10))
# dev.off()

compos<-as.data.frame.matrix(table(cur_proj$Sample_Type, cur_proj$Cell_type))
breaksList <- seq(0,1,by=.1)

compos_rel <- compos / (compos$`Epithelial(Normal)` + compos$`Epithelial(Tumor)`)
pdf(paste0(mplotPath, test_text, "compos_rel_heatmap_1", ".pdf"), width = 3.32, height = 1.42)
pheatmap(compos_rel[c(1,2),show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)))
dev.off()

pdf(paste0(mplotPath, test_text, "compos_rel_heatmap_2", ".pdf"), width = 3.32, height = 1.42)
pheatmap(compos_rel[c(1,2),show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, color = colorRampPalette(c("white", "#3690c0"))(length(breaksList)))
dev.off()



compos<-as.data.frame.matrix(table(cur_proj$Dataset, cur_proj$Cell_type))
#compos <- compos[!grepl("Liver", rownames(compos)) & !rownames(compos) %in% c("Tumor_44", "Tumor_49", "Tumor_81", "Tumor_84", "Tumor_87"),]
compos <- compos[!grepl("Liver", rownames(compos)),]
compos_rel <- compos / (compos$`Epithelial(Normal)` + compos$`Epithelial(Tumor)`)

pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap", ".pdf"), width = 3.5, height = 3.5)
pheatmap(compos_rel[,show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)
dev.off()

pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap2", ".pdf"), width = 3.5, height = 3.5)
pheatmap(compos_rel[,show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, color = colorRampPalette(c("white", "#3690c0"))(length(breaksList)), breaks = breaksList)
dev.off()

# Sort and add color bar for Tumor stage
patient_meta <- read_excel("~/Dropbox/ColonManuscript/patient_meta.xlsx")
patient_df = as.data.frame(patient_meta[match(readr::parse_number(rownames(compos_rel)), patient_meta$Patient),])
rownames(patient_df) = rownames(compos_rel)
patient_df$Sample_type = sapply(strsplit(rownames(patient_df), "_"), function(x)x[1])
patient_df$MSI = factor(patient_df$MSI, levels = c("MSS", "MSI"))
patient_df = patient_df[order(patient_df$Sample_type, patient_df$MSI, patient_df$Stage),]

show_col = c("Stage","MSI")
pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap_wtmeta", ".pdf"), width = 4.33, height = 3.5)
pheatmap(compos_rel[rownames(patient_df),show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, 
         color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList,
         annotation_row = patient_df[,show_col])
dev.off()

pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap2_wtmeta", ".pdf"), width = 4.33, height = 3.5)
pheatmap(compos_rel[rownames(patient_df),show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=6, 
         color = colorRampPalette(c("white", "#3690c0"))(length(breaksList)), breaks = breaksList,
         annotation_row = patient_df[,show_col])
dev.off()


write.xlsx(compos, paste0(msheetPath,test_text,"compos_dataset_heatmap.xlsx"), rowNames = T)
write.xlsx(patient_df, paste0(msheetPath,test_text,"compos_dataset_heatmap_patient_df.xlsx"), rowNames=T)

# Figure 2, epithelial analysis

## see analysis_epithelial.R

## DE plot see de_two_way.R



# Mouse vivo

test_text <- "mouse_invivo_"
clist_mouse <- readRDS("../mouse_colon_cello/clist.rds")
eset_mouse <- readRDS("../mouse_colon_cello/eset.rds")
cur_vis <- clist_mouse$`Mouse in vivo [IFF, low mito]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset_mouse)[cur_vis@idx,])

cur_proj$Cell_type_refined[cur_proj$Cell_type_refined %in% c("Epithelial cells [EC,TA,SC]", "EEC cells", "Goblet cells")] = "Epithelial cells (normal)"
cur_proj$Cell_type_refined[cur_proj$Cell_type_refined == "Epithelial cells (tumor)"] <- paste0("Epithelial cells (", cur_proj$Dataset[cur_proj$Cell_type_refined == "Epithelial cells (tumor)"], ")")
ctype_color <- c(
    "Epithelial cells (normal)" = "#1f78b4",
    "Epithelial cells (Tumor AOM)" = "#e31a1c",
    "Epithelial cells (Tumor ApcMin)" = "#bf5b17",
    "Epithelial cells (Tumor APKS)" = "#df65b0",
    "Fibroblast" = "#cab2d6", 
    "Endothelial"= "#33a02c",
    "Myeloid" = "#6a3d9a",
    "T cells" ="#a6cee3",
    "B & Plasma cell" = "#ff7f00"
)

cur_proj$Cell_type_refined <- factor(cur_proj$Cell_type_refined, levels = names(ctype_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_type_refined", pal=ctype_color, size = .3, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "cell_type.pdf"), g1, width = 4, height=1.7, units = "in", device = "pdf")

write.xlsx(cur_proj[,c("UMAP_1", "UMAP_2", "Cell_type_refined")], paste0(msheetPath,test_text,"cell_type.xlsx"), rowNames = T)



cur_proj <- cbind(cur_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset_mouse)[cur_vis@idx,])
cur_proj$Cell_type_refined[cur_proj$Cell_type_refined %in% c("Epithelial cells [EC,TA,SC]", "EEC cells", "Goblet cells")] = "Epithelial cells (normal)"
cur_proj$Cell_type_refined <- factor(cur_proj$Cell_type_refined, levels = c(
        "Epithelial cells (normal)",
        "Epithelial cells (tumor)",
        "Fibroblast", 
        "Endothelial",
        "Myeloid",
        "T cells",
        "B & Plasma cell"
))
compos<-as.data.frame.matrix(table(cur_proj$Dataset, cur_proj$Cell_type_refined))
compos <- compos[!grepl("Flox", rownames(compos)),]
compos_rel <- compos / (compos$`Epithelial cells (normal)` + compos$`Epithelial cells (tumor)`)
breaksList <- seq(0,1,by=.1)
pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap", ".pdf"), width = 3.3, height = 2.2)
pheatmap(compos_rel, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)
dev.off()

pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap2", ".pdf"), width = 3.3, height = 2.2)
pheatmap(compos_rel, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#3690c0"))(length(breaksList)), breaks = breaksList)
dev.off()

write.xlsx(compos, paste0(msheetPath,test_text,"compos_dataset_heatmap.xlsx"), rowNames = T)
























####### CODE BELOW THIS LINE NOT USED ##########



test_text <- "mono_hs_invivo_"
cur_vis <- clist$`Monocyte/macrophage invivo`
cur_proj <- cbind(cur_vis@proj$`PHATE-2D`, pData(eset)[cur_vis@idx,])
stype_color <- c(
    "Colon" = "#6a3d9a",   
    "Liver" = "#cab2d6", 
    "Tumor" = "#ff7f00",
    "LiverMet" = "#fdbf6f"
)
cur_proj$SampleType <- factor(cur_proj$SampleType, levels = names(stype_color))
cur_proj$PHATE1 <- - cur_proj$PHATE1
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="SampleType", pal=stype_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "sampletype.pdf"), g1, width = 3.4, height=1.6, units = "in", device = "pdf")


plot_g <- c("S100A8", "SPP1", "C1QC")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(plot_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- plot_g

glist <- list()
for(g in plot_g) {
    print(g)
    cur_proj$gene_expr <- g_exprs[[g]]
    ecut <- quantile(cur_proj$gene_expr[cur_proj$gene_expr > 0], .975)
    if(ecut < 1) ecut = 1
    cur_proj$gene_expr[cur_proj$gene_expr > ecut] = ecut
    glist[[g]]<-plotProj(cur_proj, dim_col = c(1,2), group.by="gene_expr", pal="viridis", size = .6, legend.title = g, layer_zero = T) + guides(color = guide_colorbar(barwidth = 10, barheight = 1, title = g)) +
        guides(color = guide_colorbar(
            barwidth = 4, barheight = .5, title = g,
            title.theme = element_text(size = 8),
            label.theme = element_text(size = 8))) + 
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        theme(text=element_text(family = "Helvetica", size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0))+
        monocle:::monocle_theme_opts()
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "gexpr_umap.eps"), g1, width = 6, height=2, units = "in", device = "eps")










# Receptor missing in tumoroid
test_text <- "missing_rl_myeloid_"
rec_bind <- readRDS(paste0(mstatePath, "human_rl_", "rec_bind.rds"))
lig_bind <- readRDS(paste0(mstatePath, "human_rl_", "lig_bind.rds"))

rec_myeloid <- rec_bind[rec_bind$source == "Myeloid",]
lig_myeloid <- lig_bind[lig_bind$source == "Myeloid",]

tumoroid28_de_tbl <- read_excel_allsheets(paste0(msheetPath, "Epithelial cells [wt OGN, IFF]_Tumor_28_vs_Tumoroid_28_24hr_Tumoroid_28_48hr_2020-07-20_de_all.xlsx"))
lig_myeloid$LFC_tumoroid28 <- tumoroid28_de_tbl$Tumor_28$log2fc[match(lig_myeloid$receptor, tumoroid28_de_tbl$Tumor_28$gene_short_name)]
lig_myeloid$Padj_tumoroid28 <- tumoroid28_de_tbl$Tumor_28$p_adj[match(lig_myeloid$receptor, tumoroid28_de_tbl$Tumor_28$gene_short_name)]

rec_myeloid$LFC_tumoroid28 <- tumoroid28_de_tbl$Tumor_28$log2fc[match(rec_myeloid$ligand, tumoroid28_de_tbl$Tumor_28$gene_short_name)]
rec_myeloid$Padj_tumoroid28 <- tumoroid28_de_tbl$Tumor_28$p_adj[match(rec_myeloid$ligand, tumoroid28_de_tbl$Tumor_28$gene_short_name)]

lig_myeloid<- lig_myeloid[complete.cases(lig_myeloid),]
lig_myeloid$pos = ifelse(lig_myeloid$LFC_tumoroid28 > 0, 1,0)
verts <- data.frame(name = unique(c(as.character(lig_myeloid$ligand), as.character(lig_myeloid$receptor))))
verts$type <- verts$name %in% lig_myeloid$ligand
lig_graph<-graph_from_data_frame(lig_myeloid, vertices = verts)
pdf(paste0(mplotPath, test_text, "tumor_receptor_lost.pdf"), width = 3, height = 3)
ggraph(lig_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(color = factor(pos), edge_width = LFC_tumoroid28), show.legend = F) + 
    scale_edge_color_manual(values = c("1" = "black", "0" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    scale_edge_width(range = c(0.2,1)) + 
    geom_node_point(aes(shape = type, color = type), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = type),  show.legend = F, label.size = 0.25, size= 2) + 
    scale_color_manual(values = c("TRUE" = "#6a3d9a", "FALSE" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()




rec_myeloid<- rec_myeloid[complete.cases(rec_myeloid),]
rec_myeloid$pos = ifelse(rec_myeloid$LFC_tumoroid28 > 0, 1,0)
verts <- data.frame(name = unique(c(as.character(rec_myeloid$ligand), as.character(rec_myeloid$receptor))))
verts$type <- verts$name %in% rec_myeloid$receptor
rec_graph<-graph_from_data_frame(rec_myeloid, vertices = verts)
pdf(paste0(mplotPath, test_text, "tumor_ligand_lost.pdf"), width = 3, height = 3)
ggraph(rec_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(color = factor(pos), edge_width = LFC_tumoroid28), show.legend = F) + 
    scale_edge_color_manual(values = c("1" = "black", "0" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    scale_edge_width(range = c(0.2,1)) + 
    geom_node_point(aes(shape = type, color = type), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = type),  show.legend = F, label.size = 0.25, size= 2) + 
    scale_color_manual(values = c("TRUE" = "#6a3d9a", "FALSE" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()


