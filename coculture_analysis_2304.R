

library(VisCello)
mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

savePath <- "../hcc_2023dataset_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)

scriptPath <- paste0("../hcc_cello/scripts/")
source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

clist <- readRDS("../hcc_2023dataset_cello/clist.rds")
eset <- readRDS("../hcc_2023dataset_cello/eset.rds")
library(ggrastr)

test_text <- "coculture_2023_"
cur_vis <- clist$`Coculture23 singlets [low mito]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Louvain_Cluster
#cur_proj <- cur_proj[cur_proj$Cluster != "8",]
dset_color <- get_factor_color(unique(cur_proj$Dataset))
names(dset_color) <- unique(cur_proj$Dataset)
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Dataset", pal=dset_color, size = .4, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "umap_dataset.pdf"), g1, width = 3.5, height=2.5, units = "in", device = "pdf")

write.xlsx(cur_proj[,c("UMAP_1", "UMAP_2", "Dataset", "SampleType", "Cell_type")], paste0(msheetPath, test_text, "umap_dataset.xlsx"), rowNames=T)


stype_color <- get_factor_color(unique(cur_proj$SampleType))
names(stype_color) <- unique(cur_proj$SampleType)
stype_levels = c("O", "T", "M", "F", "M+O", "M+T", "M+F", "F+M+O", "F+M+T", "M+CMO", "M+CMT")
cur_proj$SampleType <- factor(cur_proj$SampleType, levels = stype_levels)
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="SampleType", pal=stype_color, size = .4, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "umap_stype.pdf"), g1, width = 3, height=2.2, units = "in", device = "pdf")



use_color <- get_factor_color(unique(cur_proj$Patient))
names(use_color) <- unique(cur_proj$Patient)
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Patient", pal=use_color, size = .4, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "umap_patient.pdf"), g1, width = 3, height=2.2, units = "in", device = "pdf")



mac_proj <- cur_proj[which(cur_proj$Cell_type == "Macrophage"),]
g1<-plotProj(mac_proj, dim_col = c(1,2), group.by="Dataset", pal=dset_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "maczoom_dataset.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")

st2_mapping <- c("M"="M", "M+F"="M(+F)", "M+T" = "M(+T)",  "M+O" = "M(+O)", "F+M+T" = "M(+F+T)", "F+M+O" = "M(+F+O)","M+CMT" = "M+(CMT)", "M+CMO" ="M+(CMO)")
st2_color <- stype_color[names(st2_mapping)]
names(st2_color) = as.character(st2_mapping)
mac_proj$SampleType2 <- factor(st2_mapping[as.character(mac_proj$SampleType)], levels = as.character(st2_mapping))
g1<-plotProj(mac_proj, dim_col = c(1,2), group.by="SampleType2", pal=st2_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "maczoom_stype.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")

nbins=5
glist <- list()

for(st in as.character(st2_mapping)) {
    print(st)
    glist[[st]]<-ggplot(mac_proj, aes_string("UMAP_1", "UMAP_2")) + 
        geom_point_rast(color = 'lightgrey', size = .1) + 
        stat_density2d(data = mac_proj[mac_proj$SampleType2 %in% st,], aes(fill = SampleType2, color = SampleType2), geom="polygon", linewidth=.4, alpha = .3, bins = nbins) +
        scale_color_manual(values = st2_color)+
        scale_fill_manual(values = st2_color)+
        xlab("UMAP 1") +
        ylab("UMAP 2") + 
        theme(text=element_text(family = "Helvetica", size=8),
              legend.text=element_text(size=8),
              legend.key.size = unit(.3, 'cm'),
              axis.text = element_text(size=9),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0),
              #legend.box.margin=margin(-10,-10,-10,-10),
              plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
        monocle:::monocle_theme_opts() 
}
plot_cbn <- arrangeGrob(grobs = glist, ncol= 2)
ggsave(paste0(mplotPath, test_text, "maczoom_stype_distribution.pdf"), plot_cbn, width = 4, height=3.3, units = "in", device = "pdf")


# signal enrichment

de_list <- read_excel_allsheets(paste0(msheetPath, "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))
max_g <- 100
reactomeSets <- lapply(de_list, function(x) {
    sig_g<-x$gene_short_name
    sig_g[1:min(max_g, length(sig_g))]
})
library(AUCell)
cur_eset<- eset[, rownames(mac_proj)]
cur_eset <- cur_eset[rowMeans(exprs(cur_eset) > 0) > .01,]
expr_hc <- as.matrix(exprs(cur_eset))
rownames(expr_hc) <- fData(cur_eset)$gene_short_name
cells_rankings_hc <- AUCell_buildRankings(expr_hc)
saveRDS(cells_rankings_hc, paste0(rdsPath, test_text, "cells_rankings_hc.rds"))

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(reactomeSets, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
#saveRDS(cells_AUC, paste0(rdsPath, test_text, cur_thresh, "cells_AUC", "_wt_ifg_cut75.rds"))
saveRDS(cells_AUC,paste0(mstatePath, test_text, "aucell_reactome_0.2.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "aucell_reactome_0.2.rds"))
matrixAUC <- getAUC(cells_AUC)
plot_meta <-  mac_proj
saveRDS(plot_meta, paste0(rdsPath, test_text, "plot_meta.rds"))

plot_meta <- readRDS(paste0(rdsPath, test_text, "plot_meta.rds"))
plot_meta$Patient <- factor(plot_meta$Patient, levels = c("8", "24"))

save_proj <- cbind(mac_proj, t(matrixAUC))
write.xlsx(save_proj[c("UMAP_1", "UMAP_2","Dataset", "SampleType2", "Cell_type", "Patient", "C1QC", "IL1B", "SPP1")], paste0(msheetPath,test_text, "mac_save_proj.xlsx"), rowNames=T)

glist <- list()
for(show_p in names(reactomeSets)) {
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    cut_qt = .975
    plot_df$activity_score[plot_df$activity_score > quantile(plot_df$activity_score, cut_qt)] <-  quantile(plot_df$activity_score, cut_qt)
    glist[[show_p]] <- plotProj(plot_df, dim_col = which(colnames(plot_df) %in% c("UMAP_1", "UMAP_2")), group.by = "activity_score", pal = "BlueGreenRed", size = .7) +
        guides(color = guide_colorbar(
            title = NULL,
            barwidth = 6.5, barheight = .7, 
            title.theme = element_text(size = 8),
            label.theme = element_text(size = 8))) + 
        theme(text=element_text(family = "Helvetica", size=9),
              axis.text = element_text(family = "Helvetica", size=8),
              legend.position="top",
              legend.title=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0))+
        ggtitle(show_p) +
        monocle:::monocle_theme_opts()
}


pdf(paste0(mplotPath, test_text, "activity_umap_cbn", ".pdf"),  width = 8, height=2.2)
grid.draw(arrangeGrob(grobs = glist, nrow = 1))
dev.off()


glist <- list()

for(show_p in names(reactomeSets)) {
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    #plot_df <- plot_df[plot_df$Dataset2 != "T96+M_24hr",]
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    plot_df$activity_score_scaled <- scale(plot_df$activity_score)
    dup_df <- plot_df[plot_df$SampleType2 %in% c("M", "M(+F)"),]
    dup_df1 <- dup_df; dup_df1$Patient = "8"
    dup_df2 <- dup_df; dup_df2$Patient = "24"
    plot_df <- rbind(plot_df, dup_df1, dup_df2)
    plot_df <- plot_df[!is.na(plot_df$Patient),]
    # 
    # plot_df$activity_score_scaled[which(plot_df$Patient == 2)] <- scale(plot_df$activity_score[which(plot_df$Patient == 2)])
    # plot_df$activity_score_scaled[which(plot_df$Patient == 6)] <- scale(plot_df$activity_score[which(plot_df$Patient == 6)])
    
    glist[[show_p]] <- ggplot(plot_df, aes_string(x = "SampleType2", y = "activity_score_scaled", fill = "SampleType2")) +
        geom_boxplot(outlier.colour = NA)+
        scale_fill_manual(values = st2_color)+
        ggtitle(show_p)+
        ylim(c(-2.5,2.5)) + 
        xlab(NULL)+
        ylab("Activity score") +
        theme_bw() + 
        theme(text=element_text(family = "Helvetica", size=8),
              axis.text = element_text(family = "Helvetica", size=8),
              axis.ticks.x=element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1),
              strip.text.x = element_text(size = 9),
              legend.key.size = unit(.4, 'cm'),
              legend.text=element_text(size=8),
              legend.spacing.y = unit(.3, 'cm'),
              legend.margin=margin(0,0,0,0))+
        facet_wrap(~Patient)
}
pdf(paste0(mplotPath, test_text, "macactivity_stype", ".pdf"),  width =4, height=2.4)
invisible(lapply(glist, print))
dev.off()

pres_all <- list()
pres_all_star <- list()
compare_groups <- list(
    c("M", "M(+F)"),
    c("M(+T)",  "M(+O)"),
    c("M(+F+T)", "M(+F+O)"),
    c("M+(CMT)", "M+(CMO)")
)
for(show_p in names(reactomeSets)) {
    print(show_p)
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    #plot_df <- plot_df[plot_df$Dataset2 != "T96+M_24hr",]
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    plot_df$activity_score_scaled <- scale(plot_df$activity_score)
    dup_df <- plot_df[plot_df$SampleType2 %in% c("M", "M(+F)"),]
    dup_df1 <- dup_df; dup_df1$Patient = "8"
    dup_df2 <- dup_df; dup_df2$Patient = "24"
    plot_df <- rbind(plot_df, dup_df1, dup_df2)
    plot_df <- plot_df[!is.na(plot_df$Patient),]

    pres_state <- sapply(levels(plot_df$Patient), function(cur_patient) {
        print(cur_patient)
        cur_tbl = plot_df[plot_df$Patient == cur_patient,]
        pres_patient <- list()
        for(cur_group in compare_groups) {
            print(cur_group)
            pres_patient[[paste0(make.names(cur_group), collapse = "_vs_")]] <- t.test(cur_tbl$activity_score[cur_tbl$SampleType2 == cur_group[1]], cur_tbl$activity_score[cur_tbl$SampleType2 == cur_group[2]])
        }
        sapply(pres_patient, function(x) x$p.value)
    })

    pres_state_star <- pres_state
    pres_state_star[pres_state > 0.05] = "ns"
    pres_state_star[pres_state <= 0.05] = "*"
    pres_state_star[pres_state <= 0.01] = "**"
    pres_state_star[pres_state <= 0.001] = "***"
    pres_state_star[pres_state <= 0.0001] = "****"
    
    pres_all[[show_p]] <- pres_state
    pres_all_star[[show_p]] <- pres_state_star
}

write.xlsx(pres_all, paste0(msheetPath,test_text, "macrophage_state_t_test_res.xlsx"), rowNames=T)
write.xlsx(pres_all_star, paste0(msheetPath,test_text, "macrophage_state_t_test_res_star.xlsx"), rowNames=T)





# Pathways
set.seed(2020)
dep_vis <- clist$`Coculture23 Macrophage [low mito, seurat]`
cur_eset <- eset[, dep_vis@idx]
test_pairs <- list(
    c("M+T", "M+CMT")
)
eset_de <- cur_eset[, which(cur_eset$SampleType %in% unlist(test_pairs))]

g1_count <- 1000
g2_count <- 1000

eset_de$de_group <- eset_de$SampleType
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

pdeg_list<- list()
for(test_pair in test_pairs){
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), g1_count))
    group2_idx <-  which(pData(eset_de)$de_group == test_pair[2])
    group2_idx <- sample(group2_idx, min(length(group2_idx), g2_count))
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
    pdeg_list[[test_pair[1]]] <- prioritized_genes[[1]]
    pdeg_list[[test_pair[2]]] <- prioritized_genes[[2]]
}

saveRDS(pdeg_list, paste0(mstatePath, test_text, "pdeg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))

lapply(deg_list, nrow)


plot_meta <- pData(eset_de)
plot_meta <- plot_meta[, c("Patient", "SampleType"), drop=F]
plot_meta[["SampleType"]] <- factor(plot_meta[["SampleType"]], levels = st_order)
max_cell = 500
sp_idx <- lapply(st_order, function(x) {sample(which(plot_meta[["SampleType"]] == x), min(max_cell, sum(plot_meta[["SampleType"]] == x, na.rm = T)))})
plot_meta <- plot_meta[unlist(sp_idx),, drop=F]
plot_meta <- plot_meta[order(plot_meta[["SampleType"]], plot_meta$Patient),]
non_na_cells_ordered <- rownames(plot_meta)

deg_sheets <- deg_list
max_gene <- 15
plotg <- unique(unlist(lapply(deg_list, function(x) {
    x = x %>% filter(common_mean >=.1)
    use_g <- as.character(x$gene_id)
    use_g[1:max_gene]
})))

value <- eset_de@assayData$norm_exprs[plotg,non_na_cells_ordered]
rownames(value) <- fData(eset_de)$gene_short_name[match(rownames(value), rownames(eset_de))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

pdf(paste0(mplotPath, test_text, "mac_M_T_v_M_CMT_deg_", max_gene,".pdf"), width = 6, height=4)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = T, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 8)
dev.off()




show_group = "M+T" # show_group = "M+CMT"
res_list <- lapply(deg_list[names(deg_list) == show_group], function(x) {
    x$gene_id
})
gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
gene_background <- rownames(eset_de)
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)

include_g.df <- bitr(unlist(res_list), fromType = gene_id_type,
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = orgdb)
go_res<- enrichGO(gene        = include_g.df$ENTREZID,
                  universe      = bg.df$ENTREZID,
                  OrgDb         = orgdb,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

saveRDS(go_res, paste0(mstatePath, test_text, show_group, "go_res.rds"))
go_res_tbl <- go_res@result
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, show_group, "go_res.xlsx"))

# GO plot
include_go_simplified <- clusterProfiler::simplify(go_res, cutoff = .5)
plot_df<- include_go_simplified@result[1:min(nrow(include_go_simplified@result),10),]
plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
g1 <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
    geom_bar(stat = "identity", aes(fill = -log10(qvalue))) + 
    #scale_x_discrete(position = 'top') + 
    coord_flip() +
    xlab("") + 
    theme(
        text=element_text(family = "Helvetica", size=8,color="black"),
        axis.text = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
    ggtitle(paste0("Up-regulated in ", show_group)) + 
    monocle:::monocle_theme_opts() 

pdf(paste0(mplotPath, test_text, "macrophagede_", show_group, "go_bar.pdf"), width = 6, height = 2)
g1
dev.off()





show_patient = "8"
#show_patient = "24"

show_type =  "tumoroid_"
st_order <- c("T", "M+T", "F+M+T")
plot_proj <- cur_proj[which(cur_proj$Cell_type == "Epithelial(Tumoroid)" & cur_proj$Patient == show_patient & cur_proj$SampleType != "M+CMT"),]
# plot_proj <- plot_proj[!plot_proj$Cluster %in% names(which(table(plot_proj$Cluster) < 50)),]

show_type =  "organoid_"
st_order <- c("O", "M+O", "F+M+O")
plot_proj <- cur_proj[which(cur_proj$Cell_type == "Epithelial(Organoid)" & cur_proj$Patient == show_patient & cur_proj$SampleType != "M+CMO"),]
# plot_proj <- plot_proj[!plot_proj$Cluster %in% names(which(table(plot_proj$Cluster) < 50)),]

g1<-plotProj(plot_proj, dim_col = c(1,2), group.by="SampleType", pal=stype_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text,show_type, show_patient, "_zoom_stype.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")

nbins=5
glist <- list()

for(st in st_order) {
    glist[[st]]<-ggplot(plot_proj, aes_string("UMAP_1", "UMAP_2")) + 
        geom_point_rast(color = 'lightgrey', size = .1) + 
        stat_density2d(data = plot_proj[plot_proj$SampleType %in% st,], aes(fill = SampleType, color = SampleType), geom="polygon", linewidth=.4, alpha = .3, bins = nbins) +
        scale_color_manual(values = stype_color)+
        scale_fill_manual(values = stype_color)+
        xlim(c(min(plot_proj$UMAP_1) - .5, max(plot_proj$UMAP_1) + .5))+
        ylim(c(min(plot_proj$UMAP_2) - .5, max(plot_proj$UMAP_2) + .5))+
        xlab("UMAP 1") +
        ylab("UMAP 2") + 
        theme(text=element_text(family = "Helvetica", size=8),
              legend.text=element_text(size=8),
              legend.key.size = unit(.3, 'cm'),
              axis.text = element_text(size=9),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0),
              #legend.box.margin=margin(-10,-10,-10,-10),
              plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
        monocle:::monocle_theme_opts() 
}
plot_cbn <- arrangeGrob(grobs = glist, ncol= 1)
ggsave(paste0(mplotPath, test_text, show_type, show_patient, "_zoom_stype_distribution.pdf"), plot_cbn, width = 2, height=3.3, units = "in", device = "pdf")


# cc_frac<-as.data.frame.matrix(table(plot_proj$SampleType, plot_proj$Cell_cycle_phase))
# cc_frac <- cc_frac/rowSums(cc_frac)
# cc_frac$G2M_S <- cc_frac$G2M+ cc_frac$S
# cc_frac$SampleType <- factor(rownames(cc_frac), levels = c("T","M+T", "F+M+T"))
# g1 <- ggplot(cc_frac, aes_string(x = "SampleType", y = "G2M_S", fill = "SampleType")) +
#     geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
#     scale_fill_manual(values = stype_color)+
#     theme_bw() +
#     theme(text=element_text(family = "Helvetica", size=9),
#           axis.text = element_text(family = "Helvetica", size=9),
#           axis.ticks.x=element_blank(),
#           axis.text.x = element_text(angle = 90, hjust = 1),
#           legend.margin=margin(0,0,0,0))+
#     ggtitle("G2M/S fraction")+
#     xlab(NULL)
# 
# pdf(paste0(mplotPath, test_text, "tumoroid_", show_patient, "G2M_S_frac", ".pdf"),  width = 2, height=1.5)
# g1
# dev.off()
# 















cur_vis <- clist$`Coculture23 singlets [low mito]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Louvain_Cluster
cur_proj <- cur_proj[cur_proj$Cluster != "8",]
epi_proj <- cur_proj[grepl("Epithelial", cur_proj$Cell_type) & (!cur_proj$Dataset %in% c("M+CMO_8","M+CMO_24","M+CMT_8","M+CMT_24")),]
epi_proj$SampleType2 <- epi_proj$SampleType

epi_proj$Patient <- factor(epi_proj$Patient, c("8","24"))

st2_mapping <- c("O" = "O", "M+O" = "O(+M)", "F+M+O" = "O(+M+F)", "T" = "T", "M+T" = "T(+M)", "F+M+T" = "T(+M+F)")
epi_proj$SampleType2 <- factor(st2_mapping[epi_proj$SampleType], levels = as.character(st2_mapping))

show_genes <- c("VIM", "CD44")
cur_gene <- show_genes[[1]]
#cur_gene <- show_genes[[2]]
epi_proj$gene_expr <- eset@assayData$norm_exprs[which(fData(eset)$gene_short_name == cur_gene),rownames(epi_proj)]

st2_color <- stype_color[names(st2_mapping)]
names(st2_color) = as.character(st2_mapping)
plot_df <- epi_proj
g1 <- ggplot(plot_df, aes_string(x = "SampleType2", y = "gene_expr", fill = "SampleType2")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = st2_color)+
    theme_bw() + 
    ggtitle(cur_gene)+
    xlab(NULL)+
    ylab("Expression level")+
    ylim(c(0,3)) + 
    theme(text=element_text(family = "Helvetica", size=7),
          axis.text = element_text(family = "Helvetica", size=7),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x = element_text(size = 7, margin = margin(.1,0,.1,0, "cm")),
          legend.key.size = unit(.4, 'cm'),
          legend.text=element_text(size=8),
          legend.spacing.y = unit(.3, 'cm'),
          legend.margin=margin(0,0,0,0))+
    facet_wrap(~Patient)

pdf(paste0(mplotPath, test_text, cur_gene , "_expr", ".pdf"),  width = 3, height=1.9)
g1
dev.off()

# t test on VIM/CD44
compare_groups <- list(
    c("O", "O(+M)"), 
    c("O", "O(+M+F)"), 
    c("O(+M)", "O(+M+F)"), 
    c("T", "T(+M)"), 
    c("T", "T(+M+F)"), 
    c("T(+M)", "T(+M+F)")
)

test_group = "gene_expr"

pres_state <- sapply(levels(epi_proj$Patient), function(cur_patient) {
    print(cur_patient)
    cur_tbl = epi_proj[epi_proj$Patient == cur_patient,]
    pres_patient <- list()
    for(cur_group in compare_groups) {
        print(cur_group)
        pres_patient[[paste0(make.names(cur_group), collapse = "_vs_")]] <- t.test(cur_tbl[[test_group]][cur_tbl$SampleType2 == cur_group[1]], cur_tbl[[test_group]][cur_tbl$SampleType2 == cur_group[2]], alternative = "less")
    }
    sapply(pres_patient, function(x) x$p.value)
})

pres_state_star <- pres_state
pres_state_star[pres_state > 0.05] = "ns"
pres_state_star[pres_state <= 0.05] = "*"
pres_state_star[pres_state <= 0.01] = "**"
pres_state_star[pres_state <= 0.001] = "***"
pres_state_star[pres_state <= 0.0001] = "****"


write.xlsx(list(test_res = pres_state), paste0(msheetPath, test_text, cur_gene, "_expr_ttest_res.xlsx"), rowNames = T)
write.xlsx(list(test_res = pres_state_star), paste0(msheetPath, test_text, cur_gene, "_expr_ttest_res_star.xlsx"), rowNames = T)


## Cell cycle analysis for epithelial cells
cc_frac<-as.data.frame.matrix(table(epi_proj$Dataset, epi_proj$Cell_cycle_phase))
cc_frac <- cc_frac[!rownames(cc_frac) %in% c("M+CMO_8","M+CMO_24","M+CMT_8","M+CMT_24"),]
cc_frac <- cc_frac/rowSums(cc_frac)
cc_frac$G2M_S <- cc_frac$G2M+ cc_frac$S
cc_frac$SampleType2 <- epi_proj$SampleType2[match(rownames(cc_frac), epi_proj$Dataset)]
cc_frac$Patient <- epi_proj$Patient[match(rownames(cc_frac), epi_proj$Dataset)]
g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G2M_S", fill = "SampleType2")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = st2_color)+
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.margin=margin(0,0,0,0))+
    ggtitle("G2M/S fraction")+
    xlab(NULL)+
    theme_bw() 
pdf(paste0(mplotPath, test_text, "G2M_S_frac", ".pdf"),  width = 4, height=2.7)
g1
dev.off()
g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G1", fill = "SampleType2")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = st2_color)+
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.margin=margin(0,0,0,0))+
    ggtitle("G1 fraction")+
    xlab(NULL)+
    theme_bw() 
pdf(paste0(mplotPath, test_text, "G1_frac", ".pdf"),  width = 4, height=2.7)
g1
dev.off()

cc_epi <- as.data.frame(c(
    "M+O_8" = cc_frac["M+O_8","G2M_S"] - cc_frac["O_8","G2M_S"],
    "M+O_24" = cc_frac["M+O_24","G2M_S"] - cc_frac["O_24","G2M_S"],
    "M+T_8" = cc_frac["M+T_8","G2M_S"] - cc_frac["T_8","G2M_S"],
    "M+T_24" = cc_frac["M+T_24","G2M_S"] - cc_frac["T_24","G2M_S"],
    "F+M+O_8" = cc_frac["F+M+O_8","G2M_S"] - cc_frac["O_8","G2M_S"],
    "F+M+O_24" = cc_frac["F+M+O_24","G2M_S"] - cc_frac["O_24","G2M_S"],
    "F+M+T_8" = cc_frac["F+M+T_8","G2M_S"] - cc_frac["T_8","G2M_S"],
    "F+M+T_24" = cc_frac["F+M+T_24","G2M_S"] - cc_frac["T_24","G2M_S"]
))
colnames(cc_epi) = "cc_diff"
cc_epi$Patient <- epi_proj$Patient[match(rownames(cc_epi), epi_proj$Dataset)]
cc_epi$Patient <- factor(cc_epi$Patient, levels = c("8", "24"))
cc_epi$SampleType2 <- epi_proj$SampleType2[match(rownames(cc_epi), epi_proj$Dataset)]
g1<-ggplot(cc_epi, aes_string(x = "SampleType2", y = "cc_diff", fill = "SampleType2")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = st2_color) + 
    xlab(NULL) + 
    theme_bw() + 
    facet_wrap(~Patient) + 
    theme(text=element_text(family = "Helvetica", size=8),
          axis.text = element_text(family = "Helvetica", size=8),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x = element_text(size = 9),
          legend.key.size = unit(.4, 'cm'),
          legend.text=element_text(size=8),
          legend.spacing.y = unit(.3, 'cm'),
          legend.margin=margin(0,0,0,0))

pdf(paste0(mplotPath, test_text, "cc_reduce", ".pdf"),  width = 3.5, height=2)
g1
dev.off()






# Cell cycle analysis for macrophage
# mac_proj <- cur_proj[which(cur_proj$Cell_type == "Macrophage"),]
# 
# cc_frac<-as.data.frame.matrix(table(mac_proj$Dataset, mac_proj$Cell_cycle_phase))
# 
# cc_frac <- cc_frac/rowSums(cc_frac)
# cc_frac$G2M_S <- cc_frac$G2M+ cc_frac$S
# cc_frac$SampleType <- mac_proj$SampleType[match(rownames(cc_frac), mac_proj$Dataset)]
# cc_frac$Patient <- mac_proj$Patient[match(rownames(cc_frac), mac_proj$Dataset)]
# 
# g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G2M_S", fill = "SampleType")) +
#     geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
#     scale_fill_manual(values = stype_color)+
#     theme(text=element_text(family = "Helvetica", size=9),
#           axis.text = element_text(family = "Helvetica", size=9),
#           axis.ticks.x=element_blank(),
#           axis.text.x=element_blank(),
#           legend.margin=margin(0,0,0,0))+
#     ggtitle("G2M/S fraction")+
#     xlab(NULL)+
#     theme_bw() 
# pdf(paste0(mplotPath, test_text, "mac_G2M_S_frac", ".pdf"),  width = 4, height=2.7)
# g1
# dev.off()
# g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G1", fill = "SampleType")) +
#     geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
#     scale_fill_manual(values = stype_color)+
#     theme(text=element_text(family = "Helvetica", size=9),
#           axis.text = element_text(family = "Helvetica", size=9),
#           axis.ticks.x=element_blank(),
#           axis.text.x=element_blank(),
#           legend.margin=margin(0,0,0,0))+
#     ggtitle("G1 fraction")+
#     xlab(NULL)+
#     theme_bw() 
# pdf(paste0(mplotPath, test_text, "mac_G1_frac", ".pdf"),  width = 4, height=2.7)
# g1
# dev.off()
# 
# cc_mac <- as.data.frame(c(
#     "M+O_8" = cc_frac["M+O_8","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+O_24" = cc_frac["M+O_24","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+T_8" = cc_frac["M+T_8","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+T_24" = cc_frac["M+T_24","G2M_S"] - cc_frac["M","G2M_S"],
#     "F+M+O_8" = cc_frac["F+M+O_8","G2M_S"] - cc_frac["M","G2M_S"],
#     "F+M+O_24" = cc_frac["F+M+O_24","G2M_S"] - cc_frac["M","G2M_S"],
#     "F+M+T_8" = cc_frac["F+M+T_8","G2M_S"] - cc_frac["M","G2M_S"],
#     "F+M+T_24" = cc_frac["F+M+T_24","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+CMO_8" = cc_frac["M+CMO_8","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+CMO_24" = cc_frac["M+CMO_24","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+CMT_8" = cc_frac["M+CMT_8","G2M_S"] - cc_frac["M","G2M_S"],
#     "M+CMT_24" = cc_frac["M+CMT_24","G2M_S"] - cc_frac["M","G2M_S"]
# ))
# 
# colnames(cc_mac) = "cc_diff"
# cc_mac$Patient <- epi_proj$Patient[match(rownames(cc_mac), epi_proj$Dataset)]
# cc_mac$SampleType <- epi_proj$SampleType[match(rownames(cc_mac), epi_proj$Dataset)]
# 
# g1<-ggplot(cc_mac, aes_string(x = "SampleType", y = "cc_diff", fill = "SampleType")) +
#     geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
#     scale_fill_manual(values = stype_color) + 
#     xlab(NULL) + 
#     theme_bw() + 
#     facet_wrap(~Patient) + 
#     theme(text=element_text(family = "Helvetica", size=8),
#           axis.text = element_text(family = "Helvetica", size=8),
#           axis.ticks.x=element_blank(),
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           strip.text.x = element_text(size = 9),
#           legend.key.size = unit(.4, 'cm'),
#           legend.text=element_text(size=8),
#           legend.spacing.y = unit(.3, 'cm'),
#           legend.margin=margin(0,0,0,0))
# 
# pdf(paste0(mplotPath, test_text, "cc_mac", ".pdf"),  width = 3, height=2)
# g1
# dev.off()







# EMT
library(AUCell)
emt_genes <- read.delim("~/Documents/CHOP/NingProj/VisCello.cc/data-raw/public_resource/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v7.5.1.tsv", row.names=1)
emt_genes <- strsplit(emt_genes["MAPPED_SYMBOLS",], ",")
names(emt_genes) <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"

cur_eset <- eset[,rownames(epi_proj)]
expressed_g <- which(rowMeans(exprs(cur_eset) > 0) > .01)
expr_hc <- as.matrix(exprs(cur_eset[expressed_g,]))
rownames(expr_hc) <- make.names(fData(cur_eset)$gene_short_name[match(rownames(expr_hc), fData(cur_eset)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_hc)

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(emt_genes, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "emt_0.2cells_AUC.rds"))
cells_AUC <- readRDS(paste0(mstatePath, test_text, "emt_0.2cells_AUC.rds"))
plot_df <- as.data.frame(t(getAUC(cells_AUC)))
epi_proj$HALLMARK_EMT <- plot_df$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
g1 <- ggplot(epi_proj, aes(x = SampleType2, y = HALLMARK_EMT)) +
    geom_boxplot(aes(fill = SampleType2), outlier.shape=NA) + 
    scale_fill_manual(values = st2_color) +
    xlab(NULL)+
    theme_bw() +
    theme(text=element_text(family = "Helvetica", size=7),
          axis.text = element_text(family = "Helvetica", size=7),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x = element_text(size = 7, margin = margin(.1,0,.1,0, "cm")),
          legend.key.size = unit(.4, 'cm'),
          legend.text=element_text(size=8),
          legend.spacing.y = unit(.3, 'cm'),
          legend.margin=margin(0,0,0,0))+
    ylim(c(0, 0.2))+
    facet_wrap(~Patient)
ggsave(paste0(mplotPath, test_text, "epi_emt_boxplot.pdf"), g1, width = 3, height=1.8, units = "in", device = "pdf")

compare_groups <- list(
    c("O", "O(+M)"), 
    c("O", "O(+M+F)"), 
    c("O(+M)", "O(+M+F)"), 
    c("T", "T(+M)"), 
    c("T", "T(+M+F)"), 
    c("T(+M)", "T(+M+F)")
)
test_group = "HALLMARK_EMT"

pres_state <- sapply(levels(epi_proj$Patient), function(cur_patient) {
    print(cur_patient)
    cur_tbl = epi_proj[epi_proj$Patient == cur_patient,]
    pres_patient <- list()
    for(cur_group in compare_groups) {
        print(cur_group)
        pres_patient[[paste0(make.names(cur_group), collapse = "_vs_")]] <- t.test(cur_tbl[[test_group]][cur_tbl$SampleType2 == cur_group[1]], cur_tbl[[test_group]][cur_tbl$SampleType2 == cur_group[2]], alternative = "less")
    }
    sapply(pres_patient, function(x) x$p.value)
})

pres_state_star <- pres_state
pres_state_star[pres_state > 0.05] = "ns"
pres_state_star[pres_state <= 0.05] = "*"
pres_state_star[pres_state <= 0.01] = "**"
pres_state_star[pres_state <= 0.001] = "***"
pres_state_star[pres_state <= 0.0001] = "****"


write.xlsx(list(test_res = pres_state), paste0(msheetPath, test_text, cur_gene, "_emt_sig_ttest_res.xlsx"), rowNames = T)
write.xlsx(list(test_res = pres_state_star), paste0(msheetPath, test_text, cur_gene, "_emt_sig_ttest_res_star.xlsx"), rowNames = T)



