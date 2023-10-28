


library(VisCello)
library(ggrastr)

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



test_text <- "myeloid_hs_invivo_0413_"
cur_vis <- clist$`Myeloid cells in vivo [20210411]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [10PC, batch regressed]`, pData(eset)[cur_vis@idx,])
cur_proj <- cur_proj[!cur_proj$Cell_type_myeloid %in% c("Unannotated", "Macrophage-T cell doublet"), ]
mctype_color <- get_factor_color(unique(cur_proj$Cell_type_myeloid), "Set2")
names(mctype_color) <- unique(cur_proj$Cell_type_myeloid)
cur_proj$Cell_type_myeloid <- factor(cur_proj$Cell_type_myeloid, levels = c("Macrophage", "cDC1 (BATF3+)", "cDC2 (CD1C+)", "LAMP3+ DC (LAMP3+)", "pDC (LILRA4+)"))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_type_myeloid", pal=mctype_color, size = .3, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr=T) + 
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
ggsave(paste0(mplotPath, test_text, "myeloid_ctype.pdf"), g1, width = 4, height=2, units = "in", device = "pdf")

stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)
cur_proj$SampleType <- factor(cur_proj$SampleType, levels = names(stype_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="SampleType", pal=stype_color, size = .5, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5) + 
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

nbin <- 6
g1<-ggplot(cur_proj, aes_string("UMAP_1", "UMAP_2")) + 
    stat_density2d(aes(fill = SampleType, color = SampleType), geom="polygon", alpha = .2, bins = nbin) +
    # stat_density2d(data = cur_proj[cur_proj$SampleType == "Colon",], fill = stype_color["Colon"], color = stype_color["Colon"], geom="polygon", alpha = .1, bins = nbin) +
    # stat_density2d(data = cur_proj[cur_proj$SampleType == "Liver",], fill = stype_color["Liver"], color = stype_color["Liver"], geom="polygon", alpha = .1, bins = nbin) + 
    # stat_density2d(data = cur_proj[cur_proj$SampleType == "Tumor",], fill = stype_color["Tumor"], color = stype_color["Tumor"], geom="polygon", alpha = .1, bins = nbin) + 
    # stat_density2d(data = cur_proj[cur_proj$SampleType == "LiverMet",], fill = stype_color["LiverMet"], color = stype_color["LiverMet"], geom="polygon", alpha = .1, bins = nbin) + 
    xlim(c(-11,5)) + ylim(c(-8,6))+
    # scale_color_identity(name = "Sample type",
    #                      breaks = c(stype_color["Colon"], stype_color["Liver"], stype_color["Tumor"], stype_color["LiverMet"]),
    #                      labels = c("Colon", "Liver", "Tumor", "LiverMet"),
    #                      guide = "legend") + 
    scale_color_manual(values = stype_color)+
    scale_fill_manual(values = stype_color)+
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
ggsave(paste0(mplotPath, test_text, "sampletype_distribution.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")


write.xlsx(cur_proj[, c("UMAP_1", "UMAP_2", "Cell_type_myeloid", "SampleType")], paste0(msheetPath, test_text, "sampletype_distribution.xlsx"), rowNames = T)



cell_count<-as.data.frame.matrix(table(cur_proj$SampleType, cur_proj$Cell_type_myeloid))
# cell_frac <- as.data.frame(cell_count / rowSums(cell_count))
cell_count$Sample_type <- rownames(cell_count)
plot_df <- reshape2::melt(cell_count)
colnames(plot_df) <- c("Sample_type", "Cell_type", "Count")
ggplot(plot_df, aes(x = Sample_type, y = Count, fill = Cell_type)) + 
    geom_bar(stat = "identity")



#test_text <- "tam_hs_invivo_0413_"
test_text <- "tam_hs_invivo_0607_"
cur_vis <- clist$`Monocytes/macrophages in vivo [20210411]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [zoom]`, pData(eset)[cur_vis@idx,])

plot_g <- c("IL1B", "SPP1", "C1QC")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(plot_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- plot_g
glist <- list()
for(g in plot_g) {
    print(g)
    cur_proj$gene_expr <- g_exprs[[g]]
    ecut <- quantile(cur_proj$gene_expr[cur_proj$gene_expr > 0], .975)
    if(ecut < 1) ecut = 1
    cur_proj$gene_expr[cur_proj$gene_expr > ecut] = ecut
    glist[[g]]<-plotProj(cur_proj, dim_col = c(1,2), group.by="gene_expr", pal="viridis", size = .5, legend.title = g, layer_zero = T, rastr=F) + guides(color = guide_colorbar(barwidth = 10, barheight = 1, title = g, raster = T)) +
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
ggsave(paste0(mplotPath, test_text, "gexpr_umap.pdf"), g1, width = 3.3, height=2.3, units = "in", device = "pdf")


write.xlsx(cbind(cur_proj[, c("UMAP_1", "UMAP_2")], g_exprs), paste0(msheetPath, test_text, "gexpr_umap.xlsx"), rowNames = T)



#### Trajectory analysis using slingshot ####
library(slingshot)
# First define states using key genes
mac_state_g <-  c("IL1B", "SPP1", "C1QC")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(mac_state_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- mac_state_g
pass_qt <- apply(g_exprs, 2, function(x) {
    x <- scale(x, center = T, scale = T)
    ifelse(x > quantile(x[x > 0], .75), x, 0)
})
state_cell <- apply(pass_qt, 1, function(x) {
    if(all(x == 0)) {
        return(NA)
    } else {
        mac_state_g[which.max(x)[1]]
    }
})

cur_proj$Cell_state <- state_cell
state_color <- c("IL1B" = "#33a02c", "SPP1" = "#e31a1c", "C1QC" = "#ff7f00")
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_state", pal=state_color, size = 1, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5) + 
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
ggsave(paste0(mplotPath, test_text, "mac_state.pdf"), g1, width = 3.4, height=1.6, units = "in", device = "pdf")

eset$Macrophage_state <- cur_proj$Cell_state[match(colnames(eset), rownames(cur_proj))]
saveRDS(eset, paste0("../hcc_final_cello/eset.rds"))
# l2 <- slingshot(cur_proj[,c(1,2)],cur_proj$Cell_state, start.clus = "S100A8", end.clus=c("IL1B", "SPP1", "C1QC"), stretch=5)
# saveRDS(l2, paste0(mstatePath, test_text, "l2.rds"))
# 
# l2 <- readRDS( paste0(mstatePath, test_text, "l2.rds"))
# postscript(paste0(mplotPath, test_text, "slingshot_trajectory.eps"), width = 9, height = 4)
# plot(reducedDims(l2), col = state_color[cur_proj$Cell_state], pch=16, asp =1, cex = .5)
# lines(l2, lwd=2)
# dev.off()

# Plot heatmap and pathway enrichment
de_list <- read_excel_allsheets(paste0(msheetPath, "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))
gnum = 10
cnum = 100
plot_deg<- unname(unlist(lapply(de_list[names(state_color)], function(x) x$gene_short_name[1:gnum])))
plot_cell <- unlist(lapply(names(state_color), function(x) {sample(which(eset$Macrophage_state == x), cnum)}))
eset_plot <- eset[match(plot_deg, fData(eset)$gene_short_name),plot_cell]

#eset_plot <- eset_plot[rowMeans(exprs(eset_plot) > 0) > .2, ]
plot_meta <- pData(eset_plot)
plot_meta <- plot_meta[, c("Macrophage_state", "SampleType", "Dataset"), drop=F]
plot_meta$Macrophage_state <- factor(plot_meta$Macrophage_state, levels = names(state_color))
plot_meta <- plot_meta[order(plot_meta$Macrophage_state, plot_meta$SampleType, plot_meta$Dataset),]
non_na_cells_ordered <- rownames(plot_meta)
value <- eset_plot@assayData$norm_exprs[,non_na_cells_ordered]
rownames(value) <- fData(eset_plot)$gene_short_name[match(rownames(value), rownames(eset_plot))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

pdf(paste0(mplotPath, test_text, "state_deg", ".pdf"), width = 6, height=5)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = F, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta[,"Macrophage_state",drop=F], annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   annotation_colors = list("Macrophage_state" = state_color),
                   fontsize = 9)
dev.off()



gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
cur_eset <- eset[, rownames(cur_proj)]
cur_eset <- cur_eset[rowSums(exprs(cur_eset)) > 0,]
gene_background <- rownames(cur_eset)
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)

go_list <- lapply(de_list, function(x) {
    res_g <- x$gene_id
    include_g.df <- bitr(res_g, fromType = gene_id_type,
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
    return(go_res)
})
saveRDS(go_list, paste0(mstatePath, test_text,  "_macstate_go_res.rds"))

# GO plot
go_list <- readRDS(paste0(mstatePath, "tam_hs_invivo_0413_",  "_macstate_go_res.rds"))
plot_list <- list()
num_pathway <- 6
for(state in names(state_color)) {
    go_res <- go_list[[state]]
    include_go_simplified <- clusterProfiler::simplify(go_res, cutoff = .5)
    plot_df<- include_go_simplified@result[1:min(nrow(include_go_simplified@result),num_pathway),]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    plot_df$state <- state
    plot_list[[state]] <- plot_df 
}
plot_df <- do.call(rbind, plot_list)
plot_df$state <- factor(plot_df$state, levels = names(state_color))
g1 <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
    geom_bar(stat = "identity", aes(fill = -log10(qvalue))) + 
    coord_flip() +
    xlab("") + 
    theme(
        text=element_text(family = "Helvetica", size=8,color="black"),
        axis.text = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
    monocle:::monocle_theme_opts() + 
    facet_wrap(~state, nrow=1)

pdf(paste0(mplotPath, test_text, "_macstate_go_bar.pdf"), width = 6, height = 3)
g1
dev.off()

plot_mtx <- reshape2::dcast(plot_df[, c("state", "qvalue", "Description")], Description ~ state, value.var = "qvalue")
plot_mtx <- plot_mtx[order(plot_mtx$IL1B, plot_mtx$SPP1, plot_mtx$C1QC), ]
plot_mtx[is.na(plot_mtx)] = 1
rownames(plot_mtx) <- plot_mtx$Description
plot_mtx$Description <- NULL
plot_mtx <- -log10(plot_mtx)
plot_mtx[plot_mtx > 8] = 8
saveRDS(plot_mtx, paste0(mstatePath, test_text, "macstate_go_heatmap_mtx.rds"))
pdf(paste0(mplotPath, test_text, "_macstate_go_heatmap.pdf"), width = 3.35, height = 2.6)
pheatmap(plot_mtx, cluster_rows = F, cluster_cols = F, color = get_numeric_color("viridis"), fontsize = 8)
dev.off()


# Gene signature score
de_list <- read_excel_allsheets(paste0(msheetPath, "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))
max_gnum = 100
sig_deg <- lapply(de_list, function(x) {
    x$gene_short_name[1:min(max_gnum, nrow(x))]
})
lapply(sig_deg, length)

cur_vis <- clist$`Monocytes/macrophages in vivo [20210411]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [zoom]`, pData(eset)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Cluster

library(AUCell)
cur_eset<- eset[, rownames(cur_proj)]
cur_eset <- cur_eset[rowMeans(exprs(cur_eset) > 0) > .01,]
expr_hc <- as.matrix(exprs(cur_eset))
rownames(expr_hc) <- fData(cur_eset)$gene_short_name
cells_rankings_hc <- AUCell_buildRankings(expr_hc)
cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(sig_deg, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "aucell_reactome_0.2.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "aucell_reactome_0.2.rds"))
matrixAUC <- getAUC(cells_AUC)
plot_meta <-  cur_proj

glist <- list()
#sig_deg <- sig_deg[c("IL1B", "SPP1", "C1QC")]
for(show_p in c("IL1B", "SPP1", "C1QC")) {
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    cut_qt = .975
    plot_df$activity_score[plot_df$activity_score > quantile(plot_df$activity_score, cut_qt)] <-  quantile(plot_df$activity_score, cut_qt)
    glist[[show_p]] <- ggplot(plot_df, aes_string("UMAP_1", "UMAP_2")) +
        geom_point(aes(color = activity_score),size = .5, stroke = 0) + 
        scale_color_gradientn(colors = get_numeric_color("BlueGreenRed")) + 
        guides(color=guide_colorbar(barwidth = 4, barheight = .5, title=paste0(show_p, " signature score")))+
        theme(text=element_text(family = "Helvetica", size=8),
              legend.position="top",
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              legend.text=element_text(size=7),
              axis.text = element_text(size=8),
              legend.margin=margin(0,0,0,0))+
        ggtitle(NULL) +
        monocle:::monocle_theme_opts()
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "activity_umap", ".pdf"), g1, width = 6, height=2, units = "in", device = "pdf")


plot_proj <- cbind(plot_meta[, "SampleType",drop=F], t(matrixAUC[c("IL1B", "SPP1", "C1QC"), rownames(plot_meta), drop=F]))
write.xlsx(plot_proj, paste0(msheetPath, test_text, "sig_boxplot.xlsx"), rowNames=T)

plot_proj <- reshape2::melt(plot_proj)
colnames(plot_proj) <- c("SampleType", "State",  "Signature")

stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)

plot_proj$SampleType <- factor(plot_proj$SampleType, levels = names(stype_color))
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Signature", fill = "SampleType")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = stype_color)+
    theme_bw()+
    theme(text=element_text(size=9),
          axis.text = element_text(size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_sig_boxplot.pdf")),  width = 3.5, height=1.5)
g1
dev.off()

# t test
pres_all <- sapply(c("IL1B", "SPP1", "C1QC"), function(state){
    cur_tbl = plot_proj[plot_proj$State == state,]
    res_tumor <- t.test(cur_tbl$Signature[cur_tbl$SampleType == "Tumor"], cur_tbl$Signature[cur_tbl$SampleType == "Colon"])
    res_liver <- t.test(cur_tbl$Signature[cur_tbl$SampleType == "Liver"], cur_tbl$Signature[cur_tbl$SampleType == "Colon"])
    res_mets <- t.test(cur_tbl$Signature[cur_tbl$SampleType == "LiverMet"], cur_tbl$Signature[cur_tbl$SampleType == "Colon"])
    pres<-data.frame(tumor = res_tumor$p.value, liver = res_liver$p.value, mets = res_mets$p.value)
    rownames(pres) <- state
    return(pres)
})
pres_all
pres_all[pres_all <= 0.0001] = "****"
pres_all[pres_all <= 0.001] = "***"
pres_all[pres_all <= 0.01] = "**"
pres_all[pres_all <= 0.05] = "*"
pres_all[pres_all > 0.05] = "ns"
pres_all

# Stratify by tumor stage
patient_meta <- read_excel("C:/Users/qinzh/Dropbox/ColonManuscript/patient_meta.xlsx")
plot_proj <- cbind(plot_meta[, c("SampleType", "Patient"),drop=F], t(matrixAUC[c("IL1B", "SPP1", "C1QC"), rownames(plot_meta), drop=F]))
msi_patient <- patient_meta$Patient[patient_meta$MSI == "MSI"]
plot_proj <- plot_proj[plot_proj$SampleType == "Tumor" & !plot_proj$Patient %in% c("40", "86", "87"),]
plot_proj$Stage = patient_meta$Stage[match(as.numeric(plot_proj$Patient), patient_meta$Patient)]
plot_proj <- reshape2::melt(plot_proj[c("IL1B", "SPP1", "C1QC", "Stage")], id.var="Stage")
colnames(plot_proj) <- c("Stage", "State",  "Signature")
plot_proj$Stage <- factor(plot_proj$Stage)
use_color = get_factor_color(levels(plot_proj$Stage), "BuGn")
names(use_color) = levels(plot_proj$Stage)
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Signature", fill = "Stage")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = use_color)+
    theme_bw()+
    theme(text=element_text(size=9),
          axis.text = element_text(size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_signature_stage_boxplot.pdf")),  width = 3, height=1.75)
g1
dev.off()

# FateID
library(FateID)
state_color <- c("IL1B" = "#33a02c", "SPP1" = "#e31a1c", "C1QC" = "#ff7f00")
cur_vis <- clist$`Monocytes/macrophages in vivo [20210411]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [zoom]`, pData(eset)[cur_vis@idx,])
maxg <- 300
use_deg <- unlist(lapply(de_list, function(x) {
    x$gene_id[1:min(maxg, nrow(x))]
}))
cur_eset$Macrophage_state[cur_eset$Macrophage_state == "S100A12"] = "IL1B"
cur_eset <- eset[use_deg, rownames(cur_proj)]
x = as.matrix(cur_eset@assayData$norm_exprs)[use_deg,]
y = as.numeric(factor(cur_eset$Macrophage_state, levels = names(state_color)))
y[is.na(y)] <- 4
tar = c(1, 2,3)
fb <- fateBias(x, y, tar, z=NULL, minnr=5, minnrh=10, adapt=TRUE, confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)
saveRDS(fb, paste0(mstatePath, test_text, "fate_id.rds"))

fb <- readRDS(paste0(mstatePath, "myeloid_cell_invivo_20210411_",  "fate_id.rds"))
# dr <- list()
# dr$umap <- list()
# dr$umap$D2=cur_proj[,c(1,2)]
# plotFateMap(y,dr,k=2,m="umap",fb=fb,g="t2")

fb_prob <- fb$probs
colnames(fb_prob) <- paste0("Prob_", names(state_color))
cur_proj <- cbind(cur_proj, fb_prob[rownames(cur_proj), ])
saveRDS(cur_proj, paste0(mstatePath, test_text, "cur_proj.rds"))

cur_proj <- readRDS(paste0(mstatePath, "myeloid_cell_invivo_20210411_", "cur_proj.rds"))
target_groups <- c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC")
tg_color <- c("Prob_IL1B"="#33a02c", "Prob_SPP1"="#e31a1c", "Prob_C1QC"="#ff7f00")
glist <- list()
for(target_group in target_groups) {
    use_color <- colorRampPalette(c("white", tg_color[target_group]))(10)
    glist[[target_group]]<-plotProj(cur_proj, dim_col = c(1,2), group.by=target_group, pal=use_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "top", legend.title = NULL, keyheight=.5, legend_barwidth = 3, legend_barheight = .5, rastr = T) + 
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
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "fate_umap.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")


write.xlsx(cur_proj[, c("UMAP_1", "UMAP_2", "Prob_IL1B", "Prob_SPP1", "Prob_C1QC")], paste0(msheetPath, test_text, "fate_umap.xlsx"), rowNames = T)


fb_prob = cur_proj[match(colnames(eset), rownames(cur_proj)),target_groups]
pData(eset) <- cbind(pData(eset), fb_prob)
saveRDS(eset, paste0("../hcc_final_cello/eset.rds"))


plot_proj <- cur_proj[, c("SampleType", "Prob_IL1B", "Prob_SPP1", "Prob_C1QC")]
plot_proj <- reshape2::melt(plot_proj)
colnames(plot_proj) <- c("SampleType", "State",  "Probability")

use_color <- c("Prob_IL1B" = "#33a02c", "Prob_SPP1" = "#e31a1c", "Prob_C1QC" = "#ff7f00")
stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)

plot_proj$SampleType <- factor(plot_proj$SampleType, levels = names(stype_color))
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Probability", fill = "SampleType")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = stype_color)+
    theme(text=element_text(family = "Helvetica", size=8),
        axis.text = element_text(family = "Helvetica", size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)+
    theme_bw()
pdf(paste0(paste0(mplotPath, test_text, "_fate_boxplot.pdf")),  width = 3.3, height=1.3)
g1
dev.off()

# pres_all <- sapply(c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC"), function(state){
#     cur_tbl = plot_proj[plot_proj$State == state,]
#     res_tumor <- t.test(cur_tbl$Probability[cur_tbl$SampleType == "Tumor"], cur_tbl$Probability[cur_tbl$SampleType == "Colon"])
#     res_liver <- t.test(cur_tbl$Probability[cur_tbl$SampleType == "Liver"], cur_tbl$Probability[cur_tbl$SampleType == "Colon"])
#     res_mets <- t.test(cur_tbl$Probability[cur_tbl$SampleType == "LiverMet"], cur_tbl$Probability[cur_tbl$SampleType == "Colon"])
#     pres<-data.frame(tumor = res_tumor$p.value, liver = res_liver$p.value, mets = res_mets$p.value)
#     rownames(pres) <- state
#     return(pres)
# })
pres_all <- sapply(c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC"), function(state){
    cur_tbl = plot_proj[plot_proj$State == state,]
    res_tumor <- tryCatch(t.test(cur_tbl$Probability[cur_tbl$SampleType == "Tumor"], cur_tbl$Probability[cur_tbl$SampleType == "Colon"]), error = function(err)return(data.frame(p.value = NA)))
    res_mets <- tryCatch(t.test(cur_tbl$Probability[cur_tbl$SampleType == "LiverMet"], cur_tbl$Probability[cur_tbl$SampleType == "Liver"]), error = function(err)return(data.frame(p.value = NA)))
    pres<-data.frame(tumor = res_tumor$p.value, mets = res_mets$p.value)
    rownames(pres) <- state
    return(pres)
})
pres_all

write.xlsx(as.data.frame.matrix(table(plot_proj$SampleType, plot_proj$State)),  paste0(msheetPath, test_text, "mac_sptype_cellcount.xlsx"), rowNames=T)
write.xlsx(as.data.frame(pres_all), paste0(msheetPath, test_text, "mac_sptype_pres_all.xlsx"), rowNames=T)
pres_all[pres_all <= 0.0001] = "****"
pres_all[pres_all <= 0.001] = "***"
pres_all[pres_all <= 0.01] = "**"
pres_all[pres_all <= 0.05] = "*"
pres_all[pres_all > 0.05] = "ns"
pres_all

# Stratify by tumor stage
patient_meta <- read_excel("~/Dropbox/ColonManuscript/patient_meta.xlsx")
plot_proj <- cur_proj[cur_proj$SampleType == "Tumor" & !cur_proj$Patient %in% c("40", "86", "87"),]
plot_proj$Stage = patient_meta$Stage[match(as.numeric(plot_proj$Patient), patient_meta$Patient)]
plot_proj <- reshape2::melt(plot_proj[c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC", "Stage")], id.var="Stage")
colnames(plot_proj) <- c("Stage", "State",  "Probability")
plot_proj$Stage <- factor(plot_proj$Stage)
use_color = get_factor_color(levels(plot_proj$Stage), "BuGn")
names(use_color) = levels(plot_proj$Stage)
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Probability", fill = "Stage")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = use_color)+
    theme_bw()+
    theme(text=element_text(size=9),
          axis.text = element_text(size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_prob_stage_boxplot.pdf")),  width = 2.7, height=1.3)
g1
dev.off()


pres_all <- sapply(c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC"), function(state){
    cur_tbl = plot_proj[plot_proj$State == state,]
    res_2 <- t.test(cur_tbl$Probability[cur_tbl$Stage == 2], cur_tbl$Probability[cur_tbl$Stage == 1])
    res_3 <- t.test(cur_tbl$Probability[cur_tbl$Stage == 3], cur_tbl$Probability[cur_tbl$Stage == 1])
    res_4 <- t.test(cur_tbl$Probability[cur_tbl$Stage == 4], cur_tbl$Probability[cur_tbl$Stage == 1])
    pres<-data.frame(II = res_2$p.value, III = res_3$p.value, IV = res_4$p.value)
    rownames(pres) <- state
    return(pres)
})
pres_all
write.xlsx(as.data.frame(pres_all), paste0(msheetPath, test_text, "mac_stage_type_pres_all.xlsx"), rowNames=T)

pres_all[pres_all <= 0.0001] = "****"
pres_all[pres_all <= 0.001] = "***"
pres_all[pres_all <= 0.01] = "**"
pres_all[pres_all <= 0.05] = "*"
pres_all[pres_all > 0.05] = "ns"
pres_all

# Add normal 
plot_proj <- cur_proj[cur_proj$SampleType %in% c("Colon", "Tumor") & !cur_proj$Patient %in% c("40", "86", "87"),]
plot_proj$Stage = patient_meta$Stage[match(as.numeric(plot_proj$Patient), patient_meta$Patient)]
plot_proj$Stage[plot_proj$SampleType == "Colon"] = "Colon"
plot_proj$Stage <- factor(plot_proj$Stage, levels = c("Colon", 1:4))


plot_proj <- plot_proj[c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC", "Stage")]
colnames(plot_proj) <- c("IL1B", "SPP1", "C1QC", "Stage")

plot_proj <- reshape2::melt(plot_proj, id.var="Stage")
colnames(plot_proj) <- c("Stage", "State",  "Probability")
plot_proj$Stage <- factor(plot_proj$Stage)
use_color = get_factor_color(levels(plot_proj$Stage), "BuGn")
names(use_color) = levels(plot_proj$Stage)
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Probability", fill = "Stage")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = use_color)+
    theme_bw()+
    theme(text=element_text(size=9),
          axis.text = element_text(size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_prob_stage_wt_normal_boxplot.pdf")),  width = 3, height=1.3)
g1
dev.off()

write.xlsx(as.data.frame.matrix(table(plot_proj$Stage, plot_proj$State)),  paste0(msheetPath, test_text, "mac_stage_type_cellcount.xlsx"), rowNames=T)



# Stratify by MSI
patient_meta <- read_excel("~/Dropbox/ColonManuscript/patient_meta.xlsx")
save_proj <- cur_proj[!cur_proj$Patient %in% c("40"),]
save_proj$Stage = patient_meta$Stage[match(as.numeric(save_proj$Patient), patient_meta$Patient)]
save_proj$MSI = patient_meta$MSI[match(as.numeric(save_proj$Patient), patient_meta$Patient)]
write.xlsx(save_proj[,c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC", "SampleType", "Stage", "MSI")], paste0(msheetPath, test_text, "prob_stage_wt_normal_boxplot.xlsx"), rowNames=T)


plot_proj <- cur_proj[cur_proj$SampleType == "Tumor" & !cur_proj$Patient %in% c("40"),]
plot_proj$Stage = patient_meta$Stage[match(as.numeric(plot_proj$Patient), patient_meta$Patient)]
plot_proj$MSI = patient_meta$MSI[match(as.numeric(plot_proj$Patient), patient_meta$Patient)]

plot_proj = plot_proj[plot_proj$Stage == 3,]
plot_proj <- reshape2::melt(plot_proj[c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC", "MSI")], id.var="MSI")
colnames(plot_proj) <- c("MSI", "State",  "Probability")
# use_color = get_factor_color(levels(plot_proj$Stage), "BuGn")
# names(use_color) = levels(plot_proj$Stage)
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Probability", fill = "MSI")) +
    geom_boxplot(alpha = .5, outlier.colour = NA)+
    geom_point(position = position_jitterdodge(jitter.width = 0.2),
               aes_string(color = "MSI", group="MSI"), size = .5, stroke = 0)+
    #scale_fill_manual(values = use_color)+
    theme_bw()+
    theme(text=element_text(size=8),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_prob_msi_boxplot.pdf")),  width = 2.7, height=1.3)
g1
dev.off()

pres_all <- sapply(c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC"), function(state){
    cur_tbl = plot_proj[plot_proj$State == state,]
    res <- t.test(cur_tbl$Probability[cur_tbl$MSI == "MSI"], cur_tbl$Probability[cur_tbl$MSI == "MSS"])
    pres<-data.frame(res = res$p.value)
    rownames(pres) <- state
    return(pres)
})
pres_all
write.xlsx(as.data.frame(pres_all), paste0(msheetPath, test_text, "mac_msistate_type_pres_all.xlsx"), rowNames=T)
write.xlsx(as.data.frame.matrix(table(plot_proj$MSI, plot_proj$State)),  paste0(msheetPath, test_text, "mac_msistate_type_cellcount.xlsx"), rowNames=T)

pres_all[pres_all <= 0.0001] = "****"
pres_all[pres_all <= 0.001] = "***"
pres_all[pres_all <= 0.01] = "**"
pres_all[pres_all <= 0.05] = "*"
pres_all[pres_all > 0.05] = "ns"
pres_all







# Association with T activity
t_vis <- clist$`T cells in vivo [20221224, IFF]`
t_proj <- cbind(t_vis@proj$`UMAP-2D [30PC, batch regressed]`, pData(eset)[t_vis@idx,])
t_count<-as.data.frame.matrix(table(t_proj$Dataset, t_proj$Cell_type_T))
t_frac <- as.data.frame(t_count / rowSums(t_count))

m_vis <- clist$`Monocytes/macrophages in vivo [20210411]`
m_proj <- cbind(m_vis@proj$`UMAP-2D [zoom]`, pData(eset)[m_vis@idx,])

m_score = as.data.frame(as.data.table(m_proj)[ , lapply(.SD, mean, na.rm=TRUE), by=Dataset, .SDcols=c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC")])
rownames(m_score) = m_score$Dataset
m_score$Dataset <- NULL
shared_dset = intersect(rownames(t_frac), rownames(m_score))
#shared_dset = grep("Tumor_", shared_dset, value= T)
t_frac_ordered = as.data.frame(t_frac[shared_dset,])
m_score_ordered = m_score[shared_dset,]

cor_dset = grep("Tumor_", shared_dset, value= T)
show_ctype = colnames(t_frac_ordered)[!colnames(t_frac_ordered) %in% c("T-Carcinoma cell doublet")]
cor_mtx = cor(t_frac_ordered[cor_dset,show_ctype], m_score_ordered[cor_dset,])
#cor_mtx = cor(t_frac_ordered, m_score_ordered)
#pdf(paste0(mplotPath, test_text, "tm_analysis_cor_heatmap_full.pdf"), width = 3.7, height =4)
breaksList <- seq(-.7,.7,by=.1)
pdf(paste0(mplotPath, test_text, "tm_analysis_cor_heatmap.pdf"), width = 3.3, height =2)
pheatmap(cor_mtx, cluster_rows = T, cluster_cols = F, show_colnames = T, show_rownames = T,
         clustering_method = "ward.D2",
         display_numbers =T, number_format ="%.2f",fontsize_number=6,fontsize = 8,
         breaks = breaksList,
         get_numeric_bin_color(bins = 1:length(breaksList), "RdBu"))
dev.off()


write.xlsx(t_frac_ordered, paste0(msheetPath, "tm_analysis_cor_heatmap_t_frac_ordered.xlsx"), rowNames=T)
write.xlsx(m_score_ordered, paste0(msheetPath, "tm_analysis_cor_heatmap_m_score_ordered.xlsx"), rowNames=T)





cor_res<-psych::corr.test(t_frac_ordered[cor_dset,show_ctype], m_score_ordered[cor_dset,], adjust = "none")
cor_res$stars
cor_res$p


# Plot individual correlation
tm_cbn = cbind(t_frac[shared_dset,], m_score[shared_dset,])
tm_cbn$SampleType = sapply(strsplit(rownames(tm_cbn), "_"), function(x)x[1])
stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)
tm_cbn$SampleType <- factor(tm_cbn$SampleType, levels = names(stype_color))
library(ggpmisc)

plot_regression <- function(data, x, y, point_size = 1, label.y.dist = .2, label.size = 2, ylab = TRUE,
                               xlim = c(0,.7), ylim = c(0,.7), plot_pvalue = FALSE, color_legend = F) {
    g1 <- ggplot(data, aes(x=.data[[x]], y=.data[[y]])) +
        geom_point(aes(color = SampleType), size = point_size, stroke = 0) +
        scale_color_manual(values = stype_color)+
        stat_smooth(color = "black", method = "lm", alpha = .1) +

        stat_smooth(data = data[data$SampleType == "Tumor",], method = "lm",
                    color = "#e31a1c", fill= "#e31a1c", alpha = .1) +
        xlim(xlim) + 
        ylim(ylim) + 
        theme_bw() + 
        theme(text=element_text(size=8),
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              legend.margin=margin(0,0,0,0),
              #legend.box.margin=margin(-10,-10,-10,-10),
              #plot.margin = unit(c(0,0.4,0,0.1), "cm")
              )
    if(plot_pvalue) {
        g1 = g1 + 
            stat_poly_eq(formula = y ~ x,
                         eq.with.lhs = "italic(hat(y))~`=`~",
                         aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label..,sep = "~~~")),
                         #aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE, size = label.size,
                         coef.digits = 2, rr.digits = 2, p.digits = 3)+
            stat_poly_eq(data = data[data$SampleType == "Tumor",], 
                         formula = y ~ x,
                         eq.with.lhs = "italic(hat(y))~`=`~",
                         aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
                         #aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                         parse = TRUE, size = label.size, 
                         color = "#e31a1c", 
                         label.y = 1-label.y.dist,
                         coef.digits = 2, rr.digits = 2, p.digits = 3)
    } else {
        g1 = g1 + 
            stat_poly_eq(formula = y ~ x,
                         eq.with.lhs = "italic(hat(y))~`=`~",
                         #aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label..,sep = "~~~")),
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE, size = label.size,
                         coef.digits = 2, rr.digits = 2, p.digits = 3)+
            stat_poly_eq(data = data[data$SampleType == "Tumor",], 
                         formula = y ~ x,
                         eq.with.lhs = "italic(hat(y))~`=`~",
                         #aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                         parse = TRUE, size = label.size, 
                         color = "#e31a1c", 
                         label.y = 1-label.y.dist,
                         coef.digits = 2, rr.digits = 2, p.digits = 3)
    }
    if(!ylab) {
        g1 = g1 + ylab(NULL) 
    }
    if(!color_legend) {
        g1 = g1 + guides(color = F) 
    }
    return(g1)
}

pdf(paste0(mplotPath, test_text, "tm_analysis_regression", ".pdf"), width = 3.4, height=1.65)
for(ctype in colnames(t_frac_ordered)) {
    glist = list(
        plot_regression(tm_cbn, x = "Prob_SPP1", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 2),
        plot_regression(tm_cbn, x = "Prob_C1QC", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 2, ylab= FALSE)
    )
    g1<-do.call(grid.arrange,c(glist, nrow = 1))
    print(g1)
    #ggsave(paste0(mplotPath, test_text, "tm_analysis_", make.names(ctype), ".pdf"), g1, width = 4, height=2, units = "in", device = "pdf")
}
dev.off()

pdf(paste0(mplotPath, test_text, "tm_analysis_regression_wtp", ".pdf"), width = 3.4, height=1.65)
for(ctype in colnames(t_frac_ordered)) {
    glist = list(
        plot_regression(tm_cbn, x = "Prob_SPP1", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 1.2, plot_pvalue = T),
        plot_regression(tm_cbn, x = "Prob_C1QC", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 1.2, ylab= FALSE, plot_pvalue = T)
    )
    g1<-do.call(grid.arrange,c(glist, nrow = 1))
    print(g1)
    #ggsave(paste0(mplotPath, test_text, "bm_analysis_", make.names(ctype), ".pdf"), g1, width = 4, height=2, units = "in", device = "pdf")
}
dev.off()

ctype = "GZMA/GZMB+ Teff/Tem"
df = data.frame(x = tm_cbn$Prob_C1QC, y=tm_cbn[[ctype]])
lm(y~x, df)
summary(lm(y~x, df))
pdf(paste0(mplotPath, test_text, "tm_analysis_regression_wt_legend", ".pdf"), width = 2, height=1.65)
plot_regression(tm_cbn, x = "Prob_SPP1", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 1.2, plot_pvalue = T, color_legend = T)+ theme(legend.position="top")
dev.off()



# Association with B activity
b_vis <- clist$`B cells in vivo [20221224, IFF]`
b_proj <- cbind(b_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset)[b_vis@idx,])
b_count<-as.data.frame.matrix(table(b_proj$Dataset, b_proj$Cell_type_B))
b_frac <- as.data.frame(b_count / rowSums(b_count))

m_vis <- clist$`Monocytes/macrophages in vivo [20210411]`
m_proj <- cbind(m_vis@proj$`UMAP-2D [zoom]`, pData(eset)[m_vis@idx,])

m_score = as.data.frame(as.data.table(m_proj)[ , lapply(.SD, mean, na.rm=TRUE), by=Dataset, .SDcols=c("Prob_IL1B", "Prob_SPP1", "Prob_C1QC")])
rownames(m_score) = m_score$Dataset
m_score$Dataset <- NULL
shared_dset = intersect(rownames(b_frac), rownames(m_score))
#shared_dset = grep("Tumor_", shared_dset, value= T)
b_frac_ordered = as.data.frame(b_frac[shared_dset,])
m_score_ordered = m_score[shared_dset,]

cor_dset = grep("Tumor_", shared_dset, value= T)
show_ctype = colnames(b_frac_ordered)[!colnames(b_frac_ordered) %in% c(#"IGHA+ Plasma cell (Tumor)", 
                                                                       "IGHA+ Plasma cell (Normal)")]
cor_mtx = cor(b_frac_ordered[cor_dset,show_ctype], m_score_ordered[cor_dset,])
breaksList <- seq(-.7,.7,by=.1)
#cor_mtx = cor(b_frac_ordered, m_score_ordered)
#pdf(paste0(mplotPath, test_text, "bm_analysis_cor_heatmap_full.pdf"), width = 3.7, height =4)
pdf(paste0(mplotPath, test_text, "bm_analysis_cor_heatmap.pdf"), width = 3.5, height =1.8)
pheatmap(cor_mtx, cluster_rows = T, cluster_cols = F, show_colnames = T, show_rownames = T,
         clustering_method = "ward.D2",
         display_numbers =T, number_format ="%.2f",fontsize_number=6,fontsize = 8,
         breaks = breaksList,
         get_numeric_bin_color(bins = 1:length(breaksList), "RdBu"))
dev.off()

cor_res<-psych::corr.test(b_frac_ordered[cor_dset,show_ctype], m_score_ordered[cor_dset,], adjust = "none")
cor_res$stars
cor_res$p

write.xlsx(b_frac_ordered, paste0(msheetPath, "bm_analysis_cor_heatmap_b_frac_ordered.xlsx"), rowNames=T)
write.xlsx(m_score_ordered, paste0(msheetPath, "bm_analysis_cor_heatmap_m_score_ordered.xlsx"), rowNames=T)


# Plot individual correlation
bm_cbn = cbind(b_frac[shared_dset,], m_score[shared_dset,])
bm_cbn$SampleType = sapply(strsplit(rownames(bm_cbn), "_"), function(x)x[1])
stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)
bm_cbn$SampleType <- factor(bm_cbn$SampleType, levels = names(stype_color))
library(ggpmisc)

pdf(paste0(mplotPath, test_text, "bm_analysis_regression", ".pdf"), width = 3.4, height=1.65)
for(ctype in colnames(b_frac_ordered)) {
    glist = list(
        plot_regression(bm_cbn, x = "Prob_SPP1", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 2),
        plot_regression(bm_cbn, x = "Prob_C1QC", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 2, ylab= FALSE)
    )
    g1<-do.call(grid.arrange,c(glist, nrow = 1))
    print(g1)
    #ggsave(paste0(mplotPath, test_text, "bm_analysis_", make.names(ctype), ".pdf"), g1, width = 4, height=2, units = "in", device = "pdf")
}
dev.off()

pdf(paste0(mplotPath, test_text, "bm_analysis_regression_wtp", ".pdf"), width = 3.4, height=1.65)
for(ctype in colnames(b_frac_ordered)) {
    glist = list(
        plot_regression(bm_cbn, x = "Prob_SPP1", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 1.2, plot_pvalue = T),
        plot_regression(bm_cbn, x = "Prob_C1QC", y = ctype, point_size = 1.5, label.y.dist = .14, label.size = 1.2, ylab= FALSE, plot_pvalue = T)
    )
    g1<-do.call(grid.arrange,c(glist, nrow = 1))
    print(g1)
    #ggsave(paste0(mplotPath, test_text, "bm_analysis_", make.names(ctype), ".pdf"), g1, width = 4, height=2, units = "in", device = "pdf")
}
dev.off()





test_text <- "mono_mm_invivo_"
eset_mouse <- readRDS("../mouse_colon_cello/eset.rds")
clist_mouse <- readRDS("../mouse_colon_cello/clist.rds")
cur_vis <- clist_mouse$`Mouse myeloid`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D Zoom`, pData(eset_mouse)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Cluster
cur_proj <- cur_proj[!cur_proj$Cluster %in% c(6,8,9),]
cur_proj$Dataset <- as.character(cur_proj$Dataset)
cur_proj$Dataset[grepl("Colon",cur_proj$Dataset)] <- "Colon"
use_color <- c(
    "Colon" = "#6a3d9a",   
    "Tumor AOM" = "#c51b7d",
    "Tumor ApcMin" = "#7fbc41",
    "Tumor APKS" = "#fdae61"
)
cur_proj$Dataset <- factor(cur_proj$Dataset, levels = names(use_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Dataset", pal=use_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5) + 
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
ggsave(paste0(mplotPath, test_text, "dataset.pdf"), g1, width = 3.6, height=1.6, units = "in", device = "pdf")


plot_g <- c("S100a8", "Spp1", "C1qc")
g_exprs <- as.data.frame(t(as.matrix(eset_mouse@assayData$norm_exprs[match(plot_g,fData(eset_mouse)$gene_short_name), match(rownames(cur_proj), colnames(eset_mouse))])))
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



# Compute signature score for human-derived macrophage substates signature on mouse macrophage
de_list <- read_excel_allsheets(paste0("~/Dropbox/ColonManuscript/sheets/", "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))
max_gnum = 100
sig_deg <- lapply(de_list, function(x) {
    x$gene_short_name[1:min(max_gnum, nrow(x))]
})
sig_deg_mouse <- lapply(sig_deg, function(x) {
    y <- mouse_to_human_symbol(gene_symbol = x, in.type = "hs",HMD_HumanPhenotype = HMD_HumanPhenotype)
    return(y[!is.na(y)])
})
lapply(sig_deg_mouse, length)

cur_vis <- clist_mouse$`Mouse myeloid`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D Zoom`, pData(eset_mouse)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Cluster
cur_proj <- cur_proj[!cur_proj$Cluster %in% c(6,8,9),]
cur_proj$Dataset <- as.character(cur_proj$Dataset)
cur_proj$Dataset[grepl("Colon",cur_proj$Dataset)] <- "Colon"
library(AUCell)
cur_eset<- eset[, rownames(cur_proj)]
cur_eset <- cur_eset[rowMeans(exprs(cur_eset) > 0) > .01,]
expr_hc <- as.matrix(exprs(cur_eset))
rownames(expr_hc) <- fData(cur_eset)$gene_short_name
cells_rankings_hc <- AUCell_buildRankings(expr_hc)
cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(sig_deg_mouse, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
#saveRDS(cells_AUC, paste0(rdsPath, test_text, cur_thresh, "cells_AUC", "_wt_ifg_cut75.rds"))
saveRDS(cells_AUC,paste0(mstatePath, test_text, "aucell_reactome_0.2.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "aucell_reactome_0.2.rds"))
matrixAUC <- getAUC(cells_AUC)
plot_meta <-  cur_proj
saveRDS(plot_meta, paste0("../mouse_colon_cello/rds/", test_text, "plot_meta.rds"))

plot_meta <- readRDS(paste0("../mouse_colon_cello/rds/", "mono_mm_invivo_", "plot_meta.rds"))
glist <- list()
sig_deg_mouse <- sig_deg_mouse[c("IL1B", "SPP1", "C1QC")]
for(show_p in names(sig_deg_mouse)) {
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    cut_qt = .975
    plot_df$activity_score[plot_df$activity_score > quantile(plot_df$activity_score, cut_qt)] <-  quantile(plot_df$activity_score, cut_qt)
    glist[[show_p]] <- ggplot(plot_df, aes_string("UMAP_1", "UMAP_2")) +
        geom_point(aes(color = activity_score),size = .5, stroke = 0) + 
        scale_color_gradientn(colors = get_numeric_color("BlueGreenRed")) + 
        guides(color=guide_colorbar(barwidth = 4, barheight = .5, title=paste0(show_p, " signature score")))+
        theme(text=element_text(family = "Helvetica", size=8),
              legend.position="top",
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              legend.text=element_text(size=7),
              axis.text = element_text(size=8),
              legend.margin=margin(0,0,0,0))+
        ggtitle(NULL) +
        monocle:::monocle_theme_opts()
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "activity_umap", ".pdf"), g1, width = 6, height=2, units = "in", device = "pdf")


plot_proj <- cbind(plot_meta[, "Dataset",drop=F], t(matrixAUC[c("IL1B", "SPP1", "C1QC"), rownames(plot_meta), drop=F]))
write.xlsx(plot_proj, paste0(msheetPath, test_text, "sig_boxplot.xlsx"), rowNames = T)

plot_proj <- reshape2::melt(plot_proj)
colnames(plot_proj) <- c("Dataset", "State",  "Probability")

use_color <- c("IL1B" = "#33a02c", "SPP1" = "#e31a1c", "C1QC" = "#ff7f00")
dataset_color <- c(
    "Colon" = "#6a3d9a",   
    "Tumor AOM" = "#c51b7d",
    "Tumor ApcMin" = "#7fbc41",
    "Tumor APKS" = "#fdae61"
)
plot_proj$Dataset <- factor(plot_proj$Dataset, levels = names(dataset_color))
g1 <- ggplot(plot_proj, aes_string(x = "State", y = "Probability", fill = "Dataset")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = dataset_color)+
    theme_bw()+
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_sig_boxplot.pdf")),  width = 3.5, height=1.5)
g1
dev.off()



pres_all <- sapply(c("IL1B", "SPP1", "C1QC"), function(state){
    cur_tbl = plot_proj[plot_proj$State == state,]
    res_AOM <- t.test(cur_tbl$Probability[cur_tbl$Dataset == "Tumor AOM"], cur_tbl$Probability[cur_tbl$Dataset == "Colon"])
    res_ApcMin <- t.test(cur_tbl$Probability[cur_tbl$Dataset == "Tumor ApcMin"], cur_tbl$Probability[cur_tbl$Dataset == "Colon"])
    res_APKS <- t.test(cur_tbl$Probability[cur_tbl$Dataset == "Tumor APKS"], cur_tbl$Probability[cur_tbl$Dataset == "Colon"])
    pres<-data.frame(AOM = res_AOM$p.value, ApcMin = res_ApcMin$p.value, APKS = res_APKS$p.value)
    rownames(pres) <- state
    return(pres)
})
pres_all
write.xlsx(as.data.frame.matrix(table(plot_proj$Dataset, plot_proj$State)),  paste0(msheetPath, test_text, "mm_sptype_cellcount.xlsx"), rowNames=T)
write.xlsx(as.data.frame(pres_all), paste0(msheetPath, test_text, "mm_sptype_pres_all.xlsx"), rowNames=T)

pres_all[pres_all <= 0.0001] = "****"
pres_all[pres_all <= 0.001] = "***"
pres_all[pres_all <= 0.01] = "**"
pres_all[pres_all <= 0.05] = "*"
pres_all[pres_all > 0.05] = "ns"
pres_all


