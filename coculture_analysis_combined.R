

library(VisCello)
mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

savePath <- "../hcc_final_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

scriptPath <- paste0("../hcc_cello/scripts/")
source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")
library(ggrastr)

test_text <- "coculture_combine_"
cur_vis <- clist$`Coculture experiment final set [IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Cluster
mapping_name <- c(
    "Tumoroid_28_24hr" = "T28_24hr",
    "Tumoroid_28_48hr" = "T28_48hr",
    "Tumoroid+Macrophage_28_24hr" = "T28+M_24hr",
    "Tumoroid+Macrophage_28_48hr" = "T28+M_48hr",
    "Macrophage_24hr" = "M(28)_24hr", 
    "Macrophage_48hr" = "M(28)_48hr", 
    
    "Tumoroid_86_36hr" = "T86_36hr",
    "Tumoroid+Macrophage_86_36hr" = "T86+M_36hr",
    "Macrophage_36hr" = "M(86)_36hr",   
    "Organoid_86_36hr" = "O86_36hr",
    "Organoid+Macrophage_86_36hr" = "O86+M_36hr(1)", 
    "Organoid+Macrophage_86_36hr(2)" = "O86+M_36hr(2)",
    
    "Tumoroid_92_24hr" = "T92_24hr",
    #"Tumoroid_96_24hr" = "T96_24hr",
    "Tumoroid+Macrophage_92_24hr" = "T92+M_24hr",
    #"Tumoroid+Macrophage_96_24hr" = "T96+M_24hr",
    "Macrophage_b1_24hr" = "M(92)_24hr",   
    #"Macrophage_b2_24hr" = "M(92)_24hr(2)",
    "Organoid_92_24hr" = "O92_24hr",
    "Organoid+Macrophage_92_24hr" = "O92+M_24hr"
)

stype_mapping = c("Macrophage" = "M",  "Organoid+Macrophage" = "O+M",  "Tumoroid+Macrophage" = "T+M")
stype_color <- c("M" = "#54278f", "O+M" = "#41ab5d", "T+M" = "#e7298a")
dset_color = c(
    "M(28)_24hr" = "#8c96c6",   
    "M(28)_48hr" = "#88419d",   
    "T28_24hr" = "#fb9a99",  
    "T28_48hr" = "#e31a1c",
    "T28+M_24hr" = "#66c2a4", 
    "T28+M_48hr" = "#238b45",
    
    "M(86)_36hr" = "#3f007d",   
    "O86_36hr" = "#d9f0a3",
    "O86+M_36hr(1)" = "#78c679", 
    "O86+M_36hr(2)" = "#41ab5d",
    "T86_36hr" = "#fec44f",
    "T86+M_36hr" = "#ec7014",
    
    "M(92)_24hr" = "#9e9ac8",   
    "T92_24hr" = "#fc9272",
    #"T96_24hr" = "#f768a1",
    "O92_24hr" = "#7fcdbb",
    "O92+M_24hr" = "#1d91c0",
    "T92+M_24hr" = "#ef3b2c"
    #"T96+M_24hr" = "#ae017e"
)

cur_proj$Dataset2 <- mapping_name[cur_proj$Dataset]
cur_proj$Dataset2 <- factor(cur_proj$Dataset2, levels = names(dset_color))
saveRDS(cur_proj, paste0(mstatePath, test_text, "cur_proj.rds"))

g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Dataset2", pal=dset_color, size = .4, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "umap_dataset.pdf"), g1, width = 4.5, height=2.5, units = "in", device = "pdf")

mac_proj <- cur_proj[cur_proj$Cluster %in% c(3,4,6) & grepl("Macrophage", cur_proj$Dataset),]
saveRDS(mac_proj, paste0(mstatePath, test_text, "mac_proj.rds"))
g1<-plotProj(mac_proj, dim_col = c(1,2), group.by="Dataset2", pal=dset_color, size = .8, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
glist <- list()
for(show_p in names(reactomeSets)) {
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    cut_qt = .975
    plot_df$activity_score[plot_df$activity_score > quantile(plot_df$activity_score, cut_qt)] <-  quantile(plot_df$activity_score, cut_qt)
    glist[[show_p]] <- plotProj(plot_df, dim_col = which(colnames(plot_df) %in% c("UMAP_1", "UMAP_2")), group.by = "activity_score", pal = "BlueGreenRed", size = .7) +
        theme(text=element_text(family = "Helvetica", size=9),
              legend.position="top",
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              # axis.ticks.y=element_blank(),
              # axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0))+
        ggtitle(show_p) +
        monocle:::monocle_theme_opts()
}

pdf(paste0(mplotPath, test_text, "activity_umap", ".pdf"),  width = 4, height=4)
invisible(lapply(glist, print))
dev.off()


glist <- list()
plot_meta$SampleType2 <- stype_mapping[plot_meta$SampleType]
plot_meta$Patient2 <- plot_meta$Patient
plot_meta$Patient2[plot_meta$Dataset2 %in% c("M(28)_24hr", "M(28)_48hr")] <- "28"
plot_meta$Patient2[plot_meta$Dataset2 %in% c("M(86)_36hr")] <- "86"
plot_meta$Patient2[plot_meta$Dataset2 %in% c("M(92)_24hr")] <- "92"
plot_meta$Patient2 <- paste0("P",plot_meta$Patient2)
for(show_p in names(reactomeSets)) {
    plot_df <- cbind(plot_meta, t(matrixAUC[show_p, rownames(plot_meta), drop=F]))
    #plot_df <- plot_df[plot_df$Dataset2 != "T96+M_24hr",]
    colnames(plot_df)[colnames(plot_df) == show_p] <- "activity_score"
    plot_df$activity_score_scaled <- plot_df$activity_score
    plot_df$activity_score_scaled[which(plot_df$Patient2 == "P28")] <- scale(plot_df$activity_score[which(plot_df$Patient2 == "P28")])
    plot_df$activity_score_scaled[which(plot_df$Patient2 == "P86")] <- scale(plot_df$activity_score[which(plot_df$Patient2 == "P86")])
    plot_df$activity_score_scaled[which(plot_df$Patient2 == "P92")] <- scale(plot_df$activity_score[which(plot_df$Patient2 == "P92")])
    
    glist[[show_p]] <- ggplot(plot_df, aes_string(y = "activity_score_scaled", fill = "SampleType2")) +
        geom_boxplot(outlier.colour = NA)+
        scale_fill_manual(values = stype_color)+
        theme(text=element_text(family = "Helvetica", size=8),
              axis.text = element_text(family = "Helvetica", size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              legend.margin=margin(0,0,0,0))+
        ggtitle(show_p)+
        ylim(c(-2.5,2.5)) + 
        xlab(NULL)+
        ylab("Activity score") +
        theme_bw() + 
        facet_wrap(~Patient2)
}
pdf(paste0(mplotPath, test_text, "activity_stype", ".pdf"),  width = 4, height=2.2)
invisible(lapply(glist, print))
dev.off()

write.xlsx(as.data.frame.matrix(table(plot_meta$SampleType2, plot_meta$Patient2)),  paste0(msheetPath, test_text, "mac_sig_cellcount.xlsx"), rowNames=T)


# Anova and post-hoc test
for(pp in c("P28", "P86", "P92")) {
    anova_list<-list()
    tukey_list <- list()
    plot_df <- cbind(plot_meta, t(matrixAUC[, rownames(plot_meta), drop=F]))
    plot_df <- plot_df[which(plot_df$Patient2 == pp),]
    for(show_p in names(reactomeSets)) {
        df <- plot_df[,c(show_p, "SampleType2")]
        colnames(df) <- c("activity_score", "sample_type")
        lm_res <- lm(activity_score ~ sample_type, data = df)
        aov_res <- aov(lm_res)
        anova_list[[show_p]] <- data.frame(unclass(summary(aov_res)), check.names = FALSE, stringsAsFactors = FALSE)
        tukey_res <- TukeyHSD(aov_res)
        tukey_res <- as.data.frame(tukey_res$sample_type)
        tukey_res$significance <- ifelse(tukey_res[["p adj"]] <= 0.0001, "****", ifelse(tukey_res[["p adj"]] <= 0.001, "***", ifelse(tukey_res[["p adj"]] <= 0.01, "**", ifelse(tukey_res[["p adj"]] <= 0.05, "*", "ns"))))
        tukey_list[[show_p]] <- tukey_res
    }
    write.xlsx(anova_list, paste0(msheetPath, test_text, "anova_results_", pp, ".xlsx"),  row.names = T)
    write.xlsx(tukey_list, paste0(msheetPath, test_text, "tukey_list_", pp, ".xlsx"), row.names = T)
}

plot_df <- cbind(plot_meta, t(matrixAUC[, rownames(plot_meta), drop=F]))
write.xlsx(plot_df[,c("Patient", "SampleType", "C1QC", "IL1B", "SPP1")],  paste0(msheetPath, test_text, "mac_sig.xlsx"), rowNames=T)



# Plot for individual patient
plot_vis <- clist$`Coculture experiment 1 [IFF]`
plot_proj <- cbind(plot_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[plot_vis@idx,])
plot_proj$Dataset2 <- mapping_name[plot_proj$Dataset]
plot_proj$Dataset2 <- factor(plot_proj$Dataset2, levels = names(dset_color))
plot_proj$Cluster <- plot_vis@pmeta$Cluster
g1<-plotProj(plot_proj, dim_col = c(1,2), group.by="Dataset2", pal=dset_color, size = .4, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "coculture1_umap_dataset.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")


plot_proj$Spp1_expr <- eset@assayData$norm_exprs[fData(eset)$gene_short_name == "SPP1",rownames(plot_proj)]
write.xlsx(plot_proj[,c("UMAP_1", "UMAP_2", "Dataset2", "Spp1_expr")], paste0(msheetPath, test_text, "coculture1_umap_dataset.xlsx"), row.names = T)


# Plot distribution shift of macrophage
mac_proj <- plot_proj[plot_proj$Cluster == 2 & plot_proj$SampleType != "Tumoroid",]
mac_proj$SampleType2 <- c("Macrophage" = "M", "Tumoroid+Macrophage" = "T+M")[mac_proj$SampleType]
sset_color = c(
    "M" = "#88419d",   
    "T+M" = "#238b45"
)
nbins=7
g1<-ggplot(mac_proj, aes_string("UMAP_1", "UMAP_2")) + 
    # stat_density2d(data = mac_proj[mac_proj$SampleType2 == "M",], fill = "#88419d", aes(color = "#88419d"), geom="polygon", alpha = .1, bins = 5) +
    # stat_density2d(data = mac_proj[mac_proj$SampleType2 == "T+M",], fill = "#238b45", aes(color = "#238b45"), geom="polygon", alpha = .1, bins = 5) + 
    stat_density2d(aes(fill = SampleType2, color = SampleType2), geom="polygon", alpha = .2, bins = nbins) +
    scale_color_manual(values = sset_color)+
    scale_fill_manual(values = sset_color)+
    xlim(c(-14,-7)) + ylim(c(-.5,4))+
    # scale_color_identity(name = "Sample type",
    #                      breaks = c("#88419d", "#238b45"),
    #                      labels = c("M", "T+M"),
    #                      guide = "legend") + 
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
ggsave(paste0(mplotPath, test_text, "coculture1_macrophage_distribution.pdf"), g1, width = 2.8, height=2, units = "in", device = "pdf")


# Expression of Spp1
mac_proj$Spp1 <- eset@assayData$norm_exprs[fData(eset)$gene_short_name == "SPP1",rownames(mac_proj)]
expr_cut <- quantile(mac_proj$Spp1, .975)
mac_proj$Spp1[mac_proj$Spp1 > expr_cut] = expr_cut
g1<-ggplot(mac_proj, aes_string("UMAP_1", "UMAP_2")) + 
    geom_point(aes(color = Spp1), size = .4) +
    scale_color_gradientn(colors = get_numeric_color("viridis")) +
    xlim(c(-14,-7)) + ylim(c(-0.5,4))+
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
ggsave(paste0(mplotPath, test_text, "coculture1_macrophage_Spp1.pdf"), g1, width = 2.8, height=2, units = "in", device = "pdf")


# Plot for patient 86
plot_vis <- clist$`Coculture experiment 20201216 [IFF, final]`
plot_proj <- cbind(plot_vis@proj$`UMAP-2D [15PC, IFG]`, pData(eset)[plot_vis@idx,])
plot_proj$Dataset2 <- mapping_name[plot_proj$Dataset]
plot_proj$Dataset2 <- factor(plot_proj$Dataset2, levels = names(dset_color))
plot_proj$Cluster <- plot_vis@pmeta$Cluster
g1<-plotProj(plot_proj, dim_col = c(1,2), group.by="Dataset2", pal=dset_color, size = .4, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "top", legend.title = NULL, keyheight=.5, rastr = T) + 
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
ggsave(paste0(mplotPath, test_text, "coculture86_umap_dataset.pdf"), g1, width = 3, height=3, units = "in", device = "pdf")


plot_proj$S100P_expr <- eset@assayData$norm_exprs[fData(eset)$gene_short_name == "S100P",rownames(plot_proj)]
write.xlsx(plot_proj[,c("UMAP_1", "UMAP_2", "Dataset2", "S100P_expr")], paste0(msheetPath, test_text, "coculture86_umap_dataset.xlsx"), row.names = T)



# Plot distribution shift of macrophage
mac_proj <- plot_proj[plot_proj$Cluster == 2 & !plot_proj$SampleType %in% c("Tumoroid", "Organoid"),]
mac_proj$SampleType2 <- c("Macrophage" = "M", "Organoid+Macrophage" = "O+M", "Tumoroid+Macrophage" = "T+M")[mac_proj$SampleType]
sset_color = c(
    "M" = "#3f007d",   
    "T+M" = "#ec7014",
    "O+M" = "#41ab5d"
)
nbins=7
g1<-ggplot(mac_proj, aes_string("UMAP_1", "UMAP_2")) + 
    # stat_density2d(data = mac_proj[mac_proj$SampleType2 == "M",], fill = sset_color["M"], aes(color = sset_color["M"]), geom="polygon", alpha = .1, bins = 6) +
    # stat_density2d(data = mac_proj[mac_proj$SampleType2 == "O+M",], fill = sset_color["O+M"], aes(color = sset_color["O+M"]), geom="polygon", alpha = .1, bins = 6) + 
    # stat_density2d(data = mac_proj[mac_proj$SampleType2 == "T+M",], fill = sset_color["T+M"], aes(color = sset_color["T+M"]), geom="polygon", alpha = .1, bins = 6) + 
    stat_density2d(aes(fill = SampleType2, color = SampleType2), geom="polygon", alpha = .2, bins = nbins) +
    scale_color_manual(values = sset_color)+
    scale_fill_manual(values = sset_color)+
    xlim(c(-7,3)) + ylim(c(-10.5,-4.5))+
    # scale_color_identity(name = "Sample type",
    #                      breaks = c(sset_color["M"], sset_color["O+M"], sset_color["T+M"]),
    #                      labels = c("M", "O+M", "T+M"),
    #                      guide = "legend") + 
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
ggsave(paste0(mplotPath, test_text, "coculture86_macrophage_distribution.pdf"), g1, width = 2.8, height=2, units = "in", device = "pdf")

organoid_proj <- plot_proj[plot_proj$Cluster %in% c(1) & plot_proj$SampleType %in% c("Organoid", "Organoid+Macrophage"),]
organoid_proj$SampleType2 <- c("Organoid" = "O", "Organoid+Macrophage" = "O+M")[organoid_proj$SampleType]
sset_color = c(
    "O" = "#3690c0", 
    "O+M" = "#41ab5d"
)
nbins=6
g1<-ggplot(organoid_proj, aes_string("UMAP_1", "UMAP_2")) + 
    stat_density2d(aes(fill = SampleType2, color = SampleType2), geom="polygon", alpha = .2, bins = nbins) +
    xlim(c(-9,-1)) + ylim(c(3,10))+
    scale_color_manual(values = sset_color)+
    scale_fill_manual(values = sset_color)+
    # scale_color_identity(name = "Sample type",
    #                      breaks = c(sset_color["T"], sset_color["T+M"]),
    #                      labels = c("T",  "T+M"),
    #                      guide = "legend") + 
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
ggsave(paste0(mplotPath, test_text, "coculture86_organoid_distribution.pdf"), g1, width = 2.8, height=2, units = "in", device = "pdf")

tumoroid_proj <- plot_proj[plot_proj$Cluster %in% c(3,4) & !plot_proj$SampleType %in% c("Macrophage", "Organoid", "Organoid+Macrophage"),]
tumoroid_proj$SampleType2 <- c("Tumoroid" = "T", "Tumoroid+Macrophage" = "T+M")[tumoroid_proj$SampleType]
sset_color = c(
    "T" = "#fec44f",   
    "T+M" = "#ec7014"
)
g1<-ggplot(tumoroid_proj, aes_string("UMAP_1", "UMAP_2")) + 
    stat_density2d(aes(fill = SampleType2, color = SampleType2), geom="polygon", alpha = .2, bins = nbins) +
    xlim(c(7,11)) + ylim(c(-5,7))+
    scale_color_manual(values = sset_color)+
    scale_fill_manual(values = sset_color)+
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
ggsave(paste0(mplotPath, test_text, "coculture86_tumoroid_distribution.pdf"), g1, width = 2.8, height=2, units = "in", device = "pdf")

tumoroid_proj$S100P <- eset@assayData$norm_exprs[fData(eset)$gene_short_name == "S100P",rownames(tumoroid_proj)]
expr_cut <- quantile(tumoroid_proj$S100P, .975)
tumoroid_proj$S100P[tumoroid_proj$S100P > expr_cut] = expr_cut
g1<-ggplot(tumoroid_proj, aes_string("UMAP_1", "UMAP_2")) + 
    geom_point(aes(color = S100P), size = .4) +
    scale_color_gradientn(colors = get_numeric_color("viridis")) +
    xlim(c(7,11)) + ylim(c(-5,7))+
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
ggsave(paste0(mplotPath, test_text, "coculture1_tumoroid_S100P.pdf"), g1, width = 2.8, height=2.2, units = "in", device = "pdf")




cur_proj <- readRDS(paste0(mstatePath, test_text, "cur_proj.rds"))
epi_proj <- cur_proj[!cur_proj$Cluster %in% c(3,4,6) & cur_proj$SampleType != "Macrophage",]
epi_proj$CD44_expr <- eset@assayData$norm_exprs[which(fData(eset)$gene_short_name == "CD44"),rownames(epi_proj)]
stype_mapping = c("Organoid" = "O",  "Organoid+Macrophage" = "O+M",  "Tumoroid" = "T", "Tumoroid+Macrophage" = "T+M")
epi_proj$SampleType2 <- stype_mapping[epi_proj$SampleType]
stype_color <- c("O" = "#3690c0", "O+M" = "#41ab5d", "T" = "#fd8d3c", "T+M" = "#e7298a")

plot_df <- epi_proj
plot_df$SampleType2 <- factor(plot_df$SampleType2, levels = names(stype_color))
g1 <- ggplot(plot_df, aes_string(x = "SampleType2", y = "CD44_expr", fill = "SampleType2")) +
    geom_boxplot(outlier.colour = NA)+
    scale_fill_manual(values = stype_color)+
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.margin=margin(0,0,0,0))+
    ggtitle("CD44")+
    xlab(NULL)+
    ylim(c(0,2.5)) + 
    theme_bw() + 
    facet_wrap(~Patient)

pdf(paste0(mplotPath, test_text, "CD44_expr", ".pdf"),  width = 4, height=2.7)
g1
dev.off()

write.xlsx(plot_df[,c("SampleType2", "Patient", "CD44_expr")],  paste0(msheetPath, test_text, "epi_cd44_expr.xlsx"), rowNames=T)
write.xlsx(as.data.frame.matrix(table(plot_df$SampleType2, plot_df$Patient)),  paste0(msheetPath, test_text, "epi_cd44_cellcount.xlsx"), rowNames=T)

# Anova and post-hoc test
for(pp in c("28", "86", "92")) {
    anova_list<-list()
    tukey_list <- list()
    plot_df <- epi_proj[which(epi_proj$Patient == pp),]
    for(show_p in "CD44_expr") {
        df <- plot_df[,c(show_p, "SampleType2")]
        colnames(df) <- c("expr", "sample_type")
        lm_res <- lm(expr ~ sample_type, data = df)
        aov_res <- aov(lm_res)
        anova_list[[show_p]] <- data.frame(unclass(summary(aov_res)), check.names = FALSE, stringsAsFactors = FALSE)
        tukey_res <- TukeyHSD(aov_res)
        tukey_res <- as.data.frame(tukey_res$sample_type)
        tukey_res$significance <- ifelse(tukey_res[["p adj"]] <= 0.0001, "****", ifelse(tukey_res[["p adj"]] <= 0.001, "***", ifelse(tukey_res[["p adj"]] <= 0.01, "**", ifelse(tukey_res[["p adj"]] <= 0.05, "*", "ns"))))
        tukey_list[[show_p]] <- tukey_res
    }
    write.xlsx(anova_list, paste0(msheetPath, test_text, "CD44_anova_results_", pp, ".xlsx"),  row.names = T)
    write.xlsx(tukey_list, paste0(msheetPath, test_text, "CD44_tukey_list_", pp, ".xlsx"), row.names = T)
}

## Cell cycle analysis for epithelial cells

epi_proj$Dataset2 <- as.character(epi_proj$Dataset2)
epi_proj$Dataset2[epi_proj$Dataset2 %in% c("O86+M_36hr(1)", "O86+M_36hr(2)")] = "O86+M_36hr"
epi_proj$Dataset2[epi_proj$Dataset2 %in% c("T28+M_24hr", "T28+M_48hr")] = "T28+M"
epi_proj$Dataset2[epi_proj$Dataset2 %in% c("T28_24hr", "T28_48hr")] = "T28"
epi_proj$Cell_cycle_phase <- eset$Cell_cycle_phase[match(rownames(epi_proj), colnames(eset))] ####
cc_frac<-as.data.frame.matrix(table(epi_proj$Dataset2, epi_proj$Cell_cycle_phase))
cc_frac <- cc_frac[!rownames(cc_frac) %in% c("M(28)_24hr", "M(28)_48hr", "M(86)_36hr", "M(92)_24hr"),]
cc_frac <- cc_frac/rowSums(cc_frac)
cc_frac$G2M_S <- cc_frac$G2M+ cc_frac$S
cc_frac$SampleType2 <- epi_proj$SampleType2[match(rownames(cc_frac), epi_proj$Dataset2)]
cc_frac$Patient <- epi_proj$Patient[match(rownames(cc_frac), epi_proj$Dataset2)]
g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G2M_S", fill = "SampleType2")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = stype_color)+
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
    scale_fill_manual(values = stype_color)+
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

cc_decrease <- as.data.frame(c(
    "O86+M" = cc_frac$G2M_S[2] - cc_frac$G2M_S[1],
    "O92+M" = cc_frac$G2M_S[4] - cc_frac$G2M_S[3],
    "T28+M" = cc_frac$G2M_S[6] - cc_frac$G2M_S[5],
    "T86+M" = cc_frac$G2M_S[8] - cc_frac$G2M_S[7],
    "T92+M" = cc_frac$G2M_S[10] - cc_frac$G2M_S[9]
))
colnames(cc_decrease) = "cc_diff"
cc_decrease$Patient <- c("86", "92", "28", "86", "92")
cc_decrease$comparison <- c("O+M vs O", "O+M vs O", "T+M vs T","T+M vs T","T+M vs T")
use_color = c(
    "T+M vs T" = "#e7298a",
    "O+M vs O" = "#41ab5d"
)
g1<-ggplot(cc_decrease, aes_string(x = "Patient", y = "cc_diff", fill = "comparison")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = use_color) + 
    theme(text=element_text(family = "Helvetica", size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL) + 
    theme_bw()

pdf(paste0(mplotPath, test_text, "cc_reduce", ".pdf"),  width = 3, height=1)
g1
dev.off()

write.xlsx(cc_frac, paste0(msheetPath,test_text, "epi_cc_frac.xlsx"), rowNames=T)


deg_tbls <- list(T28 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid_28_24hr_Tumoroid_28_48hr_vs_Tumoroid+Macrophage_28_24hr_Tumoroid+Macrophage_28_48hr_2021-02-22_de_significant.xlsx"))[[2]],
                  O86 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Organoid+Macrophage_86_36hr_vs_Organoid_86_36hr_2021-02-22_de_significant.xlsx"))[[2]],
                  T86 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid+Macrophage_86_36hr_vs_Tumoroid_86_36hr_2021-02-22_de_significant.xlsx"))[[2]],
                  O92 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Organoid+Macrophage_92_24hr_vs_Organoid_92_24hr_2021-02-22_de_significant.xlsx"))[[2]],
                  T92 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid+Macrophage_92_24hr_vs_Tumoroid_92_24hr_2021-02-22_de_significant.xlsx"))[[2]]
                  #T96 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid+Macrophage_96_24hr_vs_Tumoroid_96_24hr_2021-02-22_de_significant.xlsx"))[[2]]
)

deg_count <- data.frame(deg_count = sapply(deg_tbls, nrow))
deg_count$comparison <- factor(rownames(deg_count), levels = rownames(deg_count))
deg_count$Patient = ifelse(grepl("28", deg_count$comparison), "28", 
                           ifelse(grepl("86", deg_count$comparison), "86",
                                  "92"))
deg_count$SampleType = ifelse(grepl("T", deg_count$comparison), "T+M vs T", 
                           ifelse(grepl("O", deg_count$comparison), "O+M vs O",
                                  ""))
use_color = c(
    "T+M vs T" = "#e7298a",
    "O+M vs O" = "#41ab5d"
)
g1<-ggplot(deg_count, aes_string(x = "Patient", y = "deg_count", fill = "SampleType")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = use_color) + 
    theme(text=element_text(family = "Helvetica", size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL) + 
    theme_bw()

pdf(paste0(mplotPath, test_text, "deg_count", ".pdf"),  width = 2.8, height=1)
g1
dev.off()

write.xlsx(deg_count, paste0(msheetPath,test_text, "epi_deg_count.xlsx"), rowNames=T)




# Cell cycle analysis for macrophage
cur_vis <- clist$`Coculture experiment final set [IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Cluster
cur_proj$Dataset2 <- mapping_name[cur_proj$Dataset]
mac_proj <- cur_proj[cur_proj$Cluster %in% c(3,4,6),]
mac_proj$Dataset2 <- as.character(mac_proj$Dataset2)
mac_proj$SampleType2 <- stype_mapping[mac_proj$SampleType]
mac_proj$Dataset2[mac_proj$Dataset2 %in% c("O86+M_36hr(1)", "O86+M_36hr(2)")] = "O86+M_36hr"
mac_proj$Dataset2[mac_proj$Dataset2 %in% c("T28+M_24hr", "T28+M_48hr")] = "T28+M"
mac_proj$Dataset2[mac_proj$Dataset2 %in% c("M(28)_24hr", "M(28)_48hr")] = "M(28)"
mac_proj$Dataset2[mac_proj$Dataset2 %in% c("T28_24hr", "T28_48hr")] = "T28"

cc_frac<-as.data.frame.matrix(table(mac_proj$Dataset2, mac_proj$Cell_cycle_phase))
cc_frac <- cc_frac[!rownames(cc_frac) %in% c("O86_36hr", "T28", "T86_36hr"),]
cc_frac <- cc_frac/rowSums(cc_frac)
cc_frac$G2M_S <- cc_frac$G2M+ cc_frac$S
cc_frac$SampleType2 <- mac_proj$SampleType2[match(rownames(cc_frac), mac_proj$Dataset2)]
cc_frac$Patient <- mac_proj$Patient[match(rownames(cc_frac), mac_proj$Dataset2)]
cc_frac$Patient[rownames(cc_frac) == "M(28)"] = "28"
cc_frac$Patient[rownames(cc_frac) == "M(86)_36hr"] = "86"
cc_frac$Patient[rownames(cc_frac) == "M(92)_24hr"] = "92"
g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G2M_S", fill = "SampleType2")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = stype_color)+
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.margin=margin(0,0,0,0))+
    ggtitle("G2M/S fraction")+
    xlab(NULL)+
    theme_bw() 
pdf(paste0(mplotPath, test_text, "mac_G2M_S_frac", ".pdf"),  width = 4, height=2.7)
g1
dev.off()
g1 <- ggplot(cc_frac, aes_string(x = "Patient", y = "G1", fill = "SampleType2")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = stype_color)+
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.margin=margin(0,0,0,0))+
    ggtitle("G1 fraction")+
    xlab(NULL)+
    theme_bw() 
pdf(paste0(mplotPath, test_text, "mac_G1_frac", ".pdf"),  width = 4, height=2.7)
g1
dev.off()

cc_increase <- as.data.frame(c(
    "O86+M" = cc_frac$G2M_S[which(rownames(cc_frac)=="O86+M_36hr")] - cc_frac$G2M_S[which(rownames(cc_frac)=="M(86)_36hr")],
    "O92+M" = cc_frac$G2M_S[which(rownames(cc_frac)=="O92+M_24hr")] - cc_frac$G2M_S[which(rownames(cc_frac)=="M(92)_24hr")],
    "T28+M" = cc_frac$G2M_S[which(rownames(cc_frac)=="T28+M")] - cc_frac$G2M_S[which(rownames(cc_frac)=="M(28)")],
    "T86+M" = cc_frac$G2M_S[which(rownames(cc_frac)=="T86+M_36hr")] - cc_frac$G2M_S[which(rownames(cc_frac)=="M(86)_36hr")],
    "T92+M" = cc_frac$G2M_S[which(rownames(cc_frac)=="T92+M_24hr")] - cc_frac$G2M_S[which(rownames(cc_frac)=="M(92)_24hr")]
))
colnames(cc_increase) = "cc_diff"
cc_increase$Patient <- c("86", "92", "28", "86", "92")
cc_increase$comparison <- c("O+M vs O", "O+M vs O", "T+M vs T","T+M vs T","T+M vs T")
use_color = c(
    "T+M vs T" = "#e7298a",
    "O+M vs O" = "#41ab5d"
)
g1<-ggplot(cc_increase, aes_string(x = "Patient", y = "cc_diff", fill = "comparison")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = use_color) + 
    theme(text=element_text(family = "Helvetica", size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL) + 
    theme_bw()

pdf(paste0(mplotPath, test_text, "mac_cc_increase", ".pdf"),  width = 2.9, height=1)
g1
dev.off()


write.xlsx(cc_frac, paste0(msheetPath,test_text, "mac_cc_frac.xlsx"), rowNames=T)


deg_tbls <- list(
    O_M86 = read_excel_allsheets(paste0(msheetPath, "CC_mac_OM86_vs_M86_2022-08-15_de_significant.xlsx"))[[2]],
    O_M92 = read_excel_allsheets(paste0(msheetPath, "CC_mac_OM92_vs_M92_2022-08-15_de_significant.xlsx"))[[2]],
    T_M28 = read_excel_allsheets(paste0(msheetPath, "CC_mac_TM28_vs_M28_2022-08-15_de_significant.xlsx"))[[2]],
    T_M86 = read_excel_allsheets(paste0(msheetPath, "CC_mac_TM86_vs_M86_2022-08-15_de_significant.xlsx"))[[2]],
    T_M92 = read_excel_allsheets(paste0(msheetPath, "CC_mac_TM92_vs_M92_2022-08-15_de_significant.xlsx"))[[2]]
)

deg_count <- data.frame(deg_count = sapply(deg_tbls, nrow))
deg_count$comparison <- factor(rownames(deg_count), levels = rownames(deg_count))
deg_count$Patient = ifelse(grepl("28", deg_count$comparison), "28", 
                           ifelse(grepl("86", deg_count$comparison), "86",
                                  "92"))
deg_count$SampleType = ifelse(grepl("T", deg_count$comparison), "T+M vs M", 
                              ifelse(grepl("O", deg_count$comparison), "O+M vs M",
                                     ""))
use_color = c(
    "T+M vs M" = "#e7298a",
    "O+M vs M" = "#41ab5d"
)
g1<-ggplot(deg_count, aes_string(x = "Patient", y = "deg_count", fill = "SampleType")) +
    geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity") +
    scale_fill_manual(values = use_color) + 
    theme(text=element_text(family = "Helvetica", size=9),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL) + 
    theme_bw()

pdf(paste0(mplotPath, test_text, "mac_deg_count", ".pdf"),  width = 2.8, height=1)
g1
dev.off()


write.xlsx(deg_count, paste0(msheetPath,test_text, "mac_deg_count.xlsx"), rowNames=T)




library(ggraph)

# Receptor missing in tumoroid
test_text <- "missing_rl_recovered_in_coculture_"
rec_bind <- readRDS(paste0(mstatePath, "human_rl_", "rec_bind.rds"))
lig_bind <- readRDS(paste0(mstatePath, "human_rl_", "lig_bind.rds"))
rec_myeloid <- rec_bind[rec_bind$source == "Myeloid",]
lig_myeloid <- lig_bind[lig_bind$source == "Myeloid",]

deg_tbls <- list(T28 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid_28_24hr_Tumoroid_28_48hr_vs_Tumoroid+Macrophage_28_24hr_Tumoroid+Macrophage_28_48hr_2021-02-22_de_significant.xlsx"))[[2]],
                 T86 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid+Macrophage_86_36hr_vs_Tumoroid_86_36hr_2021-02-22_de_significant.xlsx"))[[2]],
                 T92 = read_excel_allsheets(paste0(msheetPath, "CC_epi_Tumoroid+Macrophage_92_24hr_vs_Tumoroid_92_24hr_2021-02-22_de_significant.xlsx"))[[2]],
                 M28 = read_excel_allsheets(paste0(msheetPath, "CC_mac_Tumoroid+Macrophage_28_24hr_Tumoroid+Macrophage_28_48hr_vs_Macrophage_24hr_Macrophage_48hr_2021-03-19_de_significant.xlsx"))[[2]],
                 M86 = read_excel_allsheets(paste0(msheetPath, "CC_mac_Tumoroid+Macrophage_86_36hr_vs_Macrophage_36hr_2021-03-19_de_significant.xlsx"))[[2]],
                 M92 = read_excel_allsheets(paste0(msheetPath, "CC_mac_Tumoroid+Macrophage_92_24hr_vs_Macrophage_b1_24hr_2021-03-19_de_significant.xlsx"))[[2]]
)
deg_list <- lapply(deg_tbls, function(x)x$gene_short_name)

# Any RL pairs up in co-culture
lig_myeloid$ligand <- as.character(lig_myeloid$ligand); lig_myeloid$receptor <- as.character(lig_myeloid$receptor);
rec_myeloid$receptor <-as.character(rec_myeloid$receptor); rec_myeloid$ligand <-as.character(rec_myeloid$ligand)

lig_myeloid$lig_up_in_coculture <- lig_myeloid$ligand %in% c(deg_list$M28, deg_list$M86, deg_list$M92)
lig_myeloid$rec_up_in_coculture <- lig_myeloid$receptor %in% c(deg_list$T28, deg_list$T86, deg_list$T92)
lig_myeloid$either_up_in_coculture <- lig_myeloid$lig_up_in_coculture | lig_myeloid$rec_up_in_coculture
rec_myeloid$rec_up_in_coculture <- rec_myeloid$receptor %in% c(deg_list$M28, deg_list$M86, deg_list$M92) 
rec_myeloid$lig_up_in_coculture <- rec_myeloid$ligand %in% c(deg_list$T28, deg_list$T86, deg_list$T92)
rec_myeloid$either_up_in_coculture <- rec_myeloid$rec_up_in_coculture | rec_myeloid$lig_up_in_coculture
    
lig_myeloid<- lig_myeloid[complete.cases(lig_myeloid),]
verts <- data.frame(name = unique(c(as.character(lig_myeloid$ligand), as.character(lig_myeloid$receptor))))
verts$type <- verts$name %in% lig_myeloid$ligand
verts$size <- ifelse(verts$name %in% c(lig_myeloid$ligand[lig_myeloid$lig_up_in_coculture], lig_myeloid$receptor[lig_myeloid$rec_up_in_coculture]), 2, 1)
lig_graph<-igraph::graph_from_data_frame(lig_myeloid, vertices = verts)
pdf(paste0(mplotPath, test_text, "rl_tumoroidR_recover.pdf"), width = 5.5, height = 5)
ggraph(lig_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(color = factor(either_up_in_coculture)), show.legend = F) + 
    scale_edge_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    #scale_edge_width_manual(values = c("TRUE" = "1", "FALSE" = ".2")) + 
    geom_node_point(aes(shape = type, color = type), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = type),  show.legend = F, label.size = 0.25, size = verts$size * 2) + 
    scale_color_manual(values = c("TRUE" = "#6a3d9a", "FALSE" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()


rec_myeloid<- rec_myeloid[complete.cases(rec_myeloid),]
verts <- data.frame(name = unique(c(as.character(rec_myeloid$ligand), as.character(rec_myeloid$receptor))))
verts$type <- verts$name %in% rec_myeloid$receptor
rec_graph<-graph_from_data_frame(rec_myeloid, vertices = verts)
pdf(paste0(mplotPath, test_text, "tumor_ligand_lost.pdf"), width = 5.5, height = 2.8)
ggraph(rec_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(edge_linetype = factor(lig_down_in_tumoroid), color = factor(rec_up_in_coculture)), show.legend = F) + 
    scale_edge_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    scale_edge_width(range = c(0.2,1)) + 
    geom_node_point(aes(shape = type, color = type), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = type),  show.legend = F, label.size = 0.25, size= 2) + 
    scale_color_manual(values = c("TRUE" = "#6a3d9a", "FALSE" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()

lig_myeloid2 <- lig_myeloid[,c(1:4,7)]
colnames(lig_myeloid2) <- c("mg", "tg", "lwd", "source", "up_in_coculture")
rec_myeloid2 <- rec_myeloid[,c(1:4,7)]
colnames(rec_myeloid2) <- c("mg", "tg", "lwd", "source", "up_in_coculture")
rl_myeloid <- rbind(lig_myeloid2, rec_myeloid2)

verts <- data.frame(name = unique(c(as.character(rl_myeloid$mg), as.character(rl_myeloid$tg))))
verts$type <- ifelse(verts$name %in% lig_myeloid$receptor | verts$name %in% rec_myeloid$receptor, "receptor", "ligand")
verts$source <- ifelse(verts$name %in% rl_myeloid$tg, "T", "M")
verts <- verts[order(verts$source, verts$type),]
rl_graph<-graph_from_data_frame(rl_myeloid, vertices = verts)

pdf(paste0(mplotPath, test_text, "rl_recover.pdf"), width = 20, height = 20)
ggraph(rl_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(color = factor(up_in_coculture)), show.legend = F) + 
    scale_edge_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    #scale_edge_width_manual(values = c("TRUE" = "1", "FALSE" = ".2")) + 
    geom_node_point(aes(shape = type, color = source), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = source),  show.legend = F, label.size = 0.25) + 
    scale_color_manual(values = c("M" = "#6a3d9a", "T" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()

pdf(paste0(mplotPath, test_text, "rl_recover2.pdf"), width = 2.5, height = 2.5)
ggraph(rl_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(color = factor(up_in_coculture)), show.legend = F, width = .2) + 
    scale_edge_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    #scale_edge_width_manual(values = c("TRUE" = "1", "FALSE" = ".2")) + 
    geom_node_point(aes(shape = type, color = source), size = 1, show.legend = F) + 
    #geom_node_label(aes(label = name, color = source),  show.legend = F, label.size = 0.25) + 
    scale_color_manual(values = c("M" = "#6a3d9a", "T" = "black")) +
    scale_shape_manual(values = c("receptor" = 16, "ligand"= 17))+
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()





# Receptor/ligand missing in tumoroid for experiment 2
mac_de_list1 <- readRDS(paste0(mstatePath,"coculture_1st_plot_CoCulture_vs_Macrophage_de_sig_genes.rds"))
mac_de_list2 <- readRDS(paste0(mstatePath,"coculture_2nd_plot_CocultureTumoroid+Macrophage_vs_Macrophage_de_sig_genes.rds"))

test_text <- "missing_rl_recovered_in_coculture2_"
rec_bind <- readRDS(paste0(mstatePath, "human_rl_", "rec_bind.rds"))
lig_bind <- readRDS(paste0(mstatePath, "human_rl_", "lig_bind.rds"))
rec_myeloid <- rec_bind[rec_bind$source == "Myeloid",]
lig_myeloid <- lig_bind[lig_bind$source == "Myeloid",]

tumoroid86_de_tbl <- read_excel_allsheets(paste0(msheetPath, "G1_Tumor_86_vs_Tumoroid_86_36hr_2021-01-25_de_significant.xlsx"))
tumoroid86_coculture_de_tbl <- readRDS("~/Dropbox/ColonManuscript/states/coculture_2nd_plot_CocultureTumoroid+Macrophage_vs_Tumoroid_de_sig_genes.rds")

lig_myeloid$rec_down_in_tumoroid <- lig_myeloid$receptor %in% tumoroid86_de_tbl$Tumor_86$gene_short_name
lig_myeloid$rec_up_in_coculture <- lig_myeloid$receptor %in% tumoroid86_coculture_de_tbl$`CocultureTumoroid+Macrophage`$gene_name
lig_myeloid$lig_up_in_coculture <- lig_myeloid$ligand %in% mac_de_list2$`CocultureTumoroid+Macrophage`$gene_name
lig_myeloid$either_up_in_coculture <- lig_myeloid$rec_up_in_coculture | lig_myeloid$lig_up_in_coculture

rec_myeloid$lig_down_in_tumoroid <- rec_myeloid$ligand %in% tumoroid86_de_tbl$Tumor_86$gene_short_name
rec_myeloid$lig_up_in_coculture <- rec_myeloid$ligand %in% tumoroid86_coculture_de_tbl$`CocultureTumoroid+Macrophage`$gene_name
rec_myeloid$rec_up_in_coculture <- rec_myeloid$receptor %in% mac_de_list2$`CocultureTumoroid+Macrophage`$gene_name
rec_myeloid$either_up_in_coculture <- rec_myeloid$rec_up_in_coculture | rec_myeloid$lig_up_in_coculture

lig_myeloid<- lig_myeloid[complete.cases(lig_myeloid),]
verts <- data.frame(name = unique(c(as.character(lig_myeloid$ligand), as.character(lig_myeloid$receptor))))
verts$type <- verts$name %in% lig_myeloid$ligand
lig_graph<-igraph::graph_from_data_frame(lig_myeloid, vertices = verts)
pdf(paste0(mplotPath, test_text, "tumor_receptor_lost.pdf"), width = 5.5, height = 2.8)
ggraph(lig_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(edge_linetype = factor(rec_down_in_tumoroid), color = factor(either_up_in_coculture)), show.legend = F) + 
    scale_edge_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    #scale_edge_width_manual(values = c("TRUE" = "1", "FALSE" = ".2")) + 
    geom_node_point(aes(shape = type, color = type), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = type),  show.legend = F, label.size = 0.25, size= 2) + 
    scale_color_manual(values = c("TRUE" = "#6a3d9a", "FALSE" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()


rec_myeloid<- rec_myeloid[complete.cases(rec_myeloid),]
verts <- data.frame(name = unique(c(as.character(rec_myeloid$ligand), as.character(rec_myeloid$receptor))))
verts$type <- verts$name %in% rec_myeloid$receptor
rec_graph<-graph_from_data_frame(rec_myeloid, vertices = verts)
pdf(paste0(mplotPath, test_text, "tumor_ligand_lost.pdf"), width = 5.5, height = 2.8)
ggraph(rec_graph,  layout = 'linear', circular = T) + 
    geom_edge_arc(aes(edge_linetype = factor(lig_down_in_tumoroid), color = factor(either_up_in_coculture)), show.legend = F) + 
    scale_edge_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) + 
    #scale_edge_alpha(range = c(0.1,1)) +
    scale_edge_width(range = c(0.2,1)) + 
    geom_node_point(aes(shape = type, color = type), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = type),  show.legend = F, label.size = 0.25, size= 2) + 
    scale_color_manual(values = c("TRUE" = "#6a3d9a", "FALSE" = "black")) +
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()





# Intersection with in vivo gene list
mac_de_list1 <- readRDS(paste0(mstatePath,"coculture_1st_plot_CoCulture_vs_Macrophage_de_sig_genes.rds"))
mac_de_list2 <- readRDS(paste0(mstatePath,"coculture_2nd_plot_CocultureTumoroid+Macrophage_vs_Macrophage_de_sig_genes.rds"))
mac_invivo_list <- read_excel_allsheets(paste0(msheetPath,"macrophages in vivo_C1QC_vs_IL1B_vs_S100A8_vs_SPP1_2021-02-21_de_significant.xlsx"))

library(VennDiagram)
map_list <- list("mac_c1" = mac_de_list1$CoCulture$gene_name, "map_c2"= mac_de_list2$`CocultureTumoroid+Macrophage`$gene_name, "mac_invivo" = mac_invivo_list$SPP1$gene_short_name)
pdf(paste0(mplotPath, "spp1_macg_intersect.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(map_list, NULL, cex = 2, cat.fontface=4)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()

w3<-Reduce(intersect, map_list)
w1<-Reduce(intersect, map_list[c(1,3)])
w2 <- Reduce(intersect, map_list[c(2,3)])
wg <- Reduce(union, c(w1, w2, w3))
wtbl <- mac_invivo_list$SPP1
wtbl$in_w <- ifelse(mac_invivo_list$SPP1$gene_short_name %in% c(w1,w2), 2, ifelse(mac_invivo_list$SPP1$gene_short_name %in% w3, 3, 1))
write.xlsx(wtbl, paste0(msheetPath, "spp1_macg.xlsx"))












