



library(VisCello)
library(igraph)
savePath <- "../hcc_final_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)
scriptPath <- paste0("../hcc_cello/", "scripts/")

is_windows = FALSE
if(is_windows) {
    home_path = "C:/Users/qinzh/"
} else {
    home_path = "~"
}
mplotPath <- paste0(home_path, "/Dropbox/ColonManuscript/subplots/")
msheetPath <- paste0(home_path, "/Dropbox/ColonManuscript/sheets/")
mstatePath<- paste0(home_path, "/Dropbox/ColonManuscript/states/")


source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

eset <- readRDS("../hcc_final_cello/eset.rds")
clist <- readRDS("../hcc_final_cello/clist.rds")

library(AUCell)

test_text <- "signature_space_"

# Stem cell, differentiated signature derived based on DE in normal 
normal_degs <- readRDS(paste0(mstatePath,"normal_colon_colon_ctype_deg_sig_genes.rds"))
max_g_use <- 100
normal_degs <- lapply(normal_degs, function(x) as.character(x$gene_name[1:max_g_use]))


dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

expressed_g <- which(rowMeans(exprs(eset_epi) > 0) > .01)
expr_hc <- as.matrix(exprs(eset_epi[expressed_g,]))
rownames(expr_hc) <- make.names(fData(eset_epi)$gene_short_name[match(rownames(expr_hc), fData(eset_epi)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_hc)

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(normal_degs, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, "epi_normal_Ctype_0.2cells_AUC.rds"))

cells_AUC <- readRDS(paste0(mstatePath, "epi_normal_Ctype_0.2cells_AUC.rds"))
plot_df <- as.data.frame(t(getAUC(cells_AUC)))
plot_df <- cbind(plot_df, pData(eset_epi))
stype_color <- c("Colon" = "#1f78b4", "Organoid" = "#33a02c", "Tumor" = "#e31a1c", "Tumoroid" = "#ff7f00")
plot_df <- plot_df[plot_df$SampleType != "LiverMet", ]
show_group <- c("SC", "EC")
label_data <- plot_df %>% group_by_at("SampleType") %>% summarize_at(show_group, median)
g1<-ggplot(plot_df) +
    geom_point(aes_string(show_group[1], show_group[2], color="SampleType"),size=.5, stroke=0) +
    geom_text(data = label_data, aes_string(show_group[1], show_group[2], label="SampleType"), size = 3) + 
    scale_color_manual(values = stype_color) + 
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

# 3D scatter plot
library("scatterplot3d") 
colors <- stype_color[plot_df$SampleType]
plot_df <- as.data.frame(scale(t(getAUC(cells_AUC))))
plot_df <- cbind(plot_df, pData(eset_epi))
#scatterplot3d(plot_df[,c("SC", "TA", "EC")], pch = 16, color=colors, grid=TRUE, box=FALSE)
library(plotly)
fig <- plot_ly(plot_df, x = ~SC, y = ~TA, z = ~EC, color = ~SampleType, colors = stype_color, size=.1)
fig <- fig %>% add_markers()
fig
library(htmlwidgets)
saveWidget(fig, paste0(mplotPath,test_text,"sig_3d_scaled.html"), selfcontained = T)

# Pheatmap
plot_df <- as.data.table(t(getAUC(cells_AUC)))
plot_df$Dataset <- eset_epi$Dataset
plot_df <- plot_df[!grepl("LiverMet", plot_df$Dataset), ]
mean_score <- as.data.frame(plot_df[, lapply(.SD, mean), by = Dataset])
rownames(mean_score) <- mean_score$Dataset
mean_score$Dataset <- NULL
mean_score$PRO <- NULL

mean_score_scaled <- t(scale(t(mean_score)))

pdf(paste0(mplotPath, test_text, "epi_normal_Ctype_AUC_hmap_scaled.pdf"), width = 4, height = 2.5)
pheatmap(t(mean_score_scaled), clustering_method = "ward.D2", cluster_cols = T, cluster_rows = F, fontsize = 8, color = get_numeric_color("RdBu"))
dev.off()

mean_score_cut <- quantile(unlist(mean_score), .975)
mean_score[mean_score > mean_score_cut] = mean_score_cut
pdf(paste0(mplotPath, test_text, "epi_normal_Ctype_AUC_hmap.pdf"), width = 4, height = 2.4)
pheatmap(t(mean_score), clustering_method = "ward.D2", cluster_cols = T, cluster_rows = F,fontsize = 8, color = get_numeric_color("RdBu"))
dev.off()


# Plot as ternery
plot_df <- as.data.frame(t(getAUC(cells_AUC)))
plot_df <- cbind(plot_df, pData(eset_epi))
plot_df <- plot_df[plot_df$SampleType != "LiverMet", ]
stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Organoid" = "#4daf4a", 
    "Tumoroid" = "#ff7f00"
)
library(ggtern)
library(scales)
plot_df[1:7] <- rescale(as.matrix(plot_df[1:7]), to=c(0,100))
g1 <- ggtern(data=plot_df, aes(x=SC,y=TA, z=EC)) +
    stat_density_tern(aes(fill = SampleType, color=SampleType), geom="polygon", alpha = .15, bins = 6, bdl=0.01) +
    scale_fill_manual(values = stype_color) +
    scale_color_manual(values = stype_color) +
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=8),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))+
    theme_rgbw()
pdf(paste0(mplotPath, test_text, "epi_normal_Ctype_AUC_tern.pdf"), width = 5, height = 5)
g1
dev.off()




# Plot as 4 separate plots
glist <- list()
for(stype in names(stype_color)) {
    cur_df<- plot_df[plot_df$SampleType == stype,]
    glist[[stype]] <- ggtern(data=cur_df, aes(x=SC,y=TA, z=EC)) +
        stat_density_tern(aes(fill = SampleType, color=SampleType), geom="polygon", alpha = .15, bins = 6, bdl=0.01) +
        scale_fill_manual(values = stype_color) +
        scale_color_manual(values = stype_color) +
        theme(text=element_text(family = "Helvetica", size=8),
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0),
              #legend.box.margin=margin(-10,-10,-10,-10),
              plot.margin = unit(c(0,0.4,0,0.1), "cm"))+
        theme_rgbw() + 
        guides(fill = F, color = F)
}
g1<-do.call(grid.arrange,c(glist, ncol = 2))
ggsave(paste0(mplotPath, test_text, "epi_normal_Ctype_AUC_tern_separate.pdf"), g1, width = 5, height=5, units = "in", device = "pdf")

write.xlsx(plot_df[,1:7], paste0(msheetPath, test_text, "plot_df.xlsx"), rowNames=T)


# Define lineage plasticity as the entropy of cell states
# To compute entropy, first linearly rescale score to 0 - 1 to be a proxy of probability
use_states <- c("SC", "TA", "EC")
auc_score <- as.data.frame(t(getAUC(cells_AUC)))
cell_state_df  <- cbind(auc_score, pData(eset_epi))

shannon_entropy <- function(vec) {
    vec = vec/sum(vec)
    vec = vec[vec!=0]
    -sum(vec * log2(vec))
}
cell_state_df$entropy <- apply(cell_state_df[,use_states], 1, shannon_entropy)
library(ggbreak) 
cell_state_df <- cell_state_df[cell_state_df$SampleType!= "LiverMet",]
g1 <- ggplot(cell_state_df, aes(x = Dataset, y = entropy)) +
    geom_boxplot(aes(fill = SampleType), outlier.shape=NA, ) + 
    ylim(1,1.6) + 
    theme_classic() +
    scale_fill_manual(values = stype_color) +
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "epi_shannon_entropy.pdf"), g1, width = 5, height=5, units = "in", device = "pdf")




# Plot Signature score on umap

cells_AUC <- readRDS(paste0(mstatePath, "epi_normal_Ctype_0.2cells_AUC.rds"))
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
auc_df = as.data.frame(t(getAUC(cells_AUC)))
qt_cut = .975
auc_df = apply(auc_df, 2, function(x) {
    x[x > quantile(x,qt_cut)] = quantile(x,qt_cut)
    return(x)
})
cur_proj <- cbind(dep_vis@proj$`UMAP-2D [10PC, IFG]`, auc_df, pData(eset_epi))
library(ggrastr)

glist <- list()
show_stype = colnames(auc_df)[colnames(auc_df)!="PRO"]
for(stype in show_stype) {
    glist[[stype]] = plotProj(cur_proj, dim_col = c(1,2), 
                              group.by = stype, pal = get_numeric_color("BlueGreenRed"), 
                              size =.4, legend_barwidth = 5, legend_barheight=.8,
                              rastr = T)+
        theme(text=element_text(size=8),
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0),
              plot.margin = unit(c(0,0.4,0,0.1), "cm"))
}

# pdf(paste0(mplotPath, test_text, "colonsignature_umap", ".pdf"),  width = 3, height=3)
# invisible(lapply(glist, print))
# dev.off()
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "colonsignature_umap_cbn.pdf"), g1, width = 7, height=5, units = "in", device = "pdf")


write.xlsx(cur_proj[,c(1:12)], paste0(msheetPath, test_text, "colonsignature_umap_cbn.xlsx"), rowNames=T)


plot_proj <- cbind(auc_df, pData(eset_epi)[,c("SampleType", "Patient"), drop=F])
plot_proj <- plot_proj[plot_proj$Patient %in% c("08", "24", "28"),] # Only use those have paired info
plot_proj <- reshape2::melt(plot_proj[,colnames(plot_proj) != "Patient"], id.var = "SampleType")
colnames(plot_proj) <- c("SampleType", "cell_state", "score")
plot_proj <- plot_proj[plot_proj$SampleType!= "LiverMet" & plot_proj$cell_state != "PRO",]
stype_color <- c("Colon" = "#1f78b4", "Organoid" = "#33a02c", "Tumor" = "#e31a1c", "Tumoroid" = "#ff7f00")
g1 <- ggplot(plot_proj, aes_string(x = "cell_state", y = "score", fill = "SampleType")) +
    geom_boxplot(outlier.colour = NA, position = position_dodge(preserve = "single"))+
    #geom_jitter_rast(aes(color = SampleType),size = .1, stroke = 0, position=position_jitterdodge())+
    scale_fill_manual(values = stype_color) +
    ylim(c(0,.2))+
    theme_bw()+
    theme(text=element_text(size=8),
          axis.text = element_text(size=8),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.6, "cm"),
          #axis.text.x = element_text(angle=45, hjust=1, size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_signaturescore_boxplot_v1.pdf")),  width = 7, height=2)
g1
dev.off()


compare_groups <- list(
    c("Colon", "Organoid"), 
    c("Tumor", "Tumoroid")
)

test_group = "score"

pres_state <- sapply(show_stype, function(cur_ct) {
    print(cur_ct)
    cur_tbl = plot_proj[plot_proj$cell_state == cur_ct,]
    pres_ct <- list()
    for(cur_group in compare_groups) {
        print(cur_group)
        pres_ct[[paste0(make.names(cur_group), collapse = "_vs_")]] <- t.test(cur_tbl[[test_group]][cur_tbl$SampleType == cur_group[1]], cur_tbl[[test_group]][cur_tbl$SampleType == cur_group[2]])
    }
    sapply(pres_ct, function(x) x$p.value)
})

pres_state_star <- pres_state
pres_state_star[pres_state > 0.05] = "ns"
pres_state_star[pres_state <= 0.05] = "*"
pres_state_star[pres_state <= 0.01] = "**"
pres_state_star[pres_state <= 0.001] = "***"
pres_state_star[pres_state <= 0.0001] = "****"


write.xlsx(list(test_res = pres_state), paste0(msheetPath, test_text, "_score_stype_ttest_res.xlsx"), rowNames = T)
write.xlsx(list(test_res = pres_state_star), paste0(msheetPath, test_text, "_score_stype_ttest_res_star.xlsx"), rowNames = T)























############ NOT used #####################


# Check if signature of tumor from individual patients is preserved when cultured, SKIP
test_text <- "tumor_signature_preserve300_"
tumor_deg = read_excel_allsheets(paste0(msheetPath, "Epithelial cells [Final, wt OGN, IFF]_tumor08vstumor24vstumor28_2022-12-27_de_significant.xlsx"))

max_g_use <- 300
tumor_deg <- lapply(tumor_deg, function(x) as.character(x$gene_short_name[1:max_g_use]))

dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]


expressed_g <- which(rowMeans(exprs(eset_epi) > 0) > .01)
expr_hc <- as.matrix(exprs(eset_epi[expressed_g,]))
rownames(expr_hc) <- make.names(fData(eset_epi)$gene_short_name[match(rownames(expr_hc), fData(eset_epi)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_hc)

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(tumor_deg, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "epi_0.2cells_AUC.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "epi_0.2cells_AUC.rds"))
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
auc_df = as.data.frame(t(getAUC(cells_AUC)))
qt_cut = .99
auc_df = apply(auc_df, 2, function(x) {
    x[x > quantile(x,qt_cut)] = quantile(x,qt_cut)
    return(x)
})
cur_proj <- cbind(dep_vis@proj$`UMAP-2D [10PC, IFG]`, auc_df, pData(eset_epi))
library(ggrastr)

glist <- list()
show_stype = colnames(auc_df)
for(stype in show_stype) {
    glist[[stype]] = plotProj(cur_proj, dim_col = c(1,2), 
                              group.by = stype, pal = get_numeric_color("BlueGreenRed"), 
                              size =.4, legend_barwidth = 5, legend_barheight=.8,
                              rastr = T)+
        theme(text=element_text(size=8),
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0),
              plot.margin = unit(c(0,0.4,0,0.1), "cm"))
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "tumorsignature_umap_cbn.pdf"), g1, width = 7, height=2.5, units = "in", device = "pdf")

# Boxplot

glist <- list()
show_stype = colnames(auc_df)
plot_proj = cur_proj[cur_proj$Patient %in% c("08","24", "28") & cur_proj$SampleType %in% c("Tumor", "Tumoroid"),]
for(stype in show_stype) {
    glist[[stype]] = 
        ggplot(plot_proj, aes_string(x = "Patient", y = stype, fill = "SampleType")) +
        geom_boxplot(outlier.shape=NA, position = position_dodge(preserve = "single"), lwd=.3) + 
        scale_fill_manual(values = stype_color) + 
        theme_classic() + 
        theme(text=element_text(size=9),
              legend.text=element_text(size=9),
              axis.text = element_text(size=9),
              legend.margin=margin(0,0,0,0),
              #legend.box.margin=margin(-10,-10,-10,-10),
              plot.margin = unit(c(0,0.4,0,0.1), "cm"))
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "tumorsignature_boxplot_cbn.pdf"), g1, width = 10, height=2.3, units = "in", device = "pdf")


test_pair <- c("Tumor_08", "Tumoroid_08");test_col = "Tumor_08"
# test_pair <- c("Tumor_24", "Tumoroid_24");test_col = "Tumor_24"
# test_pair <- c("Tumor_28", "Tumoroid_28");test_col = "Tumor_28"
group1_idx <- which(plot_proj$Dataset == test_pair[1])
group2_idx <-  which(plot_proj$Dataset == test_pair[2])
g1_val <- plot_proj[[test_col]][group1_idx]
g2_val <- plot_proj[[test_col]][group2_idx]
t.test(g1_val, g2_val, alternative = "greater")

