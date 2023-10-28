




library(VisCello)
library(igraph)

scriptPath <- paste0("../hcc_cello/", "scripts/")

mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

eset_mouse <- readRDS("../mouse_colon_cello/eset.rds")
clist_mouse <- readRDS("../mouse_colon_cello/clist.rds")

library(AUCell)

test_text <- "signature_space_mouse_"

dep_vis <- clist_mouse$`Mouse normal colon epithelial [IFF]`
expr_ne <-eset_mouse[, dep_vis@idx]
expr_ne$Cluster <- dep_vis@pmeta$Cluster
expr_ne$Colon_cell_type <- NA
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(25,23,5,15,20,2,14,18,24)] = "TA"
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(3,22,10)] = "SC"
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(12,7,9)] = "PRO"
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(8,21,19,6,26)] = "CC"
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(28,27)] = "EEC"
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(13,16,17,1,4,11)] = "G"
expr_ne$Colon_cell_type[expr_ne$Cluster %in% c(29)] = "unannotated"
cc_level <- c("SC", "TA", "PRO", "CC", "G", "EEC", "unannotated")
expr_ne$Colon_cell_type <- factor(expr_ne$Colon_cell_type, levels = cc_level)
eset_mouse$Colon_cell_type <- expr_ne$Colon_cell_type[match(colnames(eset_mouse), colnames(expr_ne))]
eset_mouse$SampleType = as.character(eset_mouse$Dataset)
eset_mouse$SampleType[eset_mouse$SampleType %in% c("Colon1", "Colon2")] <- "Colon"
# saveRDS(eset_mouse, "../mouse_colon_cello/eset.rds")

cur_proj <- cbind(dep_vis@proj$`UMAP-2D [10PC, IFG]`, pData(expr_ne))
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



# Stem cell, differentiated signature derived based on DE in normal 
normal_degs <- readRDS(paste0(mstatePath,"normal_mouse_colon_colon_ctype_deg_sig_genes.rds"))
max_g_use <- 100
normal_degs <- lapply(normal_degs, function(x) as.character(x$gene_name[1:max_g_use]))

dep_vis <- clist_mouse$`Mouse epithelial [IFF]`
eset_epi <-eset_mouse[, dep_vis@idx]
expressed_g <- which(rowMeans(exprs(eset_epi) > 0) > .01)
expr_epi <- as.matrix(exprs(eset_epi[expressed_g,]))
rownames(expr_epi) <- make.names(fData(eset_epi)$gene_short_name[match(rownames(expr_epi), fData(eset_epi)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_epi)

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(normal_degs, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "epi_normal_Ctype_0.2cells_AUC.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "epi_normal_Ctype_0.2cells_AUC.rds"))
plot_df <- as.data.frame(t(getAUC(cells_AUC)))
plot_df <- cbind(plot_df, pData(eset_epi))

stype_color <- c("Colon" = "#1f78b4", "Tumor ApcMin" = "#bf5b17", "Tumor AOM" = "#e31a1c", "Tumor APKS" = "#df65b0")
show_group <- c("SC", "CC")
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
#scatterplot3d(plot_df[,c("SC", "TA", "CC")], pch = 16, color=colors, grid=TRUE, box=FALSE)
library(plotly)
fig <- plot_ly(plot_df, x = ~SC, y = ~TA, z = ~CC, color = ~SampleType, colors = stype_color, size=.1)
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
library(ggtern)
library(scales)
plot_df[1:6] <- rescale(as.matrix(plot_df[1:6]), to=c(0,100))
g1 <- ggtern(data=plot_df, aes(x=SC,y=TA, z=CC)) +
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
    glist[[stype]] <- ggtern(data=cur_df, aes(x=SC,y=TA, z=CC)) +
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


write.xlsx(plot_df[,c("SC", "TA", "PRO", "CC", "G", "EEC", "SampleType")], paste0(msheetPath, test_text, "epi_normal_Ctype_AUC_tern_separate.xlsx"), rowNames = T)



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


# Compute the same for macrophage stimulated case





