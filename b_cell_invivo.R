

library(VisCello)
library(ggrastr)

is_windows = TRUE
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


test_text <- "B_plasma_hs_invivo_221227_"
cur_vis <- clist$`B cells in vivo [20221224, IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset)[cur_vis@idx,])
#cur_proj <- cur_proj[!cur_proj$Cell_type_B %in% c("Unannotated", "Macrophage-T cell doublet"), ]
ctype_color <- get_factor_color(unique(cur_proj$Cell_type_B), "Set2")
names(ctype_color) <- c(
    "B cell","Cycling B/Plasma cell",  "IGHM+ Plasma cell" ,   
    "IGHA+ Plasma cell (Normal)",
    "IGHA+ Plasma cell (Tumor)",  "IGHG+ Plasma cell"          
)
cur_proj$Cell_type_B <- factor(cur_proj$Cell_type_B, levels = names(ctype_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_type_B", pal=ctype_color, size = .5, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr=T) + 
    guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    theme(text=element_text(size=8),
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
ggsave(paste0(mplotPath, test_text, "B_ctype.pdf"), g1, width = 4, height=2.7, units = "in", device = "pdf")


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
    theme(text=element_text(size=8),
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
ggsave(paste0(mplotPath, test_text, "sampletype.pdf"), g1, width = 3.5, height=2.7, units = "in", device = "pdf")


nbin <- 5
g1<-ggplot(cur_proj, aes_string("UMAP_1", "UMAP_2")) + 
    stat_density2d(aes(fill = SampleType, color = SampleType), geom="polygon", alpha = .1, bins = nbin) +
    geom_point(aes(color=SampleType),size = .5, stroke=0)+
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
ggsave(paste0(mplotPath, test_text, "sampletype_distribution.pdf"), g1, width = 3.2, height=2.7, units = "in", device = "pdf")


cell_count<-as.data.frame.matrix(table(cur_proj$SampleType, cur_proj$Cell_type_B))
cell_frac <- as.data.frame(cell_count / rowSums(cell_count))
cell_frac$Sample_type <- rownames(cell_count)
plot_df <- reshape2::melt(cell_frac)
colnames(plot_df) <- c("Sample_type", "Cell_type", "Count")
plot_df$Sample_type <- factor(plot_df$Sample_type, levels = names(stype_color))
g1 = ggplot(plot_df, aes(x="", y = Count, fill = Cell_type)) + 
    geom_bar(stat = "identity", width=1) +
    coord_polar("y", start=0) + 
    facet_wrap(~Sample_type) +
    scale_fill_manual(values = ctype_color)+
    theme_bw() +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          strip.text = element_text(size=8),
          axis.text = element_text(size=8),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) 
ggsave(paste0(mplotPath, test_text, "sampletype_pie.pdf"), g1, width = 4, height=3, units = "in", device = "pdf")




plot_g <- c("MS4A1", "HLA-DQA1", "MKI67", "IGHM", "IGHA1", "IGHG1")
plot_proj = cur_proj[cur_proj$Cell_type_B!="T-Carcinoma cell doublet",]
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(plot_g,fData(eset)$gene_short_name), match(rownames(plot_proj), colnames(eset))])))
colnames(g_exprs) <- plot_g
glist <- list()
for(g in plot_g) {
    print(g)
    plot_proj$gene_expr <- g_exprs[[g]]
    ecut <- quantile(plot_proj$gene_expr[plot_proj$gene_expr > 0], .975)
    if(ecut < 1) ecut = 1
    plot_proj$gene_expr[plot_proj$gene_expr > ecut] = ecut
    glist[[g]]<-plotProj(plot_proj, dim_col = c(1,2), group.by="gene_expr", pal="viridis", size = .5, legend.title = g, layer_zero = T, rastr=T) + guides(color = guide_colorbar(barwidth = 10, barheight = 1, title = g, raster = T)) +
        guides(color = guide_colorbar(
            barwidth = 4, barheight = .5, title = g,
            title.theme = element_text(size = 8),
            label.theme = element_text(size = 8))) + 
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        theme(text=element_text(size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0))+
        monocle:::monocle_theme_opts()
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "gexpr_umap.pdf"), g1, width = 4, height=3.3, units = "in", device = "pdf")






plot_g <- c("MS4A1", "HLA-DQA1", "MKI67", "IGHM", "IGHA1", "IGHG1")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(plot_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- plot_g
save_proj <- cbind(cur_proj[,c("UMAP_1", "UMAP_2", "Patient", "SampleType","Cell_type_B")], g_exprs)
write.xlsx(save_proj, paste0(msheetPath, test_text, "save_proj.xlsx"), rowNames=T)


