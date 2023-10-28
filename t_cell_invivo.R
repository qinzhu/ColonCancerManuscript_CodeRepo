

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


test_text <- "tcell_hs_invivo_221227_"
cur_vis <- clist$`T cells in vivo [20221224, IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, batch regressed]`, pData(eset)[cur_vis@idx,])
#cur_proj <- cur_proj[!cur_proj$Cell_type_T %in% c("Unannotated", "Macrophage-T cell doublet"), ]
#cur_proj <- cur_proj[!cur_proj$Cell_type_T %in% c("NK cell"), ]
tctype_color <- get_factor_color(unique(cur_proj$Cell_type_T), "Set2")
names(tctype_color) <- unique(cur_proj$Cell_type_T)
#cur_proj$Cell_type_T <- factor(cur_proj$Cell_type_T, levels = c("Macrophage", "cDC1 (BATF3+)", "cDC2 (CD1C+)", "LAMP3+ DC (LAMP3+)", "pDC (LILRA4+)"))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_type_T", pal=tctype_color, size = .5, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr=T) + 
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
ggsave(paste0(mplotPath, test_text, "T_ctype.pdf"), g1, width = 4, height=2.7, units = "in", device = "pdf")


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


cell_count<-as.data.frame.matrix(table(cur_proj$SampleType, cur_proj$Cell_type_T))
cell_frac <- as.data.frame(cell_count / rowSums(cell_count))
cell_frac$Sample_type <- rownames(cell_count)
plot_df <- reshape2::melt(cell_frac)
colnames(plot_df) <- c("Sample_type", "Cell_type", "Count")
plot_df$Sample_type <- factor(plot_df$Sample_type, levels = names(stype_color))
g1 = ggplot(plot_df, aes(x="", y = Count, fill = Cell_type)) + 
    geom_bar(stat = "identity", width=1) +
    coord_polar("y", start=0) + 
    facet_wrap(~Sample_type) +
    scale_fill_manual(values = tctype_color)+
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


# Stats test
cur_proj$Patient_ST <- paste0(cur_proj$SampleType, "_", cur_proj$Patient)
t_count<-as.data.frame.matrix(table(cur_proj$Patient_ST, cur_proj$Cell_type_T))
t_frac <- as.matrix(t_count / rowSums(t_count))
frac_tbl <- reshape2::melt(t_frac)
colnames(frac_tbl) <- c("Patient_ST", "T_subtype", "Fraction")
frac_tbl$SampleType <-  sapply(strsplit(as.character(frac_tbl$Patient_ST), "_"), function(x)x[1])
frac_tbl$Patient <-  sapply(strsplit(as.character(frac_tbl$Patient_ST), "_"), function(x)x[2])
frac_tbl$MSI <-  patient_meta$MSI[match(as.numeric(frac_tbl$Patient), patient_meta$Patient)]

frac_tbl <- frac_tbl[frac_tbl$MSI!="MSI",]
pres_all <- sapply(levels(frac_tbl$T_subtype), function(state){
    print(state)
    cur_tbl = frac_tbl[frac_tbl$T_subtype == state,]
    TvC <- t.test(cur_tbl$Fraction[cur_tbl$SampleType == "Tumor"], cur_tbl$Fraction[cur_tbl$SampleType == "Colon"], alternative = 'greater')
    CvT <- t.test(cur_tbl$Fraction[cur_tbl$SampleType == "Colon"], cur_tbl$Fraction[cur_tbl$SampleType == "Tumor"], alternative = 'greater')
    
    MvL <- t.test(cur_tbl$Fraction[cur_tbl$SampleType == "LiverMet"], cur_tbl$Fraction[cur_tbl$SampleType == "Liver"], alternative = 'greater')
    LvM <- t.test(cur_tbl$Fraction[cur_tbl$SampleType == "Liver"], cur_tbl$Fraction[cur_tbl$SampleType == "LiverMet"], alternative = 'greater')
    
    pres<-data.frame(TvC = TvC$p.value, CvT = CvT$p.value, MvL = MvL$p.value, LvM = LvM$p.value)
    rownames(pres) <- state
    return(pres)
})

pres_all
pres_all_star <- pres_all
pres_all_star[pres_all > 0.05] = "ns"
pres_all_star[pres_all <= 0.05] = "*"
pres_all_star[pres_all <= 0.01] = "**"
pres_all_star[pres_all <= 0.001] = "***"
pres_all_star[pres_all <= 0.0001] = "****"
pres_all[is.na(pres_all)] <- "NA"
pres_all_star[is.na(pres_all_star)] <- "NA"
write.xlsx(list(test_res = pres_all), paste0(msheetPath, test_text, "_t_sampletype_ttest_res.xlsx"), rowNames = T)
write.xlsx(list(test_res = pres_all_star), paste0(msheetPath, test_text, "_t_sampletype_ttest_res_star.xlsx"), rowNames = T)










plot_g <- c("CCR7", "FOXP3", "GZMB", "NKG7", "KLRG1", "MKI67")
plot_proj = cur_proj[cur_proj$Cell_type_T!="T-Carcinoma cell doublet",]
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
ggsave(paste0(mplotPath, test_text, "gexpr_umap.pdf"), g1, width = 4, height=3, units = "in", device = "pdf")



# T cell exhaustion
cur_vis <- clist$`T cells in vivo [20221224, IFF]`
cur_eset <- eset[, cur_vis@idx]

library(AUCell)
t_ext_markers <- list("T_exhaustion_markers"=c("PDCD1", "LAG3", "HAVCR2", "CTLA4", "CD244", "CD160", "TIGIT"))
feature_list = t_ext_markers

expressed_g <- which(rowMeans(exprs(cur_eset) > 0) > .01)
expr_hc <- as.matrix(exprs(cur_eset[expressed_g,]))
rownames(expr_hc) = fData(cur_eset)$gene_short_name[expressed_g]
cells_rankings_hc <- AUCell_buildRankings(expr_hc)
saveRDS(cells_rankings_hc, paste0(mstatePath, test_text, "aucell_exhaustion_rankings.rds"))
cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(feature_list, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "0.2cells_exhaustion_AUC.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "0.2cells_exhaustion_AUC.rds"))
matrixAUC <- getAUC(cells_AUC)
cur_eset$T_exhaustion_score = matrixAUC[1,]

cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, batch regressed]`, pData(cur_eset))

cur_proj$SampleType <- factor(cur_proj$SampleType, levels = names(stype_color))
#plot_proj <- cur_proj[!cur_proj$Cell_type_T %in% c("NK cell", "T-Carcinoma cell doublet", "KLRG1+ SLEC"),]
cur_proj <- cur_proj[!cur_proj$Cell_type_T %in% c("NK cell", "T-Carcinoma cell doublet"),]
T_subtype_mapping = c("CCR7+ Tn/Tm" = "Tn/Tm", "Cycling T cell" = "T(cycling)", "FOXP3+ Treg" = "Treg",
                      "GZMA/GZMB+ Teff/Tem" = "Teff/Tem", "KLRG1+ SLEC" = "SLEC", "NKT cell" = "NKT")
cur_proj$T_subtype = T_subtype_mapping[cur_proj$Cell_type_T]
cur_proj$T_subtype= factor(cur_proj$T_subtype, levels= T_subtype_mapping)

patient_meta <- readxl::read_excel("~/Dropbox/ColonManuscript/patient_meta.xlsx")
cur_proj$Stage = factor(patient_meta$Stage[match(as.numeric(cur_proj$Patient), patient_meta$Patient)])
cur_proj$MSI <-  patient_meta$MSI[match(as.numeric(cur_proj$Patient), patient_meta$Patient)]

g1 <- ggplot(cur_proj, aes_string(x = "T_subtype", y = "T_exhaustion_score", fill = "SampleType")) +
    geom_boxplot(outlier.colour = NA, position = position_dodge(preserve = "single"))+
    scale_fill_manual(values = stype_color)+
    scale_y_continuous(limits = c(0, .5))+
    theme_bw()+
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.4, "cm"),
          axis.text = element_text(size=8),
          #axis.text.x = element_text(angle=45, hjust=1, size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_exhaustion_stype_boxplot.pdf")),  width = 4.5, height=1.3)
g1
dev.off()


plot_g <- c("CCR7", "FOXP3", "GZMB", "NKG7", "KLRG1", "MKI67")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(plot_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- plot_g
save_proj <- cbind(cur_proj[,c("UMAP_1", "UMAP_2", "Patient", "SampleType", "Stage", "MSI", "T_subtype", "T_exhaustion_score")], g_exprs)
write.xlsx(save_proj, paste0(msheetPath, test_text, "save_proj.xlsx"), rowNames=T)


# t test
cur_proj$Signature = cur_proj$T_exhaustion_score
pres_all <- sapply(as.character(T_subtype_mapping), function(ctype){
    cur_tbl = cur_proj[cur_proj$T_subtype == ctype,]
    res_tumor <- tryCatch(t.test(cur_tbl$Signature[cur_tbl$SampleType == "Tumor"], cur_tbl$Signature[cur_tbl$SampleType == "Colon"], alternative ="greater"), error = function(err)return(data.frame(p.value = NA)))
    res_mets <- tryCatch(t.test(cur_tbl$Signature[cur_tbl$SampleType == "LiverMet"], cur_tbl$Signature[cur_tbl$SampleType == "Liver"], alternative ="greater"), error = function(err)return(data.frame(p.value = NA)))
    pres<-data.frame(tumor_vs_colon = res_tumor$p.value, mets_vs_liver = res_mets$p.value)
    rownames(pres) <- ctype
    return(pres)
})
pres_all
pres_all_star <- pres_all
pres_all_star[pres_all > 0.05] = "ns"
pres_all_star[pres_all <= 0.05] = "*"
pres_all_star[pres_all <= 0.01] = "**"
pres_all_star[pres_all <= 0.001] = "***"
pres_all_star[pres_all <= 0.0001] = "****"
pres_all[is.na(pres_all)] <- "NA"
pres_all_star[is.na(pres_all_star)] <- "NA"
write.xlsx(list(test_res = pres_all), paste0(msheetPath, test_text, "_t_exhaust_ttest_res.xlsx"), rowNames = T)
write.xlsx(list(test_res = pres_all_star), paste0(msheetPath, test_text, "_t_exhaust_ttest_res_star.xlsx"), rowNames = T)




plot_proj <- cur_proj
plot_proj <- plot_proj[plot_proj$SampleType == "Tumor" & !plot_proj$MSI == "MSI",]
use_color = get_factor_color(levels(plot_proj$Stage), "BuGn")
names(use_color) = levels(plot_proj$Stage)

g1 <- ggplot(plot_proj, aes_string(x = "T_subtype", y = "T_exhaustion_score", fill = "Stage")) +
    geom_boxplot(outlier.colour = NA, position = position_dodge(preserve = "single"))+
    scale_fill_manual(values = use_color)+
    scale_y_continuous(limits = c(0, .4))+
    theme_bw()+
    theme(text=element_text(size=8),
          axis.text = element_text(size=8),
          #axis.text.x = element_text(angle=45, hjust=1, size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_exhaustion_stage_boxplot.pdf")),  width = 5, height=1.5)
g1
dev.off()



t_proj <- cur_proj[cur_proj$SampleType == "Tumor" & !cur_proj$Cell_type_T %in% c("NK cell", "T-Carcinoma cell doublet"),]
t_count<-as.data.frame.matrix(table(t_proj$Patient, t_proj$T_subtype))
t_frac <- as.matrix(t_count / rowSums(t_count))
frac_tbl <- reshape2::melt(t_frac)
colnames(frac_tbl) <- c("Patient", "T_subtype", "Fraction")
frac_tbl$Stage <-  factor(patient_meta$Stage[match(as.numeric(frac_tbl$Patient), patient_meta$Patient)])
frac_tbl$MSI <-  patient_meta$MSI[match(as.numeric(frac_tbl$Patient), patient_meta$Patient)]

plot_proj <- frac_tbl[frac_tbl$MSI!="MSI",]
plot_proj$Stage <- factor(plot_proj$Stage)
use_color = get_factor_color(levels(plot_proj$Stage), "BuGn")
names(use_color) = levels(plot_proj$Stage)
g1 <- ggplot(plot_proj, aes_string(x = "T_subtype", y = "Fraction", fill = "Stage")) +
    geom_boxplot(outlier.colour = NA, position = position_dodge(preserve = "single"))+
    scale_fill_manual(values = use_color)+
    scale_y_continuous(limits = c(0, .4))+
    theme_bw()+
    theme(text=element_text(size=8),
          axis.text = element_text(size=8),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.4, "cm"),
          #axis.text.x = element_text(angle=45, hjust=1, size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_Tsubtype_fraction_stage_boxplot_v1.pdf")),  width = 4.2, height=1.3)
g1
dev.off()

pres_all <- sapply(levels(plot_proj$T_subtype), function(state){
    cur_tbl = plot_proj[plot_proj$T_subtype == state,]
    res_2 <- t.test(cur_tbl$Fraction[cur_tbl$Stage == 2], cur_tbl$Fraction[cur_tbl$Stage == 1])
    res_3 <- t.test(cur_tbl$Fraction[cur_tbl$Stage == 3], cur_tbl$Fraction[cur_tbl$Stage == 1])
    res_4 <- t.test(cur_tbl$Fraction[cur_tbl$Stage == 4], cur_tbl$Fraction[cur_tbl$Stage == 1])
    pres<-data.frame(II = res_2$p.value, III = res_3$p.value, IV = res_4$p.value)
    rownames(pres) <- state
    return(pres)
})


plot_proj <- frac_tbl

g1 <- ggplot(plot_proj, aes_string(x = "T_subtype", y = "Fraction", fill = "MSI")) +
    geom_boxplot(outlier.colour = NA, position = position_dodge(preserve = "single"))+
    scale_y_continuous(limits = c(0, .4))+
    theme_bw()+
    theme(text=element_text(size=8),
          axis.text = element_text(size=8),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.4, "cm"),
          #axis.text.x = element_text(angle=45, hjust=1, size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL)

pdf(paste0(paste0(mplotPath, test_text, "_Tsubtype_fraction_msi_boxplot_v1.pdf")),  width = 4, height=1.3)
g1
dev.off()


