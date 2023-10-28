


library(VisCello)
library(igraph)
savePath <- "../codex_211015_cp_cello/";dir.create(savePath)
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


test_text <- "CODEX_plot_211119_"
eset <- readRDS(paste0(savePath,"eset.rds"))
pData(eset)$cell_type == "Non_cell"
clist <- readRDS(paste0(savePath,"clist.rds"))
# Figure for codex cell type annotation
plot_entries <- names(clist)[!names(clist) == "all cells combine"]

ctype_color <- c(
    "Epithelial" = "#e31a1c",
    "Macrophage" = "#810f7c",
    "T_cell" ="#a6cee3",
    "B_cell" = "#fdbf6f",
    "Stromal_cell" = "#cab2d6", 
    "Muscle/Myofibroblast" = "#08519c",
    "Endothelial"= "#33a02c",
    "Unknown" = "lightgrey"
)
for(cur_name in plot_entries) {
    print(cur_name)
    cur_vis <- clist[[cur_name]]
    cur_proj <- cbind(cur_vis@proj$Image, pData(eset)[cur_vis@idx,])
    cur_proj <- cur_proj[cur_proj$cell_type!="Non_cell",]
    library(ggrastr)
    cur_proj$cell_type <- factor(cur_proj$cell_type, levels = names(ctype_color))
    cur_proj[2] <- -cur_proj[2]
    g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="cell_type", pal=ctype_color, size = .3, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr = T) + 
        guides(color=guide_legend(title=NULL, keyheight=.6, override.aes = list(size=3)))+
        xlab(NULL) +
        ylab(NULL) + 
        scale_x_continuous( expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(text=element_text(family = "Helvetica", size=8),
              legend.text=element_text(size=8),
              axis.text = element_text(size=9),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white", colour = "black", size = .2),
              #legend.box.margin=margin(-10,-10,-10,-10),
              plot.margin = unit(c(0,0.4,0,0.1), "cm"))
    ggsave(paste0(mplotPath, test_text, cur_name, "_img_cell_type.pdf"), g1, width = 4, height=1.8, units = "in", device = "pdf")
}



# Figure for fraction of cells

pmeta_filtered <- pData(eset)[pData(eset)$cell_type != "Non_cell",]

compos<-as.data.frame.matrix(table(pmeta_filtered$dataset,pmeta_filtered$cell_type))
show_ctype <- names(ctype_color)[names(ctype_color) != "Unknown"]
compos_show <- compos[,show_ctype]
compos_show <- compos_show / rowSums(compos_show)
pdf(paste0(mplotPath, test_text, "compos_heatmap", ".pdf"), width = 3.5, height = 2.5)
pheatmap(compos_show, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "red"))(10))
dev.off()

 

# Pairwise cell type enrichment
res_list <- list()
for(cur_name in plot_entries) {
    print(cur_name)
    cur_vis <- clist[[cur_name]]
    cur_proj <- cbind(cur_vis@proj$Image, pData(eset)[cur_vis@idx,])
    cur_proj <- cur_proj[!cur_proj$cell_type%in% c("Non_cell", "Unknown"),]
    library(ggrastr)
    cur_proj$cell_type <- factor(cur_proj$cell_type, levels = names(ctype_color))
    
    knn_k = 6
    knn_res <-FNN::get.knn(cur_proj[,c("X","Y")],k=knn_k)
    knn_ctype <- apply(knn_res$nn.index, 2, function(x) cur_proj$cell_type[x])
    
    sur_res<-sapply(show_ctype, function(ct) {
        suroundings <- knn_ctype[cur_proj$cell_type == ct,]
        table(unlist(suroundings))/length(unlist(suroundings))
    })
    colnames(sur_res) <-show_ctype
    sur_res <- t(sur_res[colnames(sur_res),])
    res_list[[cur_name]] <- sur_res
    pdf(paste0(mplotPath, test_text, cur_name, "_avg_surrounding_",knn_k,"nn.pdf"), width = 3, height = 3)
    pheatmap(sur_res, cluster_cols = F, cluster_rows = F, display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "red"))(10), gaps_row = 1:nrow(sur_res), legend = F)
    dev.off()
}

normal_average <- (res_list$CRC36_86_reg001 + res_list$CRC37_92_reg001 + res_list$CRC38_28_reg001)/ 3
tumor_average <- (res_list$CRC36_86_reg002 + res_list$CRC37_92_reg002 + res_list$CRC38_28_reg002 + res_list$CRC38_28_reg003)/ 4

pdf(paste0(mplotPath, test_text, "normal", "_avg_surrounding_",knn_k,"nn.pdf"), width = 3, height = 3)
pheatmap(normal_average, cluster_cols = F, cluster_rows = F, display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#3690c0"))(10), gaps_row = 1:nrow(normal_average), legend = F)
dev.off()

pdf(paste0(mplotPath, test_text, "tumor", "_avg_surrounding_",knn_k,"nn.pdf"), width = 3, height = 3)
pheatmap(tumor_average, cluster_cols = F, cluster_rows = F, display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(10), gaps_row = 1:nrow(tumor_average), legend = F)
dev.off()

write.xlsx(res_list, paste0(msheetPath, test_text, "avg_surrounding_6nn.xlsx"), rowNames = T)

# cur_ctype <- show_ctype[1]
# cur_sur <- rbind(data.frame(fraction = c(res_list$CRC36_86_reg001[cur_ctype, ], res_list$CRC37_92_reg001[cur_ctype, ], res_list$CRC38_28_reg001[cur_ctype, ]), cell_type = rep(show_ctype,3), sample_type="normal"), data.frame(fraction = c(res_list$CRC36_86_reg002[cur_ctype, ], res_list$CRC37_92_reg002[cur_ctype, ], res_list$CRC38_28_reg002[cur_ctype, ], res_list$CRC38_28_reg003[cur_ctype, ]), cell_type = rep(show_ctype,4), sample_type="tumor"))
# 
# g1<-ggplot(cur_sur, aes(x=cell_type, y=fraction,fill=sample_type)) + 
#     geom_boxplot(lwd=0.3, outlier.size = 0.1, show.legend = F) + 
#     geom_point(data = df_med, aes(x=Dataset, y=med_g), size = 0, stroke=0, shape = 22) +
#     #scale_fill_manual(values = combine_dataset_color) +
#     scale_y_continuous(limits = c(0, 8000))+
#     labs(x="", y = "#Genes/cell")+
#     guides(fill=F)+
#     theme_bw() + 
#     theme(text = element_text(size=6), 
#           legend.text=element_text(size=6),
#           axis.text = element_text(size=6),
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#           legend.position = "right", 
#           plot.title = element_text(hjust = 0.5),
#           plot.margin = unit(c(0,0.4,0,0.1), "cm"))


# Macrophage subtypes
hist(eset@assayData$norm_exprs["CD68-BX15_cyc004_ch004",])
#hist(eset@assayData$norm_exprs["CD14-BX49_cyc015_ch002",])
hist(eset@assayData$norm_exprs["SPP1-BX33_cyc006_ch004",])
#hist(eset@assayData$norm_exprs["FIBRONECTIN-BX31_cyc016_ch002",])
#hist(eset@assayData$norm_exprs["MARCO-BX04_cyc003_ch002",])
hist(eset@assayData$norm_exprs["C1QB-BX37_cyc011_ch002",])
hist(eset@assayData$norm_exprs["HLA-DQA1-BX32_cyc012_ch003",])

eset$mac_state = "Non-macrophage"
eset$mac_state[eset$cell_type == "Macrophage"] = "Macrophage"
# C1QB_filter <- eset@assayData$norm_exprs["C1QB-BX37_cyc011_ch002",] > 6 & eset$cell_type == "Macrophage"
# eset$mac_state[C1QB_filter] = paste0(eset$mac_state[C1QB_filter], ", C1QB")
HLADQA1_filter <- eset@assayData$norm_exprs["HLA-DQA1-BX32_cyc012_ch003",] > 8 & eset$cell_type == "Macrophage"
eset$mac_state[HLADQA1_filter] = paste0(eset$mac_state[HLADQA1_filter], ", HLA-DQA1")
SPP1_filter <- eset@assayData$norm_exprs["SPP1-BX33_cyc006_ch004",] > 8 & eset$cell_type == "Macrophage"
eset$mac_state[SPP1_filter] = paste0(eset$mac_state[SPP1_filter], ", SPP1")
#saveRDS(eset, paste0(savePath, "eset.rds"))
table(eset$mac_state, eset$dataset)


mac_meta <- pData(eset)[eset$cell_type == "Macrophage",]
mac_meta$sample_type <- ifelse(grepl("reg001", mac_meta$dataset), "Normal", "Tumor")
mac_compos <- as.data.frame.matrix(table(mac_meta$mac_state, mac_meta$sample_type))

mac_frac <- t(t(mac_compos) / colSums(mac_compos))
mac_frac_plot <- reshape2::melt(mac_frac)

mac_frac_plot <- mac_frac_plot[mac_frac_plot$Var1 != "Macrophage",]
mac_frac_plot$Var1 <- gsub("Macrophage, ", "", mac_frac_plot$Var1)
library(ggbreak) 
g1 <- ggplot(mac_frac_plot, aes(x=Var1, y=value,fill=Var2)) + 
    geom_bar(stat = "identity",position = "dodge2") + 
    scale_fill_manual(values = c("Normal" = "#1f78b4", "Tumor" = "#e31a1c")) +
    #scale_y_continuous(expand = c(0, 0)) + 
    labs(x="", y = "Fraction")+
    theme(text = element_text(size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=8),
          #axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) + 
    scale_y_break(c(0.02, 0.45)) +
    monocle:::monocle_theme_opts()
ggsave(paste0(mplotPath, test_text, "_macsubstate_ybreak.pdf"), g1, width = 3.5, height=2, units = "in", device = "pdf")

write.xlsx(mac_compos, paste0(msheetPath, test_text, "mac_compos.xlsx"), rowNames = T)


prop.test(mac_compos["Macrophage, SPP1", "Tumor"], sum(mac_compos[, "Tumor"]), p = mac_compos["Macrophage, SPP1", "Normal"]/sum(mac_compos[, "Normal"]), alternative = "two.sided")

prop.test(mac_compos["Macrophage, HLA-DQA1, SPP1", "Tumor"], sum(mac_compos[, "Tumor"]), p = mac_compos["Macrophage, HLA-DQA1, SPP1", "Normal"]/sum(mac_compos[, "Normal"]), alternative = "two.sided")

prop.test(mac_compos["Macrophage, HLA-DQA1", "Tumor"], sum(mac_compos[, "Tumor"]), p = mac_compos["Macrophage, HLA-DQA1", "Normal"]/sum(mac_compos[, "Normal"]), alternative = "two.sided")

