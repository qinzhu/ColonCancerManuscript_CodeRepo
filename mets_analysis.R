

library(VisCello)

mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

scriptPath <- paste0("../hcc_cello/scripts/")
source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")


# Liver mets analysis

test_text <- "mets_invivo_0421_"
cur_vis <- clist$`All in vivo cells [20210411, IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)
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
use_dataset <- c("Colon_09","Tumor_09", "Tumor_27", "Tumor_44", "Tumor_49", "Liver_09", "Liver_44", "LiverMet_09", "LiverMet_27", "LiverMet_44", "LiverMet_49")
#use_dataset <- c("Colon_09", "Colon_44", "Tumor_09", "Tumor_27", "Tumor_44", "LiverMet_09", "LiverMet_27", "LiverMet_44", "Liver_09", "Liver_44")
#use_dataset <- c("Tumor_09", "Tumor_27",  "LiverMet_09", "LiverMet_27", "Liver_09", "Liver_44")
cur_proj <- cur_proj[cur_proj$Dataset %in% use_dataset,]

cur_proj$Cell_type <- factor(cur_proj$Cell_type, levels = names(ctype_color))
cur_proj$Sample_Type <- factor(cur_proj$SampleType, levels = names(stype_color))
compos<-as.data.frame.matrix(table(cur_proj$Sample_Type, cur_proj$Cell_type))
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

# compos by fraction heatmap
compos_frac <- compos / rowSums(compos)
pdf(paste0(mplotPath, test_text, "compos_frac_heatmap", ".pdf"), width = 4, height = 1.88)
pheatmap(compos_frac[,show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "red"))(10))
dev.off()

compos_rel <- compos / (compos$`Epithelial(Normal)` + compos$`Epithelial(Tumor)`)
pdf(paste0(mplotPath, test_text, "compos_rel_heatmap", ".pdf"), width = 3.5, height = 1.58)
pheatmap(compos_rel[complete.cases(compos_rel),show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "red"))(10))
dev.off()


write.xlsx(compos, paste0(msheetPath,test_text,"compos.xlsx"), rowNames=T)


cur_proj$Dataset <- factor(cur_proj$Dataset, levels = use_dataset)
compos<-as.data.frame.matrix(table(cur_proj$Dataset, cur_proj$Cell_type))
compos_frac <- compos / rowSums(compos)
compos_rel <- compos / (compos$`Epithelial(Normal)` + compos$`Epithelial(Tumor)`)
breaksList <- seq(0,1,by=.1)
pdf(paste0(mplotPath, test_text, "compos_dataset_heatmap", ".pdf"), width = 3.8, height = 3.2)
pheatmap(compos_frac[,show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)
#pheatmap(compos_rel[complete.cases(compos_rel),show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "red"))(10))
dev.off()

# Entropy

entropy_df <- data.frame(entropy = apply(compos_frac[,show_ctype],1, function(freqs) {
    -sum(freqs * log2(freqs), na.rm = T)
}))
entropy_df$dataset <- rownames(entropy_df)
entropy_df$dataset <- factor(entropy_df$dataset, levels = rev(use_dataset))
entropy_df$SampleType <- factor(sapply(strsplit(rownames(entropy_df), "_"), function(x)x[1]), levels = names(stype_color))
g1<-ggplot(entropy_df, aes(x = dataset, y = entropy)) +
    geom_bar(stat="identity", aes(fill = SampleType)) + 
    scale_fill_manual(values = stype_color)+
    coord_flip() + 
    #scale_y_continuous(position = "right") +
    theme_classic()
ggsave(paste0(mplotPath, test_text, "entropy_bar.pdf"), g1, width = 3, height=3, units = "in", device = "pdf")

# 
t.test(entropy_df$entropy[grepl("Liver_", entropy_df$dataset)], entropy_df$entropy[grepl("LiverMet_", entropy_df$dataset)], alternative = "less") # 0.05266



# test_text <- "mets_epi_"
# library(AUCell)
# cur_vis <- clist$`Metastasis epithelial invivo [cleaned]`
# cur_eset <- eset[, cur_vis@idx]
# 
# # Test on DEG using Mann-Whitny U test
# mets_09_deg<-read_excel_allsheets(paste0(msheetPath, "Metastasis epithelial invivo [cleaned]_LiverMet_09_vs_Tumor_09_2020-08-12_de_significant.xlsx"))
# mets_27_deg<-read_excel_allsheets(paste0(msheetPath, "Metastasis epithelial invivo [cleaned]_LiverMet_27_vs_Tumor_27_2020-08-12_de_significant.xlsx"))
# mets_09_deg <- mets_09_deg$LiverMet_09 # 165
# mets_27_deg <- mets_27_deg$LiverMet_27 # 7845
# 
# # mets_09_sSeq<-read_excel_allsheets(paste0(msheetPath, "Metastasis epithelial invivo [cleaned]_LiverMet_09_vs_Tumor_09_2020-08-12_de_sSeq.xlsx"))
# # mets_27_sSeq<-read_excel_allsheets(paste0(msheetPath, "Metastasis epithelial invivo [cleaned]_LiverMet_27_vs_Tumor_27_2020-08-12_de_sSeq.xlsx"))
# 
# shared_deg <- intersect(mets_09_deg$gene_short_name, mets_27_deg$gene_short_name)
# non_epi_sig<-readRDS(paste0(mstatePath,"de_non_epi_primary_","de_sig_genes.rds"))
# b_cell_genes <- non_epi_sig$`Plasma cell`$gene_name
# shared_deg <- shared_deg[!shared_deg%in% as.character(b_cell_genes)]
# shared_deg_res <- data.frame(gene = shared_deg, 
#                              mets_09_FDR = mets_09_deg$FDR[match(shared_deg, mets_09_deg$gene_short_name)],
#                              mets_09_effectSize = mets_09_deg$effectSize[match(shared_deg, mets_09_deg$gene_short_name)], 
#                              mets_27_FDR = mets_27_deg$FDR[match(shared_deg, mets_27_deg$gene_short_name)],
#                              mets_27_effectSize = mets_27_deg$effectSize[match(shared_deg, mets_27_deg$gene_short_name)])
# 
# write.xlsx(shared_deg_res, paste0(msheetPath, test_text, "shared_deg_res.xlsx"))
# 
# 
# 
# # heatmap of DEG
# 
# max_cell = 1000
# max_g = 100
# plot_g <- rownames(fData(cur_eset))[match(shared_deg, fData(cur_eset)$gene_short_name)]
# use_sample <- c("Tumor_09", "Tumor_27", "LiverMet_09", "LiverMet_27")
# plot_sample <- unlist(lapply(use_sample, function(x) {
#     cur_idx <- which(cur_eset$Dataset == x)
#     sample(cur_idx, min(max_cell, length(cur_idx)), replace = F)
# }))
# eset_plot <- cur_eset[plot_g, plot_sample]
# eset_plot <- eset_plot[rowMeans(exprs(eset_plot) > 0) > .1, ]
# dim(eset_plot)
# plot_meta <- pData(eset_plot)
# plot_meta <- plot_meta[, c("Dataset"), drop=F]
# use_color = c(
#     "Tumor_09" = "#fb9a99",   
#     "Tumor_27" = "#fdbf6f",   
#     "LiverMet_09" = "#e31a1c", 
#     "LiverMet_27" = "#ff7f00"
# )
# plot_meta$Dataset <- factor(plot_meta$Dataset, levels = names(use_color))
# plot_meta <- plot_meta[order(plot_meta$Dataset),, drop=F]
# non_na_cells_ordered <- rownames(plot_meta)
# value <- eset_plot@assayData$norm_exprs[,non_na_cells_ordered]
# rownames(value) <- fData(eset_plot)$gene_short_name[match(rownames(value), rownames(eset_plot))]
# 
# value <- t(scale(t(as.matrix(value))))
# limits <- c(-2,2)
# if(!is.null(limits)) {
#     value[value < limits[1]] <- limits[1]
#     value[value > limits[2]] <- limits[2]
# }
# 
# pdf(paste0(mplotPath, test_text, "deg_mets", ".pdf"), width = 5, height=7)
# pheatmap::pheatmap(value,
#                    color = get_numeric_color("viridis"),
#                    cluster_rows = T, cluster_cols = F,
#                    #clustering_distance_rows = distfun1, 
#                    clustering_method = "ward.D2",
#                    show_rownames = T,
#                    show_colnames = FALSE, annotation_row = NULL,
#                    annotation_col = plot_meta, annotation_names_row = FALSE, annotation_colors = list(Dataset = use_color),
#                    annotation_names_col = FALSE, 
#                    fontsize = 6)
# dev.off()

# marker_gene <- c("CEACAM5", "CEACAM1", "SATB2", "TP53", "BRAF","PIK3CA", "RECQL4")
# marker_gene %in% fData(cur_eset)$gene_short_name
# 
# eset_plot <- cur_eset[, cur_eset$Dataset %in%use_sample]
# 
# eset_09 <- eset_plot[fData(eset_plot)$gene_short_name %in% marker_gene, eset_plot$Patient == "09"]
# de_09 <- run_mw(dat=as.matrix(exprs(eset_09)), group=eset_09$Dataset, fdata = fData(eset_09), min_fdr= 1, id_col = "gene_id", name_col = "gene_short_name",  detRate = 0)
# 
# eset_27 <- eset_plot[fData(eset_plot)$gene_short_name %in% marker_gene, eset_plot$Patient == "27"]
# de_27 <- run_mw(dat=as.matrix(exprs(eset_27)), group=eset_27$Dataset, fdata = fData(eset_27), min_fdr= 1, id_col = "gene_id", name_col = "gene_short_name",  detRate = 0)
# 
# g_exprs <- as.data.frame(t(as.matrix(eset_plot@assayData$norm_exprs[match(marker_gene,fData(eset_plot)$gene_short_name),])))
# colnames(g_exprs) <- marker_gene
# plot_meta <- pData(eset_plot)
# use_color = c(
#     "Tumor_09" = "#fb9a99",   
#     "LiverMet_09" = "#e31a1c",
#     "Tumor_27" = "#fdbf6f",   
#     "LiverMet_27" = "#ff7f00"
# )
# plot_meta$Dataset <- factor(plot_meta$Dataset, levels = names(use_color))
# 
# glist <- list()
# for(g in marker_gene) {
#     print(paste0(g, ", 09:", g %in% mets_09_deg$gene_short_name, "; ", "27:", g %in% mets_27_deg$gene_short_name))
#     gene_values <- as.matrix(g_exprs[,g,drop=F])
#     
#     # Run MW test
#     
#     ecut <- quantile(gene_values[gene_values > 0], .99)
#     #gene_values[gene_values == 0] <- NA
#     noise <- rnorm(n = length(x = gene_values[,g]))/1e+05
#     gene_values[,g] <- gene_values[,g] + noise
#     #ecut <- max(quantile(gene_values,1),1)
#     bpGroup = "Dataset"
#     colnames(gene_values) <- "expression_level"
#     df <- cbind(gene_values, plot_meta)
#     glist[[g]]<- ggplot(df, aes_string(x=bpGroup, y="expression_level")) +
#         geom_violin(aes_string(fill = bpGroup, color= bpGroup), trim = T, scale = "width", alpha = .7) + 
#         geom_jitter(size = .3, aes_string(colour = bpGroup), alpha = 1)+
#         scale_color_manual(values = use_color) +
#         scale_fill_manual(values = use_color)+
#         #ylim(0, ecut) +
#         theme_classic() + 
#         guides(alpha = F, fill=F,color=F) + 
#         xlab("") +
#         ylab(g) +
#         theme(text=element_text(family = "Helvetica", size=8),
#               axis.text.x = element_text(angle=45, hjust=1, size=8), 
#               axis.text.y = element_text(size=8), 
#               legend.text=element_text(size=8),
#               axis.text = element_text(size=8),
#               legend.margin=margin(0,0,0,0))
# }
# 
# g1<-do.call(grid.arrange,c(glist, ncol = 3))
# ggsave(paste0(mplotPath, test_text, "violin_gexpr_violin.pdf"), g1, width = 5, height=7, units = "in", device = "pdf")
# 
# 
# 
# expressed_g <- which(rowMeans(exprs(cur_eset) > 0) > .1)
# expressed_g <- expressed_g[!fData(cur_eset)$gene_short_name[expressed_g] %in% as.character(b_cell_genes)]
# expr_hc <- as.matrix(exprs(cur_eset[expressed_g,]))
# rownames(expr_hc) <- make.names(fData(cur_eset)$gene_short_name[match(rownames(expr_hc), fData(cur_eset)$gene_id)])
# cells_rankings_hc <- AUCell_buildRankings(expr_hc)
# #reactogsea_list <- readRDS(paste0(rdsPath, "hc_dissect_nohscc_terminal_terminal_reactogsea_gsea.rds"))
# reactomeSets <- gage::readList("data-raw/public_resource/ReactomePathways.gmt")
# reactomeSets_genes <- sapply(reactomeSets, function(x){paste(x, collapse = "/")})
# cur_thresh <- 0.2
# cells_AUC <- AUCell_calcAUC(reactomeSets, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
# #saveRDS(cells_AUC, paste0(rdsPath, test_text, cur_thresh, "cells_AUC", "_wt_ifg_cut75.rds"))
# saveRDS(cells_AUC, paste0(mstatePath, test_text, "aucell_reactome_0.2cells_AUC.rds"))
# cells_AUC <- readRDS(paste0(mstatePath, test_text,"aucell_reactome_0.2cells_AUC.rds"))
# matrixAUC <- getAUC(cells_AUC)
# 
# 
# 
# test_folder <- "mets_pathway/"; dir.create(paste0(mstatePath, test_folder))
# matrixAUC_test <- matrixAUC
# 
# compare_pairs <- list("Patient_09" = c("LiverMet_09", "Tumor_09"), "Patient_27" = c("LiverMet_27","Tumor_27"))
# cluster_test <- cur_eset$Dataset
# for(i in 1:length(compare_pairs)) {
#     p <- compare_pairs[[i]]
#     pn <- names(compare_pairs)[i]
#     print(pn)
#     g1_idx <- which(cluster_test == p[1])
#     g2_idx <- which(cluster_test == p[2])
#     tres<-apply(matrixAUC_test, 1, function(x){
#         g1_val <- x[g1_idx]
#         g2_val <- x[g2_idx]
#         t.test(g1_val, g2_val)
#     })
#     tres_tbl<- purrr::map_df(tres, broom::tidy)
#     tres_tbl <- tres_tbl %>% 
#         tibble::add_column(reactome_pathway= rownames(matrixAUC_test), .before=1) %>%
#         tibble::add_column(q.value = p.adjust(tres_tbl$p.value, method = "fdr"), .before=7) %>%
#         tibble::add_column(pathway_genes = reactomeSets_genes[rownames(matrixAUC_test)]) %>%
#         dplyr::arrange(q.value)
#     tres_tbl_list <- list(
#         tres_tbl %>% dplyr::filter(statistic > 0),
#         tres_tbl %>% dplyr::filter(statistic < 0)
#     )
#     names(tres_tbl_list) <- p
#     saveRDS(tres_tbl_list, paste0(mstatePath, test_folder, pn, "ttest_pathway.rds"))
#     write.xlsx(tres_tbl_list, paste0(mstatePath, test_folder, pn, "ttest_pathway.xlsx"))
# }
# 
# 
# mets_pathway_list <- list()
# q_cut <- 1e-2
# for(i in 1:length(compare_pairs)) {
#     pn <- names(compare_pairs)[i]
#     tres_tbl_list <- readRDS(paste0(mstatePath, test_folder, pn, "ttest_pathway.rds"))
#     mets_pathway_list[[pn]] <- tres_tbl_list[[grep("LiverMet", names(tres_tbl_list))]]
# }
# 
# mets_pathways<- lapply(mets_pathway_list, function(x) {
#     x$reactome_pathway[x$q.value <= q_cut]
# })
# # Shared pathway
# shared_p <- Reduce(intersect, mets_pathways)
# 
# specific_deg_list<-lapply(list(shared = shared_p), function(x) {
#     res<-reactomeSets[x]
#     res[sapply(res, length) >=5]
# })
# 
# remove_redundancy <- function(deg_list, red_cut=.5) {
#     return_list <- lapply(deg_list, function(x) {
#         for(i in 1:(length(x)-1)){
#             if(!is.na(x[[i]])) {
#                 for(j in (i+1):length(x)) {
#                     if(!is.na(x[[j]])) {
#                         x_i <- unlist(x[[i]])
#                         x_j <- unlist(x[[j]])
#                         res<-length(intersect(x[[i]], x[[j]]))/length(union(x_i, x_j)) 
#                         if(res >= red_cut) x[[j]] <- NA
#                     }
#                 }  
#             }
#         }
#         return(x)
#     })
#     return_list <-lapply(return_list, function(x) names(x)[!is.na(x)])
#     return(return_list)
# }
# plot_list<-remove_redundancy(specific_deg_list, red_cut = .3)
# 
# use_color = c(
#     "Tumor_09" = "#fb9a99",   
#     "Tumor_27" = "#fdbf6f", 
#     "LiverMet_09" = "#e31a1c",
#     "LiverMet_27" = "#ff7f00"
# )
# plot_meta <-  pData(cur_eset)[colnames(matrixAUC_test),]
# plot_meta <- plot_meta[plot_meta$Dataset %in% names(use_color),]
# plot_meta$Dataset <- factor(plot_meta$Dataset, levels = names(use_color))
# plot_meta <- plot_meta[order(plot_meta$Dataset),]
# 
# plotp <- plot_list$shared
# value  <- matrixAUC[plotp,rownames(plot_meta)]
# value <- t(scale(t(as.matrix(value))))
# limits <- c(-1.5,1.5)
# if(!is.null(limits)) {
#     value[value < limits[1]] <- limits[1]
#     value[value > limits[2]] <- limits[2]
# }
# 
# 
# pdf(paste0(mplotPath, test_text, "upinmets_pathway_", "heatmap",".pdf"), width = 8, height=2.5)
# pheatmap::pheatmap(value,
#                    color = get_numeric_color("RdBu"),
#                    cluster_rows = T, cluster_cols = F,
#                    clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
#                    show_colnames = FALSE, annotation_row = NULL,
#                    annotation_col = plot_meta[,c("Dataset"),drop=F], annotation_names_row = FALSE,
#                    annotation_names_col = FALSE, 
#                    annotation_colors = list("Dataset" = use_color), 
#                    fontsize = 8)
# dev.off()



# Epithelial DE analysis compare Mets vs Tumor
#mets_epi_deg<-read_excel_allsheets(paste0(msheetPath, "Met_epi_cleaned_LiverMet_vs_Tumor_2021-04-22_de_significant.xlsx"))
mets_nonepi_deg<- read_excel_allsheets(paste0(msheetPath, "met_only_Myeloid_vs_Plasma cell_vs_T cell_vs_B cell_vs_Mast cell_vs_Fibroblast_vs_Endothelial_vs_Myofibroblast_2021-04-22_de_significant.xlsx"))


# # Test on DEG using Mann-Whitny U test
mets_09_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_09_vs_Tumor_09_2021-04-22_de_significant.xlsx"))
mets_27_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_27_vs_Tumor_27_2021-04-22_de_significant.xlsx"))
mets_44_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_44_vs_Tumor_44_2021-04-22_de_significant.xlsx"))
mets_49_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_49_vs_Tumor_49_2021-04-22_de_significant.xlsx"))

mets_epi_deg<- list(p_09 = mets_09_deg$LiverMet_09, p_27 = mets_27_deg$LiverMet_27, p_44 = mets_44_deg$LiverMet_44, p_49 = mets_49_deg$LiverMet_49)

# Remove contaminating
mets_epi_deg <- lapply(mets_epi_deg, function(x) {
    x$gene_id[!x$gene_id %in% mets_nonepi_deg$`Plasma cell`$gene_id]
})
# Get consensus genes
union_deg<- as.data.frame(table(unlist(mets_epi_deg)))
colnames(union_deg) <- c("gene_id", "Freq")
union_deg$gene_short_name <- fData(eset)$gene_short_name[match(union_deg$gene_id, fData(eset)$gene_id)]
union_deg<- union_deg[order(union_deg$Freq, decreasing = T),]
union_deg$fdr_09 <- mets_09_deg$LiverMet_09$FDR[match(union_deg$gene_id, mets_09_deg$LiverMet_09$gene_id)]
union_deg$fdr_27 <- mets_27_deg$LiverMet_27$FDR[match(union_deg$gene_id, mets_27_deg$LiverMet_27$gene_id)]
union_deg$fdr_44 <- mets_44_deg$LiverMet_44$FDR[match(union_deg$gene_id, mets_44_deg$LiverMet_44$gene_id)]
union_deg$fdr_49 <- mets_49_deg$LiverMet_49$FDR[match(union_deg$gene_id, mets_49_deg$LiverMet_49$gene_id)]
write.csv(union_deg, paste0(msheetPath, "mets_epi_deg_union.csv"))

cur_vis <- clist$`All in vivo cells [20210411, IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
use_dataset <- c("Tumor_09", "Tumor_27", "Tumor_44", "Tumor_49", "LiverMet_09", "LiverMet_27", "LiverMet_44", "LiverMet_49")
cur_proj <- cur_proj[cur_proj$Dataset %in% use_dataset & cur_proj$Cell_type %in% c("Epithelial(Normal)", "Epithelial(Tumor)"),]

library(ggpubr)
plot_df <- cur_proj
test_genes <- c("ZNF207", "DDX5", "H3F3B", "BTG1", "IL2RG", "CEACAM1")
eset_plot <- eset[,rownames(plot_df)]
#saveRDS(eset_plot, paste0(mstatePath, "mets_epi_eset_plot.rds"))
g_exprs <- as.data.frame(t(as.matrix(eset_plot@assayData$norm_exprs[match(test_genes,fData(eset_plot)$gene_short_name),])))
colnames(g_exprs) <- test_genes
plot_meta <- pData(eset_plot)
plot_meta$Dataset <- factor(plot_meta$Dataset, levels = use_dataset)
glist <- list()
for(g in test_genes) {
    gene_values <- as.matrix(g_exprs[,g,drop=F])
    ecut <- quantile(gene_values[gene_values > 0], .99)
    noise <- rnorm(n = length(x = gene_values[,g]))/1e+05
    gene_values[,g] <- gene_values[,g] + noise
    colnames(gene_values) <- "expression_level"
    df <- cbind(gene_values, plot_meta)
    glist[[g]]<- ggplot(df, aes_string(x="Patient", y="expression_level")) +
        geom_violin(aes_string(fill = "SampleType", color= "SampleType"), trim = T, scale = "width", alpha = .7) +
        geom_jitter(size = .5, stroke=0, aes_string(colour = "SampleType"), alpha = 1, position=position_jitterdodge())+
        scale_color_manual(values = stype_color) +
        scale_fill_manual(values = stype_color)+
        #ylim(0, ecut) +
        theme_classic() +
        ggtitle(g)+
        #guides(alpha = F, fill=F,color=F) +
        xlab("") +
        ylab("") +
        theme(text=element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(angle=45, hjust=1, size=8),
              axis.text.y = element_text(size=8),
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              plot.title = element_text(hjust = 0.5),
              legend.margin=margin(0,0,0,0))
}

g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "violin_gexpr_violin.pdf"), g1, width = 5, height=3, units = "in", device = "pdf")


save_df <- cbind(g_exprs, plot_meta[, c("Patient", "SampleType")])
write.xlsx(save_df, paste0(msheetPath,test_text,"violin_gexpr.xlsx"), rowNames=T)



library(AUCell)
eset_plot <- readRDS(paste0(mstatePath, "mets_epi_eset_plot.rds"))
b_cell_genes <- mets_nonepi_deg$`Plasma cell`$gene_short_name
cur_eset <- eset_plot
expressed_g <- which(rowMeans(exprs(cur_eset) > 0) > .05)
expressed_g <- expressed_g[!fData(cur_eset)$gene_short_name[expressed_g] %in% as.character(b_cell_genes)]
expr_hc <- as.matrix(exprs(cur_eset[expressed_g,]))
rownames(expr_hc) <- make.names(fData(cur_eset)$gene_short_name[match(rownames(expr_hc), fData(cur_eset)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_hc)
#reactogsea_list <- readRDS(paste0(rdsPath, "hc_dissect_nohscc_terminal_terminal_reactogsea_gsea.rds"))
reactomeSets <- gage::readList("data-raw/public_resource/ReactomePathways.gmt")
reactomeSets_genes <- sapply(reactomeSets, function(x){paste(x, collapse = "/")})
cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(reactomeSets, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
#saveRDS(cells_AUC, paste0(rdsPath, test_text, cur_thresh, "cells_AUC", "_wt_ifg_cut75.rds"))
saveRDS(cells_AUC, paste0(mstatePath, test_text, "aucell_reactome_0.2cells_AUC.rds"))
cells_AUC <- readRDS(paste0(mstatePath, test_text,"aucell_reactome_0.2cells_AUC.rds"))
matrixAUC <- getAUC(cells_AUC)

write.xlsx(as.data.frame(matrixAUC), paste0(msheetPath,test_text, "aucell_reactome_0.2cells_AUC.xlsx"), rowNames = T)

test_folder <- "mets_pathway_0422/"; dir.create(paste0(mstatePath, test_folder))
matrixAUC_test <- matrixAUC

compare_pairs <- list(
    "Patient_09" = c("LiverMet_09", "Tumor_09"), 
    "Patient_27" = c("LiverMet_27","Tumor_27"),
    "Patient_44" = c("LiverMet_44","Tumor_44"), 
    "Patient_49" = c("LiverMet_49","Tumor_49"))
cluster_test <- cur_eset$Dataset
for(i in 1:length(compare_pairs)) {
    p <- compare_pairs[[i]]
    pn <- names(compare_pairs)[i]
    print(pn)
    g1_idx <- which(cluster_test == p[1])
    g2_idx <- which(cluster_test == p[2])
    tres<-apply(matrixAUC_test, 1, function(x){
        g1_val <- x[g1_idx]
        g2_val <- x[g2_idx]
        t.test(g1_val, g2_val)
    })
    tres_tbl<- purrr::map_df(tres, broom::tidy)
    tres_tbl <- tres_tbl %>%
        tibble::add_column(reactome_pathway= rownames(matrixAUC_test), .before=1) %>%
        tibble::add_column(q.value = p.adjust(tres_tbl$p.value, method = "fdr"), .before=7) %>%
        tibble::add_column(pathway_genes = reactomeSets_genes[rownames(matrixAUC_test)]) %>%
        dplyr::arrange(q.value)
    tres_tbl_list <- list(
        tres_tbl %>% dplyr::filter(statistic > 0),
        tres_tbl %>% dplyr::filter(statistic < 0)
    )
    names(tres_tbl_list) <- p
    saveRDS(tres_tbl_list, paste0(mstatePath, test_folder, pn, "ttest_pathway.rds"))
    write.xlsx(tres_tbl_list, paste0(mstatePath, test_folder, pn, "ttest_pathway.xlsx"))
}


mets_up_list <- list()
mets_down_list <- list()
q_cut <- 1e-2
for(i in 1:length(compare_pairs)) {
    pn <- names(compare_pairs)[i]
    tres_tbl_list <- readRDS(paste0(mstatePath, test_folder, pn, "ttest_pathway.rds"))
    mets_up_list[[pn]] <- tres_tbl_list[[grep("LiverMet", names(tres_tbl_list))]]
    mets_down_list[[pn]] <- tres_tbl_list[[grep("Tumor", names(tres_tbl_list))]]
}

mets_up<- lapply(mets_up_list, function(x) {
    x$reactome_pathway[x$q.value <= q_cut]
})
mets_down <- lapply(mets_down_list, function(x) {
    x$reactome_pathway[x$q.value <= q_cut]
})
# Shared pathway
union_up <- as.data.frame(table(unlist(mets_up)))
union_up$mlog_qval_09 <- -log10(mets_up_list$Patient_09$q.value[match(union_up$Var1, mets_up_list$Patient_09$reactome_pathway)])
union_up$mlog_qval_27 <- -log10(mets_up_list$Patient_27$q.value[match(union_up$Var1, mets_up_list$Patient_27$reactome_pathway)])
union_up$mlog_qval_44 <- -log10(mets_up_list$Patient_44$q.value[match(union_up$Var1, mets_up_list$Patient_44$reactome_pathway)])
union_up$mlog_qval_49 <- -log10(mets_up_list$Patient_49$q.value[match(union_up$Var1, mets_up_list$Patient_49$reactome_pathway)])
union_up$mlog_qval_avg <- rowMeans(union_up[,c("mlog_qval_09", "mlog_qval_27", "mlog_qval_44", "mlog_qval_49")], na.rm = T)
union_up <- union_up[order(union_up$Freq, union_up$mlog_qval_avg, decreasing = T),]


union_down <- as.data.frame(table(unlist(mets_down)))
union_down$mlog_qval_09 <- -log10(mets_down_list$Patient_09$q.value[match(union_down$Var1, mets_down_list$Patient_09$reactome_pathway)])
union_down$mlog_qval_27 <- -log10(mets_down_list$Patient_27$q.value[match(union_down$Var1, mets_down_list$Patient_27$reactome_pathway)])
union_down$mlog_qval_44 <- -log10(mets_down_list$Patient_44$q.value[match(union_down$Var1, mets_down_list$Patient_44$reactome_pathway)])
union_down$mlog_qval_49 <- -log10(mets_down_list$Patient_49$q.value[match(union_down$Var1, mets_down_list$Patient_49$reactome_pathway)])
union_down$mlog_qval_avg <- rowMeans(union_down[,c("mlog_qval_09", "mlog_qval_27", "mlog_qval_44", "mlog_qval_49")], na.rm = T)
union_down <- union_down[order(union_down$Freq, union_down$mlog_qval_avg, decreasing = T),]

specific_deg_list<-lapply(list(up = union_up$Var1[union_up$Freq >=2], down = union_down$Var1[union_down$Freq >=2]), function(x) {
    res<-reactomeSets[as.character(x)]
    res[sapply(res, length) >=5]
})

remove_redundancy <- function(deg_list, red_cut=.5) {
    return_list <- lapply(deg_list, function(x) {
        for(i in 1:(length(x)-1)){
            if(!is.na(x[[i]])) {
                for(j in (i+1):length(x)) {
                    if(!is.na(x[[j]])) {
                        x_i <- unlist(x[[i]])
                        x_j <- unlist(x[[j]])
                        res<-length(intersect(x[[i]], x[[j]]))/length(union(x_i, x_j))
                        if(res >= red_cut) x[[j]] <- NA
                    }
                }
            }
        }
        return(x)
    })
    return_list <-lapply(return_list, function(x) names(x)[!is.na(x)])
    return(return_list)
}
plot_list<-remove_redundancy(specific_deg_list, red_cut = .3)

use_color = c(
    "Tumor_09" = "#a6cee3",
    "Tumor_27" = "#b2df8a",
    "Tumor_44" = "#fb9a99",
    "Tumor_49" = "#fdbf6f",
    "LiverMet_09" = "#1f78b4",
    "LiverMet_27" = "#33a02c",
    "LiverMet_44" = "#e31a1c",
    "LiverMet_49" = "#ff7f00"
)
plot_meta <-  pData(cur_eset)[colnames(matrixAUC_test),]
plot_meta <- plot_meta[plot_meta$Dataset %in% names(use_color),]
plot_meta$Dataset <- factor(plot_meta$Dataset, levels = names(use_color))
plot_meta <- plot_meta[order(plot_meta$Dataset),]

max_num = 9
#plotp <- c(plot_list$up[1:max_num], plot_list$down[1:max_num])
plotp <- plot_list$up[1:max_num]
value  <- matrixAUC[plotp,rownames(plot_meta)]
value <- t(scale(t(as.matrix(value))))
limits <- c(-1.5,1.5)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}


pdf(paste0(mplotPath, test_text, "upinmets_pathway_", "heatmap",".pdf"), width = 7, height=1.8)
pheatmap::pheatmap(value,
                   color = get_numeric_color("RdBu"),
                   cluster_rows = F, cluster_cols = F,
                   clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta[,c("Dataset"),drop=F], annotation_names_row = FALSE,
                   annotation_names_col = FALSE,
                   annotation_colors = list("Dataset" = use_color),
                   fontsize = 8)
dev.off()

# R/L interaction

############## Receptor - Ligand ###############
PairsLigRec <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/PairsLigRec.txt")

############## Signature of non-epi cells ###############
mets_nonepi_deg<- read_excel_allsheets(paste0(msheetPath, "met_only_Myeloid_vs_Plasma cell_vs_T cell_vs_B cell_vs_Mast cell_vs_Fibroblast_vs_Endothelial_vs_Myofibroblast_2021-04-22_de_significant.xlsx"))
non_epi_sig<-lapply(mets_nonepi_deg, function(x) {
    x$is_ligand <- x$gene_short_name %in% as.character(PairsLigRec$Ligand.ApprovedSymbol)
    x$is_receptor <- x$gene_short_name %in% as.character(PairsLigRec$Receptor.ApprovedSymbol)
    return(x)
})

# mets_sig <- read.csv(paste0(msheetPath, "mets_epi_deg_union.csv"))
# mets_sig$X <- NULL
# mets_sig <- mets_sig[mets_sig$Freq >=1, ] 

mets_09_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_09_vs_Tumor_09_2021-04-22_de_significant.xlsx"))
mets_27_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_27_vs_Tumor_27_2021-04-22_de_significant.xlsx"))
mets_44_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_44_vs_Tumor_44_2021-04-22_de_significant.xlsx"))
mets_49_deg<-read_excel_allsheets(paste0(msheetPath, "mets_epi_cleaned_LiverMet_49_vs_Tumor_49_2021-04-22_de_significant.xlsx"))

mets_epi_deg<- list(p_09 = mets_09_deg$LiverMet_09, p_27 = mets_27_deg$LiverMet_27, p_44 = mets_44_deg$LiverMet_44, p_49 = mets_49_deg$LiverMet_49)

mets_sig <- lapply(mets_epi_deg, function(x) {
    x$is_ligand <- x$gene_short_name %in% as.character(PairsLigRec$Ligand.ApprovedSymbol)
    print(sum(x$is_ligand))
    x$is_receptor <- x$gene_short_name %in% as.character(PairsLigRec$Receptor.ApprovedSymbol)
    print(sum(x$is_receptor))
    return(x)
})

edge_colors <- c(
    "Tumor Epithelial" = "black",
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

# Plot interaction using ggraph
##### Ligand on non-epi talk to Receptors on tumor ####
library(ggraph)
library(igraph)
mets_receptors_list <- lapply(mets_sig, function(x) {
    x$gene_short_name[x$is_receptor]
})
mets_receptors<-as.data.frame(table(unlist(mets_receptors_list)))
mets_receptors <- mets_receptors[mets_receptors$Freq != 0,]
mets_receptors$Patient <- sapply(mets_receptors$Var1, function(x){
    x <- as.character(x)
    res<-sapply(names(mets_receptors_list), function(p) {
        if(x %in% as.character(mets_receptors_list[[p]])) return(p) else return(NA)
    })
    paste0(res[!is.na(res)], collapse = ",")
})
saveRDS(mets_receptors, paste0(mstatePath, test_text, "mets_receptors.rds"))

mets_ligands_list <- lapply(mets_sig, function(x) {
    x$gene_short_name[x$is_ligand]
})
mets_ligands<-as.data.frame(table(unlist(mets_ligands_list)))
mets_ligands <- mets_ligands[mets_ligands$Freq != 0,]
mets_ligands$Patient <- sapply(mets_ligands$Var1, function(x){
    x <- as.character(x)
    res<-sapply(names(mets_ligands_list), function(p) {
        if(x %in% as.character(mets_ligands_list[[p]])) return(p) else return(NA)
    })
    paste0(res[!is.na(res)], collapse = ",")
})
saveRDS(mets_ligands, paste0(mstatePath, test_text, "mets_ligands.rds"))

lig_non_epi<-lapply(1:length(non_epi_sig), function(i) {
    x <- non_epi_sig[[i]]
    res <- data.frame(ligand = as.character(x$gene_short_name[x$is_ligand]))
    res$receptor <- as.character(PairsLigRec$Receptor.ApprovedSymbol[match(res$ligand, PairsLigRec$Ligand.ApprovedSymbol)])
    res$source <- names(non_epi_sig[i])
    res$lwd <- mets_receptors$Freq[match(res$receptor, mets_receptors$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(non_epi_sig[i])
    return(res)
})
lig_plot <- bind_rows(lig_non_epi)
# Circular plot
#lig_plot$ligand <- paste0(lig_plot$ligand, "_", lig_plot$source)
verts <- data.frame(name = unique(c(as.character(lig_plot$ligand), as.character(lig_plot$receptor))))
verts$name <- as.character(verts$name)
verts$gene_type <- c("TRUE" = "ligand", "FALSE" = "receptor")[as.character(verts$name %in% lig_plot$ligand)]
verts$source <- lig_plot$source[match(verts$name, lig_plot$ligand)]
verts$source <- ifelse(is.na(verts$source), "Tumor Epithelial", verts$source)

lig_graph<-graph_from_data_frame(lig_plot, vertices = verts)
pdf(paste0(mplotPath, test_text, "nonepi_ligand_tumor_receptor_circ_mplot.pdf"), width = 3, height = 3)
ggraph(lig_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(color = source), size = 1, show.legend = F) + 
    #geom_node_label(aes(label = name, color = source),  show.legend = F) + 
    scale_edge_width(range = c(0.1,.5)) + 
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()

lig_degree <- degree(lig_graph, mode = "in")
lig_sheet <- lig_plot
lig_sheet$in_degree <- lig_degree[match(lig_sheet$receptor, names(lig_degree))]
colnames(lig_sheet)[colnames(lig_sheet) == "lwd"] <- "sample_overexpress"
lig_sheet <- lig_sheet[order(lig_sheet$sample_overexpress, lig_sheet$in_degree, lig_sheet$receptor, decreasing = T),]
write.xlsx(lig_sheet, paste0(msheetPath, test_text, "nonepi_lig_tumor_receptor_wt_degree.xlsx"))

top_list <- lig_sheet
top_list$ligand <- NULL;top_list$source = NULL
top_list <- top_list[!duplicated(top_list),]   
rownames(top_list) <- top_list$receptor; top_list$receptor <- NULL

pdf(paste0(mplotPath, test_text, "R_s", ".pdf"), width = 1, height = 2.6)
pheatmap(top_list[1:15,1,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#e41a1c"))(10), show_colnames=F,number_format ="%.0f",legend = F)
pheatmap(top_list[1:15,2,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#377eb8"))(10), show_colnames=F,number_format ="%.0f",legend = F)
dev.off()




##### Receptor on non-epi talk to Ligands on mets ####
rec_non_epi<-lapply(1:length(non_epi_sig), function(i) {
    x <- non_epi_sig[[i]]
    res <- data.frame(receptor = as.character(x$gene_short_name[x$is_receptor]))
    res$ligand <- as.character(PairsLigRec$Ligand.ApprovedSymbol[match(res$receptor, PairsLigRec$Receptor.ApprovedSymbol)])
    res$lwd <- mets_ligands$Freq[match(res$ligand, mets_ligands$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(non_epi_sig[i])
    return(res)
})

rec_bind <- bind_rows(rec_non_epi)
saveRDS(rec_bind, paste0(mstatePath, test_text, "rec_bind.rds"))
# Circular plot
rec_plot <- rec_bind
verts <- data.frame(name = unique(c(as.character(rec_plot$ligand), as.character(rec_plot$receptor))))
verts$name <- as.character(verts$name)
verts$gene_type <- c("TRUE" = "ligand", "FALSE" = "receptor")[as.character(verts$name %in% rec_plot$ligand)]
verts$source <- rec_plot$source[match(verts$name, rec_plot$receptor)]
verts$source <- ifelse(is.na(verts$source), "Tumor Epithelial", verts$source)

rec_graph<-graph_from_data_frame(rec_plot, vertices = verts)
pdf(paste0(mplotPath, test_text, "nonepi_receptor_mets_ligand_circ.pdf"), width = 35, height = 35)
ggraph(rec_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(shape = gene_type, color = source), size = 2, show.legend = F) + 
    geom_node_label(aes(label = name, color = source),  show.legend = F) + 
    scale_edge_width(range = c(0.2,2)) + 
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()

pdf(paste0(mplotPath, test_text, "nonepi_receptor_mets_ligand_circ_mplot.pdf"), width = 3, height = 3)
ggraph(rec_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(shape = gene_type, color = source), size = 1, show.legend = F) + 
    #geom_node_label(aes(label = name, color = source),  show.legend = F) + 
    scale_edge_width(range = c(0.1,.5)) + 
    #scale_label_size_continuous(range = c(1,5), trans = "log10")+
    theme() +
    theme_void()
dev.off()

rec_degree <- degree(rec_graph, mode = "in")
rec_sheet <- rec_plot
rec_sheet$in_degree <- rec_degree[match(rec_sheet$ligand, names(rec_degree))]
colnames(rec_sheet)[colnames(rec_sheet) == "lwd"] <- "sample_overexpress"
rec_sheet <- rec_sheet[order(rec_sheet$sample_overexpress, rec_sheet$in_degree, rec_sheet$ligand, decreasing = T),]
write.xlsx(rec_sheet, paste0(msheetPath,test_text,"nonepi_rec_tumor_ligand_wt_degree.xlsx"))


top_list <- rec_sheet
top_list$receptor <- NULL;top_list$source = NULL
top_list <- top_list[!duplicated(top_list),]   
rownames(top_list) <- top_list$ligand; top_list$ligand <- NULL

breaksList <- seq(0,3,by=1)

pdf(paste0(mplotPath, test_text, "L_s", ".pdf"), width = .85, height = 2.6)
pheatmap(top_list[1:15,1,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#e41a1c"))(length(breaksList)), show_colnames=F,number_format ="%.0f",breaks = breaksList, legend = F)
pheatmap(top_list[1:15,2,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#377eb8"))(10), show_colnames=F,number_format ="%.0f",legend = F)
dev.off()

