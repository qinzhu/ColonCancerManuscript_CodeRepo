
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

# clist <- readRDS("../hcc_cello/clist.rds")
# eset <- readRDS("../hcc_cello/eset.rds")
clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")

# First get DEGs for all non-epithelial cell types
test_text<-"de_non_epi_primary_2022_"

#dep_vis <- clist$`All in vivo cells [no OGN, IFF, low-MT]`
dep_vis <- clist$`All in vivo cells [20210411, IFF]`
eset_vivo <- eset[, dep_vis@idx]
eset_vivo <- eset_vivo[,eset_vivo$Patient!="40"]
eset_tum <- eset_vivo[,eset_vivo$SampleType == "Tumor"]
use_ctype <- names(table(eset_tum$Cell_type))
use_ctype <- use_ctype[!use_ctype%in% c("Epithelial(Normal)", "Epithelial(Tumor)")]
eset_de<- eset_tum
eset_de$de_group <- eset_de$Cell_type
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)
bg_count <- 3000
max_count <- 500

pdeg_list<- list()
for(test_ct in use_ctype){
    print(test_ct)
    test_pair <- c(test_ct, "Background")
    group1_idx <- which(pData(eset_de)$de_group == test_ct)
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <-  which(pData(eset_de)$de_group != test_ct)
    group2_idx <- sample(group2_idx, min(length(group2_idx), bg_count))
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
    pdeg_list[[test_ct]] <- prioritized_genes[[1]]
}

saveRDS(pdeg_list, paste0(mstatePath, test_text, "pdeg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))



############## Receptor - Ligand ###############
PairsLigRec <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/PairsLigRec.txt")

############## Signature of non-epi cells ###############
non_epi_sig<-readRDS(paste0(mstatePath,"de_non_epi_primary_","de_sig_genes.rds"))
non_epi_sig<-lapply(non_epi_sig, function(x) {
    x$is_ligand <- x$gene_name %in% as.character(PairsLigRec$Ligand.ApprovedSymbol)
    x$is_receptor <- x$gene_name %in% as.character(PairsLigRec$Receptor.ApprovedSymbol)
    return(x)
})

tumor_sig <- readRDS(paste0(mstatePath,"tumor_vs_colon_de_","de_sig_genes.rds"))
tumor_sig<-lapply(tumor_sig, function(x) {
    x$is_ligand <- x$gene_name %in% as.character(PairsLigRec$Ligand.ApprovedSymbol)
    print(sum(x$is_ligand))
    x$is_receptor <- x$gene_name %in% as.character(PairsLigRec$Receptor.ApprovedSymbol)
    print(sum(x$is_receptor))
    return(x)
})









test_text <- "human_rl_2022_"

edge_colors <- c(
    "Tumor Epithelial" = "black",
    "Fibroblast" = "#cab2d6", 
    "Myofibroblast" = "#fb9a99",
    "Endothelial"= "#33a02c",
    "Myeloid" = "#6a3d9a",
    "T cell" ="#a6cee3",
    "B cell" = "#fdbf6f",
    "Plasma cell" = "#ff7f00",
    "Mast cell" = "#ae017e"
)

# Plot interaction using ggraph
##### Ligand on non-epi talk to Receptors on tumor ####
library(ggraph)
library(igraph)
tumor_receptors_list <- lapply(tumor_sig, function(x) {
    x$gene_name[x$is_receptor]
})
tumor_receptors<-as.data.frame(table(unlist(tumor_receptors_list)))
tumor_receptors <- tumor_receptors[tumor_receptors$Freq != 0,]
tumor_receptors$Patient <- sapply(tumor_receptors$Var1, function(x){
    x <- as.character(x)
    res<-sapply(names(tumor_receptors_list), function(p) {
        if(x %in% as.character(tumor_receptors_list[[p]])) return(p) else return(NA)
    })
    paste0(res[!is.na(res)], collapse = ",")
})
saveRDS(tumor_receptors, paste0(mstatePath, test_text, "tumor_receptors.rds"))

tumor_ligands_list <- lapply(tumor_sig, function(x) {
    x$gene_name[x$is_ligand]
})
tumor_ligands<-as.data.frame(table(unlist(tumor_ligands_list)))
tumor_ligands <- tumor_ligands[tumor_ligands$Freq != 0,]
tumor_ligands$Patient <- sapply(tumor_ligands$Var1, function(x){
    x <- as.character(x)
    res<-sapply(names(tumor_ligands_list), function(p) {
        if(x %in% as.character(tumor_ligands_list[[p]])) return(p) else return(NA)
    })
    paste0(res[!is.na(res)], collapse = ",")
})
saveRDS(tumor_ligands, paste0(mstatePath, test_text, "tumor_ligands.rds"))

tumor_receptors <- readRDS(paste0(mstatePath, test_text, "tumor_receptors.rds"))
tumor_ligands <- readRDS(paste0(mstatePath, test_text, "tumor_ligands.rds"))

lig_non_epi<-lapply(1:length(non_epi_sig), function(i) {
    x <- non_epi_sig[[i]]
    res <- data.frame(ligand = as.character(x$gene_name[x$is_ligand]))
    res$receptor <- as.character(PairsLigRec$Receptor.ApprovedSymbol[match(res$ligand, PairsLigRec$Ligand.ApprovedSymbol)])
    res$lwd <- tumor_receptors$Freq[match(res$receptor, tumor_receptors$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(non_epi_sig[i])
    return(res)
})

lig_bind <- bind_rows(lig_non_epi)
saveRDS(lig_bind, paste0(mstatePath, test_text, "lig_bind.rds"))
lig_plot <- lig_bind
# Circular plot
lig_plot$ligand <- paste0(lig_plot$ligand, "_", lig_plot$source)
verts <- data.frame(name = unique(c(as.character(lig_plot$ligand), as.character(lig_plot$receptor))))
verts$name <- as.character(verts$name)
verts$gene_type <- c("TRUE" = "ligand", "FALSE" = "receptor")[as.character(verts$name %in% lig_plot$ligand)]
verts$source <- lig_plot$source[match(verts$name, lig_plot$ligand)]
verts$source <- ifelse(is.na(verts$source), "Tumor Epithelial", verts$source)

lig_graph<-graph_from_data_frame(lig_plot, vertices = verts)
pdf(paste0(mplotPath, test_text, "nonepi_ligand_tumor_receptor_circ.pdf"), width = 35, height = 35)
ggraph(lig_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
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


pdf(paste0(mplotPath, test_text, "nonepi_ligand_tumor_receptor_circ_mplot.pdf"), width = 3, height = 3)
ggraph(lig_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
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

pdf(paste0(mplotPath, test_text, "R_s", ".pdf"), width = 1.18, height = 2.6)
pheatmap(top_list[1:15,1,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#e41a1c"))(10), show_colnames=F,number_format ="%.0f",)
pheatmap(top_list[1:15,2,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#377eb8"))(10), show_colnames=F,number_format ="%.0f",)
dev.off()




    

##### Receptor on non-epi talk to Ligands on tumor ####
rec_non_epi<-lapply(1:length(non_epi_sig), function(i) {
    x <- non_epi_sig[[i]]
    res <- data.frame(receptor = as.character(x$gene_name[x$is_receptor]))
    res$ligand <- as.character(PairsLigRec$Ligand.ApprovedSymbol[match(res$receptor, PairsLigRec$Receptor.ApprovedSymbol)])
    res$lwd <- tumor_ligands$Freq[match(res$ligand, tumor_ligands$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(non_epi_sig[i])
    return(res)
})

rec_bind <- bind_rows(rec_non_epi)
saveRDS(rec_bind, paste0(mstatePath, test_text, "rec_bind.rds"))
# Circular plot
rec_plot <- rec_bind
rec_plot$receptor <- paste0(rec_plot$receptor, "_", rec_plot$source)
verts <- data.frame(name = unique(c(as.character(rec_plot$ligand), as.character(rec_plot$receptor))))
verts$name <- as.character(verts$name)
verts$gene_type <- c("TRUE" = "ligand", "FALSE" = "receptor")[as.character(verts$name %in% rec_plot$ligand)]
verts$source <- rec_plot$source[match(verts$name, rec_plot$receptor)]
verts$source <- ifelse(is.na(verts$source), "Tumor Epithelial", verts$source)

rec_graph<-graph_from_data_frame(rec_plot, vertices = verts)
pdf(paste0(mplotPath, test_text, "nonepi_receptor_tumor_ligand_circ.pdf"), width = 35, height = 35)
ggraph(rec_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
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

pdf(paste0(mplotPath, test_text, "nonepi_receptor_tumor_ligand_circ_mplot.pdf"), width = 3, height = 3)
ggraph(rec_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
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

pdf(paste0(mplotPath, test_text, "L_s", ".pdf"), width = 1.2, height = 2.6)
pheatmap(top_list[1:15,1,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#e41a1c"))(10), show_colnames=F,number_format ="%.0f",)
pheatmap(top_list[1:15,2,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#377eb8"))(10), show_colnames=F,number_format ="%.0f",)
dev.off()



test_text <- "hs_mm_rl_"

human_receptors <- readRDS(paste0(mstatePath, "human_rl_2022_tumor_receptors.rds"))
human_ligands <- readRDS(paste0(mstatePath, "human_rl_2022_tumor_ligands.rds"))
mouse_receptors <- readRDS(paste0(mstatePath, "mouse_rl_tumor_receptors.rds"))
mouse_ligands <- readRDS(paste0(mstatePath, "mouse_rl_tumor_ligands.rds"))

# Compare human RL interactions with mouse
human_tumor_rec <- read.xlsx(paste0(msheetPath, "human_rl_2022_nonepi_lig_tumor_receptor_wt_degree.xlsx"))
human_tumor_lig <- read.xlsx(paste0(msheetPath, "human_rl_2022_nonepi_rec_tumor_ligand_wt_degree.xlsx"))
mouse_tumor_rec <- read.xlsx(paste0(msheetPath, "mouse_rl_nonepi_lig_tumor_receptor_wt_degree2.xlsx"))
mouse_tumor_lig <- read.xlsx(paste0(msheetPath, "mouse_rl_nonepi_rec_tumor_ligand_wt_degree2.xlsx"))

human_tumor_rec$samples <- human_receptors$Patient[match(human_tumor_rec$receptor, human_receptors$Var1)]
human_tumor_rec <- human_tumor_rec %>% dplyr::select(receptor=receptor, ligand = ligand, ligand_source=source, num_sample_overexpress = sample_overexpress, sample_overexpress = samples, in_degree = in_degree)

human_tumor_lig$samples <- human_ligands$Patient[match(human_tumor_lig$ligand, human_ligands$Var1)]
human_tumor_lig <- human_tumor_lig %>% dplyr::select(ligand = ligand, receptor=receptor, receptor_source=source, num_sample_overexpress = sample_overexpress, sample_overexpress = samples, in_degree = in_degree)

HMD_HumanPhenotype <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/HMD_HumanPhenotype.rpt", header = F)
human_tumor_rec$mouse_ortholog <- as.character(HMD_HumanPhenotype$V5)[match(as.character(human_tumor_rec$receptor), as.character(HMD_HumanPhenotype$V1))]
human_tumor_lig$mouse_ortholog <- as.character(HMD_HumanPhenotype$V5)[match(as.character(human_tumor_lig$ligand), as.character(HMD_HumanPhenotype$V1))]

human_tumor_rec$mouse_ortholog_de_in_mouse_tumor <- ifelse(human_tumor_rec$mouse_ortholog %in% mouse_receptors$Var1, T,F)
human_tumor_lig$mouse_ortholog_de_in_mouse_tumor <- ifelse(human_tumor_lig$mouse_ortholog %in% mouse_ligands$Var1, T,F)

human_tumor_rec$mouse_tumor_overexpress <- mouse_receptors$Tumor_sample[match(human_tumor_rec$mouse_ortholog, mouse_receptors$Var1)]
human_tumor_lig$mouse_tumor_overexpress <- mouse_ligands$Tumor_sample[match(human_tumor_lig$mouse_ortholog, mouse_ligands$Var1)]

write.xlsx(human_tumor_rec, paste0(msheetPath, test_text, "human_tumor_receptors.xlsx"))
write.xlsx(human_tumor_lig, paste0(msheetPath, test_text, "human_tumor_ligands.xlsx"))

mouse_tumor_rec$mouse_tumor_overexpress <- mouse_receptors$Tumor_sample[match(mouse_tumor_rec$receptor, mouse_receptors$Var1)]
mouse_tumor_lig$mouse_tumor_overexpress <- mouse_ligands$Tumor_sample[match(mouse_tumor_lig$ligand, mouse_ligands$Var1)]

write.xlsx(mouse_tumor_rec, paste0(msheetPath, test_text, "mouse_tumor_receptors.xlsx"))
write.xlsx(mouse_tumor_lig, paste0(msheetPath, test_text, "mouse_tumor_ligands.xlsx"))

# tumoroid28_de_tbl <- read_excel_allsheets(paste0(msheetPath, "Epithelial cells [wt OGN, IFF]_Tumor_28_vs_Tumoroid_28_24hr_Tumoroid_28_48hr_2020-07-20_de_all.xlsx"))
# 
# human_tumor_rec$LFC_vivo_vs_vitro <- tumoroid28_de_tbl$Tumor_28$log2fc[match(human_tumor_rec$receptor, tumoroid28_de_tbl$Tumor_28$gene_short_name)]
# human_tumor_rec$Pval_vivo_vs_vitro <- tumoroid28_de_tbl$Tumor_28$p_adj[match(human_tumor_rec$receptor, tumoroid28_de_tbl$Tumor_28$gene_short_name)]
# 
# human_tumor_lig$LFC_vivo_vs_vitro <- tumoroid28_de_tbl$Tumor_28$log2fc[match(human_tumor_lig$ligand, tumoroid28_de_tbl$Tumor_28$gene_short_name)]
# human_tumor_lig$Pval_vivo_vs_vitro <- tumoroid28_de_tbl$Tumor_28$p_adj[match(human_tumor_lig$ligand, tumoroid28_de_tbl$Tumor_28$gene_short_name)]
# 
# write.xlsx(human_tumor_rec[order(human_tumor_rec$mouse_ortholog_de_in_mouse_tumor, human_tumor_rec$num_sample_overexpress, human_tumor_rec$LFC_vivo_vs_vitro, human_tumor_rec$in_degree, human_tumor_rec$receptor, decreasing = T),], paste0(msheetPath, test_text, "tumor_receptors_stringent_0722.xlsx"))
# write.xlsx(human_tumor_lig[order(human_tumor_lig$mouse_ortholog_de_in_mouse_tumor, human_tumor_lig$num_sample_overexpress, human_tumor_lig$LFC_vivo_vs_vitro, human_tumor_lig$in_degree, human_tumor_lig$ligand, decreasing = T),], paste0(msheetPath, test_text, "tumor_ligands_stringent_0722.xlsx"))
# 
# 
# write.xlsx(list(receptors =human_tumor_rec[order(human_tumor_rec$mouse_ortholog_de_in_mouse_tumor, human_tumor_rec$num_sample_overexpress, human_tumor_rec$LFC_vivo_vs_vitro, human_tumor_rec$in_degree, human_tumor_rec$receptor, decreasing = T),], ligands = human_tumor_lig[order(human_tumor_lig$mouse_ortholog_de_in_mouse_tumor, human_tumor_lig$num_sample_overexpress, human_tumor_lig$LFC_vivo_vs_vitro, human_tumor_lig$in_degree, human_tumor_lig$ligand, decreasing = T),]), paste0(msheetPath, test_text, "tumor_rl_stringent_0722.xlsx"))



# Intersection wo caring whether the paired ligand/receptor can be found
human_rec_all <- human_receptors
colnames(human_rec_all) <- c("receptor", "num_sample_overexpress", "sample_overexpress")
human_rec_all$mouse_ortholog <- as.character(HMD_HumanPhenotype$V5)[match(as.character(human_rec_all$receptor), as.character(HMD_HumanPhenotype$V1))]
human_rec_all$mouse_ortholog_de_in_mouse_tumor <- ifelse(human_rec_all$mouse_ortholog %in% mouse_receptors$Var1, T,F)
human_rec_all$num_mouse_tumor_overexpress <- mouse_receptors$Freq[match(human_rec_all$mouse_ortholog, mouse_receptors$Var1)]
human_rec_all$mouse_tumor_overexpress <- mouse_receptors$Tumor_sample[match(human_rec_all$mouse_ortholog, mouse_receptors$Var1)]

write.xlsx(human_rec_all[order(human_rec_all$mouse_ortholog_de_in_mouse_tumor, human_rec_all$num_sample_overexpress, human_rec_all$num_mouse_tumor_overexpress, human_rec_all$receptor, decreasing = T),], paste0(msheetPath, test_text, "tumor_receptors_loose2.xlsx"))


human_lig_all <- human_ligands
colnames(human_lig_all) <- c("ligand", "num_sample_overexpress", "sample_overexpress")
human_lig_all$mouse_ortholog <- as.character(HMD_HumanPhenotype$V5)[match(as.character(human_lig_all$ligand), as.character(HMD_HumanPhenotype$V1))]
human_lig_all$mouse_ortholog_de_in_mouse_tumor <- ifelse(human_lig_all$mouse_ortholog %in% mouse_ligands$Var1, T,F)
human_lig_all$num_mouse_tumor_overexpress <- mouse_ligands$Freq[match(human_lig_all$mouse_ortholog, mouse_ligands$Var1)]
human_lig_all$mouse_tumor_overexpress <- mouse_ligands$Tumor_sample[match(human_lig_all$mouse_ortholog, mouse_ligands$Var1)]

write.xlsx(human_lig_all[order(human_lig_all$mouse_ortholog_de_in_mouse_tumor, human_lig_all$num_sample_overexpress, human_lig_all$num_mouse_tumor_overexpress, human_lig_all$ligand, decreasing = T),], paste0(msheetPath, test_text, "tumor_ligands_loose2.xlsx"))


### Intersect venn diagram

# Receptor
library(VennDiagram)
human_rec_ortholog_stringent <- unique(human_tumor_rec$mouse_ortholog)
human_rec_ortholog_stringent <- human_rec_ortholog_stringent[!is.na(human_rec_ortholog_stringent)]
human_rec_ortholog_loose <- unique(human_rec_all$mouse_ortholog)
human_rec_ortholog_loose <- human_rec_ortholog_loose[!is.na(human_rec_ortholog_loose)]

rec_list <- list(human_stringent = as.character(human_rec_ortholog_stringent), human_loose = as.character(human_rec_ortholog_loose), mouse_stringent = as.character(unique(mouse_tumor_rec$receptor)), mouse_loose = as.character(unique(mouse_receptors$Var1)))

use_color <- c("human_stringent" = "#377eb8", "human_loose" = "#4daf4a", "mouse_stringent" = "#a65628", "mouse_loose" = "#ff7f00")
use_list <- c("human_stringent", "mouse_stringent")
pdf(paste0(mplotPath, test_text,"venn_receptor_stringent2.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(rec_list[use_list], NULL, cex = 2, fill=use_color[use_list], cat.fontface=4, main = "Receptor intersection")
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()

use_list <- c("human_loose", "mouse_loose")
pdf(paste0(mplotPath, test_text,"venn_receptor_loose2.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(rec_list[use_list], NULL, cex = 2, fill=use_color[use_list], cat.fontface=4, main = "Receptor intersection")
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()


# Ligand

human_lig_ortholog_stringent <- unique(human_tumor_lig$mouse_ortholog)
human_lig_ortholog_stringent <- human_lig_ortholog_stringent[!is.na(human_lig_ortholog_stringent)]
human_lig_ortholog_loose <- unique(human_lig_all$mouse_ortholog)
human_lig_ortholog_loose <- human_lig_ortholog_loose[!is.na(human_lig_ortholog_loose)]

lig_list <- list(human_stringent = as.character(human_lig_ortholog_stringent), human_loose = as.character(human_lig_ortholog_loose), mouse_stringent = as.character(unique(mouse_tumor_lig$ligand)), mouse_loose = as.character(unique(mouse_ligands$Var1)))

use_color <- c("human_stringent" = "#377eb8", "human_loose" = "#4daf4a", "mouse_stringent" = "#a65628", "mouse_loose" = "#ff7f00")
use_list <- c("human_stringent", "mouse_stringent")
pdf(paste0(mplotPath, test_text,"venn_ligand_stringent2.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(lig_list[use_list], NULL, cex = 2, fill=use_color[use_list], cat.fontface=4, main = "Ligand intersection")
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()

use_list <- c("human_loose", "mouse_loose")
pdf(paste0(mplotPath, test_text,"venn_ligand_loose2.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(lig_list[use_list], NULL, cex = 2, fill=use_color[use_list], cat.fontface=4, main = "Ligand intersection")
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()


