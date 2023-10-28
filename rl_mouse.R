




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
source(paste0(scriptPath, "compute_go.R"))

clist <- readRDS("../mouse_colon_cello/clist.rds")
eset <- readRDS("../mouse_colon_cello/eset.rds")


# First get DEGs for all non-epithelial cell types
test_text<-"mouse_de_non_epi_" 
dep_vis <- clist$`Mouse all data [IFF, low mito]`
eset_vivo <- eset[, dep_vis@idx]
eset_vivo <- eset_vivo[, !eset_vivo$Dataset %in% c("ApcFlox", "BcatFlox")]
use_ctype <- names(table(eset_vivo$Cell_type_refined))
use_ctype <- c("B & Plasma cell", "Fibroblast", "Endothelial", "Myeloid", "T cells")
eset_de<- eset_vivo[,!eset_vivo$Cell_type_refined %in% c("Epithelial cells (tumor)", "Epithelial cells [EC,TA,SC]", "Goblet cells", "EEC cells")]
eset_de$de_group <- eset_de$Cell_type_refined
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)
bg_count <- 2000
max_count <- 300

set.seed(2020)
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
names(deg_list) <- make.names(names(deg_list))
saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))






set.seed(2020)
test_text <- "mouse_tumor_vs_background_de_"
dep_vis <- clist$`Mouse all data [IFF, low mito]`
eset_vivo <- eset[, dep_vis@idx]
eset_vivo <- eset_vivo[, !eset_vivo$Dataset %in% c("ApcFlox", "BcatFlox")]
bg_count <- 2000
bg_idx <- sample(which(eset_vivo$Cell_type != "Epithelial cells (tumor)"), bg_count, replace = F)
tumor_idx <- which(eset_vivo$Cell_type == "Epithelial cells (tumor)")

max_count <- 500
eset_de <- eset_vivo[,c(tumor_idx, bg_idx)]

eset_de$de_group <- ifelse(eset_de$Cell_type == "Epithelial cells (tumor)", as.character(eset_de$Dataset), "Background")
eset_de <- eset_de[,!grepl("Colon", eset_de$de_group)]

id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

tumor_names <- grep("Tumor",names(table(eset_de$Dataset)), value = T)

pdeg_list<- list()
for(tumor in tumor_names){
    test_pair <- c(tumor, "Background")
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <-  which(pData(eset_de)$de_group == test_pair[2])
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
    pdeg_list[[tumor]] <- prioritized_genes[[1]]
}
saveRDS(pdeg_list, paste0(mstatePath, test_text, "pdeg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})

saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))





############## Receptor - Ligand ###############
PairsLigRec <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/PairsLigRec.txt")
HMD_HumanPhenotype <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/HMD_HumanPhenotype.rpt", header = F)
View(PairsLigRec)
PairsLigRec$Receptor.mouse <- as.character(HMD_HumanPhenotype$V5)[match(as.character(PairsLigRec$Receptor.ApprovedSymbol), as.character(HMD_HumanPhenotype$V1))]
PairsLigRec$Ligand.mouse <-as.character(HMD_HumanPhenotype$V5)[match(as.character(PairsLigRec$Ligand.ApprovedSymbol), as.character(HMD_HumanPhenotype$V1))]

############## Signature of non-epi cells ###############

non_epi_sig<-readRDS(paste0(mstatePath,"mouse_de_non_epi_","de_sig_genes.rds"))
non_epi_sig<-lapply(non_epi_sig, function(x) {
    x$is_ligand <- as.character(x$gene_name) %in% as.character(PairsLigRec$Ligand.mouse)
    x$is_receptor <- as.character(x$gene_name) %in% as.character(PairsLigRec$Receptor.mouse)
    return(x)
})

tumor_sig <- readRDS(paste0(mstatePath,"mouse_tumor_vs_background_de_","de_sig_genes.rds"))
tumor_sig<-lapply(tumor_sig, function(x) {
    x$is_ligand <- as.character(x$gene_name) %in% as.character(PairsLigRec$Ligand.mouse)
    print(sum(x$is_ligand))
    x$is_receptor <- as.character(x$gene_name) %in% as.character(PairsLigRec$Receptor.mouse)
    print(sum(x$is_receptor))
    return(x)
})



test_text <- "mouse_rl_"
library(ggraph)
library(igraph)
tumor_receptors_list <- lapply(tumor_sig, function(x) {
    x$gene_name[x$is_receptor]
})
tumor_receptors<-as.data.frame(table(unlist(tumor_receptors_list)))
tumor_receptors <- tumor_receptors[tumor_receptors$Freq != 0,]
tumor_receptors$Tumor_sample <- sapply(tumor_receptors$Var1, function(x){
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
tumor_ligands$Tumor_sample <- sapply(tumor_ligands$Var1, function(x){
    x <- as.character(x)
    res<-sapply(names(tumor_ligands_list), function(p) {
        if(x %in% as.character(tumor_ligands_list[[p]])) return(p) else return(NA)
    })
    paste0(res[!is.na(res)], collapse = ",")
})
saveRDS(tumor_ligands, paste0(mstatePath, test_text, "tumor_ligands.rds"))



lig_non_epi<-lapply(1:length(non_epi_sig), function(i) {
    x <- non_epi_sig[[i]]
    res <- data.frame(ligand = as.character(x$gene_name[x$is_ligand]))
    res$receptor <- as.character(PairsLigRec$Receptor.mouse[match(res$ligand, PairsLigRec$Ligand.mouse)])
    res$lwd <- tumor_receptors$Freq[match(res$receptor, tumor_receptors$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(non_epi_sig[i])
    return(res)
})

lig_bind <- bind_rows(lig_non_epi)

edge_colors <- c(
    "Tumor Epithelial" = "black",
    "Fibroblast" = "#cab2d6", 
    "Endothelial"= "#33a02c",
    "Myeloid" = "#6a3d9a",
    "T.cells" ="#a6cee3",
    "B...Plasma.cell" = "#ff7f00"
)



lig_plot <- lig_bind

# Circular plot
lig_plot$ligand <- paste0(lig_plot$ligand, "_", lig_plot$source)
verts <- data.frame(name = unique(c(as.character(lig_plot$ligand), as.character(lig_plot$receptor))))
verts$name <- as.character(verts$name)
verts$gene_type <- c("TRUE" = "ligand", "FALSE" = "receptor")[as.character(verts$name %in% lig_plot$ligand)]
verts$source <- lig_plot$source[match(verts$name, lig_plot$ligand)]
verts$source <- ifelse(is.na(verts$source), "Tumor Epithelial", verts$source)

lig_graph<-graph_from_data_frame(lig_plot, vertices = verts)
pdf(paste0(mplotPath, test_text,"nonepi_ligand_tumor_receptor_circ.pdf"), width = 35, height = 35)
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

pdf(paste0(mplotPath, test_text,"nonepi_ligand_tumor_receptor_circ_mplot.pdf"), width = 1.6, height = 1.6)
ggraph(lig_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(color = source), size = .3, show.legend = F) + 
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
write.xlsx(lig_sheet[order(lig_sheet$sample_overexpress, lig_sheet$in_degree, lig_sheet$receptor, decreasing = T),], paste0(msheetPath, test_text, "nonepi_lig_tumor_receptor_wt_degree2.xlsx"))





##### Ligand on non-epi talk to Receptors on tumor ####
tumor_ligands <- lapply(tumor_sig, function(x) {
    x$gene_name[x$is_ligand]
})
tumor_ligands<-as.data.frame(table(unlist(tumor_ligands)))
tumor_ligands <- tumor_ligands[tumor_ligands$Freq != 0,]

rec_non_epi<-lapply(1:length(non_epi_sig), function(i) {
    x <- non_epi_sig[[i]]
    res <- data.frame(receptor = as.character(x$gene_name[x$is_receptor]))
    res$ligand <- as.character(PairsLigRec$Ligand.mouse[match(res$receptor, PairsLigRec$Receptor.mouse)])
    res$lwd <- tumor_ligands$Freq[match(res$ligand, tumor_ligands$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(non_epi_sig[i])
    return(res)
})

rec_bind <- bind_rows(rec_non_epi)

# Circular plot
rec_plot <- rec_bind
rec_plot$receptor <- paste0(rec_plot$receptor, "_", rec_plot$source)
verts <- data.frame(name = unique(c(as.character(rec_plot$ligand), as.character(rec_plot$receptor))))
verts$name <- as.character(verts$name)
verts$gene_type <- c("TRUE" = "ligand", "FALSE" = "receptor")[as.character(verts$name %in% rec_plot$ligand)]
verts$source <- rec_plot$source[match(verts$name, rec_plot$receptor)]
verts$source <- ifelse(is.na(verts$source), "Tumor Epithelial", verts$source)

rec_graph<-graph_from_data_frame(rec_plot, vertices = verts)
pdf(paste0(mplotPath,  test_text,"nonepi_receptor_tumor_ligand_circ.pdf"), width = 35, height = 35)
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

pdf(paste0(mplotPath,  test_text,"nonepi_receptor_tumor_ligand_circ_mplot.pdf"), width = 1.6, height = 1.6)
ggraph(rec_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(color = source), size = .3, show.legend = F) + 
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
write.xlsx(rec_sheet[order(rec_sheet$sample_overexpress, rec_sheet$in_degree, rec_sheet$ligand, decreasing = T),], paste0(msheetPath, test_text,"nonepi_rec_tumor_ligand_wt_degree2.xlsx"))











