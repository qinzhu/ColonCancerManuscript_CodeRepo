
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

# First get DEGs for all myeloid substates
de_list <- read_excel_allsheets(paste0("~/Dropbox/ColonManuscript/sheets/", "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))


############## Receptor - Ligand ###############
PairsLigRec <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/PairsLigRec.txt")

############## Signature of non-epi cells ###############
myeloid_sig<-lapply(de_list, function(x) {
    names(x)[names(x) == "gene_short_name"] = "gene_name"
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



test_text <- "human_rl_myeloidonly_"

edge_colors <- c("IL1B" = "#33a02c", "SPP1" = "#e31a1c", "C1QC" = "#ff7f00", "Tumor Epithelial" = "black")
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

lig_myeloid<-lapply(1:length(myeloid_sig), function(i) {
    x <- myeloid_sig[[i]]
    res <- data.frame(ligand = as.character(x$gene_name[x$is_ligand]))
    res$receptor <- as.character(PairsLigRec$Receptor.ApprovedSymbol[match(res$ligand, PairsLigRec$Ligand.ApprovedSymbol)])
    res$lwd <- tumor_receptors$Freq[match(res$receptor, tumor_receptors$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(myeloid_sig[i])
    return(res)
})

lig_bind <- bind_rows(lig_myeloid)
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
pdf(paste0(mplotPath, test_text, "myeloid_ligand_tumor_receptor_circ.pdf"), width = 35, height = 35)
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

pdf(paste0(mplotPath, test_text, "myeloid_ligand_tumor_receptor_circ_mplot.pdf"), width = 3, height = 3)
ggraph(lig_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(color = source), size = 1, show.legend = F) + 
    geom_node_text(aes(label = name, color = source), size = 3, show.legend = F) + 
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
lig_sheet$ligand <- sapply(strsplit(lig_sheet$ligand, "_"), function(x)x[1])
write.xlsx(lig_sheet, paste0(msheetPath, test_text, "myeloid_lig_tumor_receptor_wt_degree.xlsx"))


top_list <- lig_sheet
top_list$ligand <- NULL;top_list$source = NULL
top_list <- top_list[!duplicated(top_list),]   
rownames(top_list) <- top_list$receptor; top_list$receptor <- NULL

pdf(paste0(mplotPath, test_text, "R_s", ".pdf"), width = 1.15, height = 2.6)
pheatmap(top_list[1:15,1,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#e41a1c"))(10), show_colnames=F,number_format ="%.0f",)
pheatmap(top_list[1:15,2,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#377eb8"))(10), show_colnames=F,number_format ="%.0f",)
dev.off()




    

##### Receptor on non-epi talk to Ligands on tumor ####
rec_myeloid<-lapply(1:length(myeloid_sig), function(i) {
    x <- myeloid_sig[[i]]
    res <- data.frame(receptor = as.character(x$gene_name[x$is_receptor]))
    res$ligand <- as.character(PairsLigRec$Ligand.ApprovedSymbol[match(res$receptor, PairsLigRec$Receptor.ApprovedSymbol)])
    res$lwd <- tumor_ligands$Freq[match(res$ligand, tumor_ligands$Var1)]
    res <- res[complete.cases(res) & res$lwd != 0, ]
    res$source <- names(myeloid_sig[i])
    return(res)
})

rec_bind <- bind_rows(rec_myeloid)
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
pdf(paste0(mplotPath, test_text, "myeloid_receptor_tumor_ligand_circ.pdf"), width = 35, height = 35)
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

pdf(paste0(mplotPath, test_text, "myeloid_receptor_tumor_ligand_circ_mplot.pdf"), width = 3, height = 3)
ggraph(rec_graph,  layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = factor(source), edge_width = as.numeric(lwd)), show.legend = F) + 
    scale_edge_color_manual(values = edge_colors) + 
    scale_color_manual(values = edge_colors) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')+
    geom_node_point(aes(color = source), size = 1, show.legend = F) + 
    geom_node_text(aes(label = name, color = source), size = 3, show.legend = F) + 
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
rec_sheet$receptor <- sapply(strsplit(rec_sheet$receptor, "_"), function(x)x[1])
write.xlsx(rec_sheet, paste0(msheetPath,test_text,"myeloid_rec_tumor_ligand_wt_degree.xlsx"))


top_list <- rec_sheet
top_list$receptor <- NULL;top_list$source = NULL
top_list <- top_list[!duplicated(top_list),]   
rownames(top_list) <- top_list$ligand; top_list$ligand <- NULL

pdf(paste0(mplotPath, test_text, "L_s", ".pdf"), width = 1.15, height = 2.6)
pheatmap(top_list[1:15,1,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#e41a1c"))(10), show_colnames=F,number_format ="%.0f",)
pheatmap(top_list[1:15,2,drop=F], cluster_rows = F, cluster_cols = F, display_numbers =T, fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#377eb8"))(10), show_colnames=F,number_format ="%.0f",)
dev.off()


