

set.seed(2020)
mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"
resultPath <- paste0(mstatePath, "de_vivo_vs_vitro/");dir.create(resultPath)

savePath <- "../hcc_final_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)
scriptPath <- paste0("../hcc_cello/scripts/")

clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")

dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

eset_epi$sample_pair <- NA
eset_epi$sample_pair[eset_epi$Dataset %in% c("Tumoroid_08", "Tumoroid_24", "Tumoroid_28")] = "Tumoroid"
eset_epi$sample_pair[eset_epi$Dataset %in% c("Tumor_08", "Tumor_24", "Tumor_28")] = "Tumor"
eset_epi$sample_pair[eset_epi$SampleType == "Colon"] = "Colon"
eset_epi$sample_pair[eset_epi$SampleType == "Organoid"] = "Organoid"

####### Compare tumoroid vs tumor ##########
test_col = "sample_pair"
test_text <- "de_tumoroid_vs_tumor_"
test_pair <- c("Tumoroid", "Tumor")
# test_text <- "de_organoid_vs_colon_"
# test_pair <- c("Organoid", "Colon")
eset_de <- eset_epi[, eset_epi$SampleType %in% test_pair]
# Derive and plot DEGs
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
max_count <- 3000
feature_data <- fData(eset_de)

eset_de$de_group <- eset_de[[test_col]]
group1_idx <- sample(which(pData(eset_de)$de_group == test_pair[1]), max_count)
group2_idx <-  sample(which(pData(eset_de)$de_group == test_pair[2]), max_count)
test_idx <- c(group1_idx, group2_idx)
test_clus <- c(rep(test_pair[1], length(group1_idx)), rep(test_pair[2], length(group2_idx)))
cur_cds <- eset_de[, test_idx]
cur_cds$test_cluster <- factor(test_clus, levels = test_pair)
#gene_idx <- Matrix::rowSums(exprs(cur_cds)) > 0
#cur_cds <- cur_cds[gene_idx,]
feature_data <- fData(cur_cds)
prioritized_genes <- runsSeq(dat=as.matrix(exprs(cur_cds)), group=cur_cds$test_cluster, fdata = feature_data, order_by="pvalue", p_cutoff= .05, min_mean = 0, min_log2fc = 0, id_col = id_col, name_col = name_col)

prioritized_genes <- lapply(prioritized_genes, function(x) {
    x$proportion1 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_pair[1]]) > 0)
    x$proportion2 <- Matrix::rowMeans(exprs(cur_cds[match(x$gene_id, rownames(cur_cds)), cur_cds$test_cluster == test_pair[2]]) > 0)
    return(x)
})

saveRDS(prioritized_genes, paste0(resultPath,test_text,"prioritized_genes.rds"))

prioritized_genes <- readRDS(paste0(resultPath,test_text,"prioritized_genes.rds"))
non_epi_sig<-readRDS(paste0(mstatePath,"de_non_epi_primary_","de_sig_genes.rds"))
#non_epi_sig_genes <- unique(unlist(lapply(non_epi_sig, function(x) x$gene_id)))
prioritized_genes<- lapply(prioritized_genes, function(x) {
    x[!x$gene_id %in% non_epi_sig$`Plasma cell`$gene_id, ]
})

# Volcano plot
library(EnhancedVolcano)
pdf(paste0(mplotPath, test_text, "deg_volcano.pdf"), width = 5, height=7)
EnhancedVolcano(prioritized_genes[[1]], 
                lab = as.character(prioritized_genes[[1]]$gene_name), x = "log2fc", y = "p_adj",
                labSize = 5,
                selectLab = as.character(prioritized_genes[[1]]$gene_name[which(prioritized_genes[[1]]$p_adj<1e-50)]),
                col=c('black', 'black', 'black', 'red3'),
                pCutoff = 1e-2, FCcutoff = 1, 
                title = "sSeq results",
                subtitle = paste0(test_pair, collapse = "_vs_"),
                raster = T)
dev.off()

                
de_list <- lapply(prioritized_genes, function(x) {
    x %>% dplyr::filter(p_adj <= .01 & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(de_list, paste0(resultPath,test_text,"de_sig_genes.rds"))
write.xlsx(de_list, paste0(resultPath,test_text,"de_sig_genes.xlsx"))



plot_meta <- pData(eset_de)
plot_meta <- plot_meta[, c(test_col, "Dataset"), drop=F]
ctype_levels = test_pair
plot_meta[[test_col]] <- factor(plot_meta[[test_col]], levels = ctype_levels)
max_cell = 500
sp_idx <- lapply(ctype_levels, function(x) {sample(which(plot_meta[[test_col]] == x), min(max_cell, sum(plot_meta[[test_col]] == x, na.rm = T)))})
plot_meta <- plot_meta[unlist(sp_idx),, drop=F]
plot_meta <- plot_meta[order(plot_meta[[test_col]], plot_meta$Dataset),]
non_na_cells_ordered <- rownames(plot_meta)

deg_sheets <- de_list
max_gene <- 20
plotg <- unique(unlist(lapply(de_list, function(x) {
    x = x %>% filter(common_mean >=.1)
    use_g <- as.character(x$gene_id)
    use_g[1:max_gene]
})))

value <- eset_de@assayData$norm_exprs[plotg,non_na_cells_ordered]
rownames(value) <- fData(eset_de)$gene_short_name[match(rownames(value), rownames(eset_de))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

pdf(paste0(resultPath, test_text, "deg_", max_gene,".pdf"), width = 6, height=4)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = F, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta[,test_col,drop=F], annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 8)
dev.off()


#test_text <- "de_tumoroid_vs_tumor_vivocleaned_"
gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
gene_background <- rownames(eset_de)
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)
go_res <- lapply(de_list, function(x) {
    include_g.df <- bitr(rownames(x), fromType = gene_id_type,
                         toType = c("SYMBOL", "ENTREZID"),
                         OrgDb = orgdb)
    include_g_go<- enrichGO(gene        = include_g.df$ENTREZID,
                            universe      = bg.df$ENTREZID,
                            OrgDb         = orgdb,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)
    return(include_g_go)
})
saveRDS(go_res, paste0(mstatePath, test_text, "go_res.rds"))
go_res <- readRDS(paste0(mstatePath, test_text, "go_res.rds"))
go_res_tbl <- lapply(go_res, function(x) {
    x@result
})
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, "go_res.xlsx"))

# GO plot
glist <- list()
for(i in 1:length(go_res)) {
    include_go_simplified <- clusterProfiler::simplify(go_res[[i]], cutoff = .7)
    plot_df<- include_go_simplified@result[1:min(nrow(include_go_simplified@result),8),]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    glist[[i]] <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
        geom_bar(stat = "identity", aes(fill = -log10(qvalue))) + 
        #scale_x_discrete(position = 'top') + 
        coord_flip() +
        xlab("") + 
        theme(
            text=element_text(family = "Helvetica", size=8,color="black"),
            axis.text = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
        ggtitle(paste0("Up-regulated in ", names(go_res)[i])) + 
        monocle:::monocle_theme_opts() 
}


pdf(paste0(mplotPath, test_text, "go_bar.pdf"), width = 5, height = 2.8)
do.call(grid.arrange, c(glist, ncol=1))
dev.off()
pdf(paste0(mplotPath, test_text, "go_bar_1.pdf"), width = 4, height = 1.7)
glist[[1]]
dev.off()
pdf(paste0(mplotPath, test_text, "go_bar_2.pdf"), width = 4, height = 1.7)
glist[[2]]
dev.off()






