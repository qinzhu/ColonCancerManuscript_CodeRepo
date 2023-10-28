


# Compare Tumoroid vs Organoid (in vitro)
# Compare Tumor vs colon (in vivo)
set.seed(2020)
mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"
resultPath <- paste0(mstatePath, "de_2way/");dir.create(resultPath)

savePath <- "../hcc_final_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)
scriptPath <- paste0("../hcc_cello/", "scripts/")

clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")

test_text <- "de2way_tumor_vs_colon_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]


####### Compare tumor vs normal colon ##########
test_col = "SampleType"
test_pair <- c("Tumor", "Colon")
eset_epi_vivo <- eset_epi[, eset_epi$SampleType %in% test_pair]
# Derive and plot DEGs
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
max_count <- 3000
eset_de <- eset_epi_vivo
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
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 8)
dev.off()


####### Compare tumoroid vs organoid ##########
test_text <- "de2way_tumoroid_vs_organoid_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]
test_col = "SampleType"
test_pair <- c("Tumoroid", "Organoid")
eset_epi_vitro <- eset_epi[, eset_epi$SampleType %in% test_pair]
# Derive and plot DEGs
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
max_count <- 3000
eset_de <- eset_epi_vitro
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
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 8)
dev.off()



# Compare in vivo wt in vitro
test_text <- "vivo_vs_vitro_"
vivo_de_res <- readRDS(paste0(resultPath, "de2way_tumor_vs_colon_prioritized_genes.rds"))
vitro_de_res <- readRDS(paste0(resultPath, "de2way_tumoroid_vs_organoid_prioritized_genes.rds"))

#shared_g <- intersect(rownames(vivo_de_res$Tumor), rownames(vitro_de_res$Tumoroid))
expressed_g <- rownames(eset_epi)[which(rowSums(exprs(eset_epi) > 0) > ncol(eset_epi) * .001)]

combine_de_res<- data.frame(vivo_lfc = vivo_de_res$Tumor[expressed_g, ]$log2fc, vivo_mlgP = -log10(vivo_de_res$Tumor[expressed_g, ]$p_adj), vitro_lfc = vitro_de_res$Tumoroid[expressed_g, ]$log2fc, vitro_mlgP = -log10(vitro_de_res$Tumoroid[expressed_g, ]$p_adj))
rownames(combine_de_res) <- expressed_g
combine_de_res$gene_name <-fData(eset_epi)$gene_short_name[match(expressed_g, rownames(eset_epi))]

max_g <- 5
# euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

# coord_df <- combine_de_res[,c("vivo_lfc", "vitro_lfc")] 
# lm_res <- lm(vitro_lfc~vivo_lfc+0, coord_df)
# plot(lm_res)
text_f <- unique(c(
    rownames(combine_de_res[order(combine_de_res$vivo_lfc + combine_de_res$vitro_lfc),][1:max_g,]),
    rownames(combine_de_res[order(combine_de_res$vivo_lfc + combine_de_res$vitro_lfc, decreasing = T),][1:max_g,])
    #rownames(combine_de_res[order(combine_de_res$vivo_lfc - combine_de_res$vitro_lfc, decreasing = T),][1:max_g,]),
    #rownames(combine_de_res[order(-combine_de_res$vivo_lfc + combine_de_res$vitro_lfc, decreasing = T),][1:max_g,])
))

# combine_de_res$mlgP_diff <- combine_de_res$vivo_mlgP - combine_de_res$vitro_mlgP
# diff_lim <- c(-10, 10)
# combine_de_res$mlgP_diff[combine_de_res$mlgP_diff < diff_lim[1]] = diff_lim[1]
# combine_de_res$mlgP_diff[combine_de_res$mlgP_diff > diff_lim[2]] = diff_lim[2]
# 
# combine_de_res$lfc_diff <- combine_de_res$vivo_lfc - combine_de_res$vitro_lfc
# diff_lim <- c(-10, 10)
# combine_de_res$lfc_diff[combine_de_res$lfc_diff < diff_lim[1]] = diff_lim[1]
# combine_de_res$lfc_diff[combine_de_res$lfc_diff > diff_lim[2]] = diff_lim[2]
# combine_de_res$avg_mlgP <- (combine_de_res$vivo_mlgP + combine_de_res$vitro_mlgP) / 2
# combine_de_res$avg_mlgP[combine_de_res$avg_mlgP > 10] = 10


library(ggrepel)
library(ggrastr)
# Remove genes that are non-significant in either comparison
#de_res_plot <- combine_de_res[combine_de_res$vivo_mlgP >= 2 | combine_de_res$vitro_mlgP >= 2, ]
de_res_plot <- combine_de_res
df <- de_res_plot[,c("vivo_lfc", "vitro_lfc")]
colnames(df) <- c("x", "y")
# m <- lm(y ~ x, df)
# Also highlight receptors and ligands
PairsLigRec <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/PairsLigRec.txt")
de_res_plot$is_receptor_ligand <- de_res_plot$gene_name %in% c(PairsLigRec$Receptor.ApprovedSymbol, PairsLigRec$Ligand.ApprovedSymbol)

g1 <- ggplot(de_res_plot, aes(vivo_lfc, vitro_lfc)) +
    geom_point_rast(size = .5, stroke = 0, alpha = .3) + 
    geom_text(data = combine_de_res[text_f,], aes(label = gene_name), size = 2) +
    #geom_text(x = -8, y = 10, label = lm_eqn(df), parse = TRUE) + 
    geom_abline(intercept = 0, slope = 1)+
    stat_smooth(method = "lm", formula = "y ~ x+0", se = F, color = "blue", size = .75) + 
    stat_smooth(data = de_res_plot[de_res_plot$is_receptor_ligand, ], method = "lm", formula = "y ~ x+0", se = F, color = "red", size = .75) + 
    stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon", alpha = .4) + 
    scale_fill_viridis_c() + 
    xlab("Tumor vs Colon LFC") + 
    ylab("Tumoroid vs Organoid LFC")+
    theme_bw()+
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=8)) 
ggsave(paste0(mplotPath, test_text, "lfc_vivo_vitro_220607.pdf"), g1, width = 3.5, height=2.5, units = "in", device = "pdf")
saveRDS(de_res_plot, paste0(mstatePath, test_text, "de_res_plot.rds"))

# Fit results separately
res1 <- lm(vitro_lfc~vivo_lfc + 0, data = de_res_plot)
res2 <- lm(vitro_lfc~vivo_lfc + 0, data = de_res_plot[de_res_plot$is_receptor_ligand, ])
# test if significantly different https://stats.stackexchange.com/questions/93540/testing-equality-of-coefficients-from-two-different-regressions
z_score = (summary(res1)$coefficients[1] - summary(res2)$coefficients[1]) / sqrt(summary(res1)$coefficients[2]^2 + summary(res2)$coefficients[2]^2)
pnorm(z_score, mean = 0, sd = 1, lower.tail = FALSE)



# number of genes in each quadrant
sum(de_res_plot$vivo_lfc >0 & de_res_plot$vitro_lfc > 0)
sum(de_res_plot$vivo_lfc >0 & de_res_plot$vitro_lfc <= 0)
sum(de_res_plot$vivo_lfc <=0 & de_res_plot$vitro_lfc <= 0)
sum(de_res_plot$vivo_lfc <=0 & de_res_plot$vitro_lfc > 0)

# enrichment test
q1 <- de_res_plot[de_res_plot$vivo_lfc >0 & de_res_plot$vitro_lfc > 0,]
q2 <- de_res_plot[de_res_plot$vivo_lfc >0 & de_res_plot$vitro_lfc <= 0,]
q3 <- de_res_plot[de_res_plot$vivo_lfc <=0 & de_res_plot$vitro_lfc <= 0,]
q4 <- de_res_plot[de_res_plot$vivo_lfc <=0 & de_res_plot$vitro_lfc > 0,]

qlist <- list(
    q1 = q1,
    q2 = q2,
    q3 = q3,
    q4 = q4
)
gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
gene_background <- rownames(de_res_plot)
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)
go_res <- lapply(qlist, function(x) {
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
saveRDS(go_res, paste0(mstatePath, test_text, "quadrant_go_res.rds"))
go_res_tbl <- lapply(go_res, function(x) {
    x@result
})
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, "quadrant_go_res.xlsx"))

# GO plot
glist <- list()
for(i in 1:length(go_res)) {
    include_go_simplified <- clusterProfiler::simplify(go_res[[i]], cutoff = .5)
    plot_df<- include_go_simplified@result[1:5,]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    glist[[i]] <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
        geom_bar(stat = "identity", aes(fill = -log10(qvalue))) + 
        #scale_x_discrete(position = 'top') + 
        coord_flip() +
        xlab("") + 
        theme(
            text=element_text(family = "Helvetica", size=8,color="black"),
            axis.text = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
        monocle:::monocle_theme_opts() 
}


pdf(paste0(mplotPath, test_text, "quadrant_go_bar.pdf"), width = 5, height = 3.5)
do.call(grid.arrange, c(glist[1:3], ncol=1))
dev.off()

# Save list and annotate TF and receptors
gene_csv <- read.csv("../hcc_cello/cellphonedb_data/gene_input.csv",stringsAsFactors = F)
receptor_csv <- read.csv("../hcc_cello/cellphonedb_data/protein_input.csv", stringsAsFactors = F)
receptor_gene <- gene_csv$gene_name[match(receptor_csv$uniprot[as.logical(receptor_csv$receptor)], gene_csv$uniprot)]

# Resource from http://humantfs.ccbr.utoronto.ca/download.php
tf_csv <- read.csv("../hcc_cello/human_tf/DatabaseExtract_v_1.01.csv", row.names = 1, stringsAsFactors = F)
human_tf<- tf_csv$HGNC.symbol[tf_csv$Is.TF. == "Yes"]

qgene_list <- lapply(qlist, function(x) {
    x$gene_id <- rownames(x)
    x <- x[order(x$vivo_mlgP, decreasing = T),]
    return(x[c("gene_id", "gene_name", "vivo_lfc", "vitro_lfc", "vivo_mlgP", "vitro_mlgP")])
})
saveRDS(qgene_list, paste0(mstatePath, test_text, "qgene_list.rds"))
write.xlsx(qgene_list, paste0(msheetPath, test_text, "quadrant_gene_list.xlsx"))

qgene_list <- readRDS(paste0(mstatePath, test_text, "qgene_list.rds"))
receptor_list <- lapply(qgene_list , function(x) {
    x[x$gene_name %in% receptor_gene, ,drop=F]
})
write.xlsx(receptor_list, paste0(msheetPath, test_text, "quadrant_receptor_list.xlsx"))

tf_list <- lapply(qgene_list , function(x) {
    x[x$gene_name %in% human_tf, ,drop=F]
})
write.xlsx(tf_list, paste0(msheetPath, test_text, "quadrant_tf_list.xlsx"))


# enzymatic gene, talk to Yuxuan for the list
drug_target <- read.delim("~/Documents/CHOP/NingProj/hcc_cello/TherapeuticTarget/version2019/AllDrugTargetInfo_forread.txt", header=FALSE, comment.char="#", stringsAsFactors = F)
drug_target <- drug_target$V3[drug_target$V2 == "UNIPROID"]
drug_target <- grep("HUMAN", drug_target, value = T)
#drug_target <- gsub("_HUMAN", "",drug_target)
drug_target <- unlist(lapply(drug_target, function(x) unlist(strsplit(x,"; |[-]|[/]"))))

# View(AllDrugTargetInfo_forread)
# 
# library(biomaRt)
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
# filters<-listFilters(ensembl)
# attributes <- listAttributes(ensembl)

write.table(data.frame(gene = drug_target), paste0(msheetPath, test_text, "drug_uniprot_toConvert.txt"), quote = F, row.names = F, col.names = F)

drug_target_converted <- read.table(paste0(msheetPath, "vivo_vs_vitro_drug_uniprot_Converted.txt"), stringsAsFactors = F, header = T)

drug_target_list <- lapply(qgene_list , function(x) {
    x[x$gene_name %in% drug_target_converted$To, ,drop=F]
})
write.xlsx(drug_target_list, paste0(msheetPath, test_text, "quadrant_drug_target_list.xlsx"))

# Distance to y=x
dist_to_expected <- de_res_plot$vivo_lfc - de_res_plot$vitro_lfc
hist(dist_to_expected, breaks = 100)

hl_cells <- dist_to_expected > 4
de_res_plot$hl_cell <- ifelse(hl_cells, T, F)
g1 <- ggplot(de_res_plot, aes(vivo_lfc, vitro_lfc)) +
    geom_point(size = 1, stroke = 0, aes(color = hl_cell)) + 
    geom_text(data = combine_de_res[text_f,], aes(label = gene_name), size = 3) +
    #stat_smooth(method = "lm", formula = "y~x+0", se = T) + 
    #stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon", bins = 100) + 
    scale_fill_viridis_c() + 
    theme_bw()+
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9)) 






# Add comparison to mouse

test_text <- "de3way_"

eset_mouse <- readRDS("../mouse_colon_cello/eset.rds")
clist_mouse <- readRDS("../mouse_colon_cello/clist.rds")

dep_vis <- clist_mouse$`Mouse epithelial [IFF]`
eset_mouse_epi <- eset_mouse[,dep_vis@idx]

eset_mouse_epi$de_group <- as.character(eset_mouse_epi$Dataset)
eset_mouse_epi$de_group[eset_mouse_epi$de_group %in% c("Colon1", "Colon2")] <- "Colon"

test_col = "de_group"
test_pairs <- list(
    "ApcMin"= c("Tumor ApcMin", "Colon"),
    "AOM" = c("Tumor AOM", "Colon"),
     "APKS" = c("Tumor APKS", "Colon")
)
# Derive and plot DEGs
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
max_count <- 500
eset_de <- eset_mouse_epi
feature_data <- fData(eset_de)
eset_de$de_group <- eset_de[[test_col]]

pdeg_list <- list()
for(i in 1:length(test_pairs)) {
    test_pair <- test_pairs[[i]]
    cur_name <- names(test_pairs)[i]
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <- which(pData(eset_de)$de_group == test_pair[2])
    group2_idx <- sample(group2_idx, min(length(group2_idx), max_count))
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
    pdeg_list[[cur_name]] <-prioritized_genes
}

saveRDS(pdeg_list, paste0(resultPath,test_text,"pdeg_list.rds"))

# de_list <- lapply(prioritized_genes, function(x) {
#     x %>% dplyr::filter(p_adj <= .01 & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
# })
# saveRDS(de_list, paste0(resultPath,test_text,"de_sig_genes.rds"))
# write.xlsx(de_list, paste0(resultPath,test_text,"de_sig_genes.xlsx"))



vivo_de_res <- readRDS(paste0(resultPath, "de3way_tumor_vs_colon_prioritized_genes.rds"))
vitro_de_res <- readRDS(paste0(resultPath, "de3way_tumoroid_vs_organoid_prioritized_genes.rds"))
mouse_pdeg_list <- readRDS(paste0(resultPath,"de3way_mousetumor_vs_colon_pdeg_list.rds"))

use_g <- rownames(vivo_de_res$Tumor)
combine_de_res<- data.frame(vivo_lfc = vivo_de_res$Tumor[use_g, ]$log2fc, vivo_mlgP = -log10(vivo_de_res$Tumor[use_g, ]$p_adj), vitro_lfc = vitro_de_res$Tumoroid[use_g, ]$log2fc, vitro_mlgP = -log10(vitro_de_res$Tumoroid[use_g, ]$p_adj))
rownames(combine_de_res) <- use_g
combine_de_res$gene_name <-fData(eset_epi)$gene_short_name[match(use_g, rownames(eset_epi))]
HMD_HumanPhenotype <- read.delim("~/Documents/CHOP/HSC/Preprocess.VisCello.eht/data-raw/public_resource/HMD_HumanPhenotype.rpt", header=FALSE)
combine_de_res$mouse_gene_name <- mouse_to_human_symbol(gene_symbol = combine_de_res$gene_name, in.type = "hs", HMD_HumanPhenotype = HMD_HumanPhenotype)
sum(!is.na(combine_de_res$mouse_gene_name)) # 17151

combine_de_res$ApcMin_lfc <- mouse_pdeg_list$ApcMin$`Tumor ApcMin`$log2fc[match(combine_de_res$mouse_gene_name, mouse_pdeg_list$ApcMin$`Tumor ApcMin`$gene_name)]
combine_de_res$ApcMin_mlgP <- -log10(mouse_pdeg_list$ApcMin$`Tumor ApcMin`$p_adj[match(combine_de_res$mouse_gene_name, mouse_pdeg_list$ApcMin$`Tumor ApcMin`$gene_name)])

combine_de_res$AOM_lfc <- mouse_pdeg_list$AOM$`Tumor AOM`$log2fc[match(combine_de_res$mouse_gene_name, mouse_pdeg_list$AOM$`Tumor AOM`$gene_name)]
combine_de_res$AOM_mlgP <- -log10(mouse_pdeg_list$AOM$`Tumor AOM`$p_adj[match(combine_de_res$mouse_gene_name, mouse_pdeg_list$AOM$`Tumor AOM`$gene_name)])

combine_de_res$APKS_lfc <- mouse_pdeg_list$APKS$`Tumor APKS`$log2fc[match(combine_de_res$mouse_gene_name, mouse_pdeg_list$APKS$`Tumor APKS`$gene_name)]
combine_de_res$APKS_mlgP <- -log10(mouse_pdeg_list$APKS$`Tumor APKS`$p_adj[match(combine_de_res$mouse_gene_name, mouse_pdeg_list$APKS$`Tumor APKS`$gene_name)])

combine_de_res$mouse_lfc <- (combine_de_res$ApcMin_lfc+ combine_de_res$AOM_lfc+ combine_de_res$APKS_lfc)/3

saveRDS(combine_de_res, paste0(mstatePath, "de3way_combine_de_res.rds"))


plot_df <- combine_de_res[complete.cases(combine_de_res),]
g1 <- ggplot(plot_df, aes(vivo_lfc, vitro_lfc)) +
    geom_point(data = combine_de_res, aes(vivo_lfc, mouse_lfc), size = .5, stroke = 0) + 
    geom_abline(intercept = 0, slope = 1)+
    stat_smooth(method = "lm", formula = "y ~ x +0", se = F, color = "blue") + 
    stat_smooth(data = combine_de_res, aes(vivo_lfc, ApcMin_lfc), method = "lm", formula = "y ~ x +0", se = F, color = "red") + 
    stat_smooth(data = combine_de_res, aes(vivo_lfc, AOM_lfc), method = "lm", formula = "y ~ x +0", se = F, color = "green") +
    stat_smooth(data = combine_de_res, aes(vivo_lfc, AOM_lfc), method = "lm", formula = "y ~ x +0", se = F, color = "purple") +
    #stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon", alpha = .4) + 
    scale_fill_viridis_c() + 
    theme_bw()+
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9)) 
ggsave(paste0(mplotPath, test_text, "lfc_vivo_vitro.pdf"), g1, width = 4, height=2.8, units = "in", device = "pdf")

df <- plot_df[,c("vivo_lfc", "vitro_lfc")]
colnames(df) <- c("x", "y")
m <- lm(y ~ x+0, df)
summary(m)

df <- plot_df[,c("vivo_lfc", "ApcMin_lfc")]
colnames(df) <- c("x", "y")
m <- lm(y ~ x+0, df)
summary(m)

df <- plot_df[,c("vivo_lfc", "AOM_lfc")]
colnames(df) <- c("x", "y")
m <- lm(y ~ x+0, df)
summary(m)

df <- plot_df[,c("vivo_lfc", "APKS_lfc")]
colnames(df) <- c("x", "y")
m <- lm(y ~ x+0, df)
summary(m)

cor_res <- data.frame(vitro = cor(plot_df$vivo_lfc, plot_df$vitro_lfc),
                      ApcMin = cor(plot_df$vivo_lfc, plot_df$ApcMin_lfc),
                      AOM = cor(plot_df$vivo_lfc, plot_df$AOM_lfc),
                      APKS = cor(plot_df$vivo_lfc, plot_df$APKS_lfc)
            )

cor_res <- reshape2::melt(cor_res)
colnames(cor_res) <- c("variable", "Correlation")
g1 <- ggplot(cor_res)+
    geom_col(aes_string(x = "variable", y = "Correlation"), fill = "darkblue") +
    theme_classic()+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9)) 
ggsave(paste0(mplotPath, test_text, "cor_bar.pdf"), g1, width = 2.2, height=1.4, units = "in", device = "pdf")

deg_res <- data.frame(vitro = sum(combine_de_res$vivo_mlgP > 2 & combine_de_res$vivo_lfc > 0 & combine_de_res$vitro_mlgP>2 & combine_de_res$vitro_lfc > 0)/sum(combine_de_res$vivo_mlgP > 2 & combine_de_res$vivo_lfc > 0),
                      ApcMin = sum(plot_df$vivo_mlgP > 2 & plot_df$vivo_mlgP > 0 & plot_df$ApcMin_mlgP>2 & plot_df$ApcMin_lfc > 0)/sum(combine_de_res$vivo_mlgP > 2 & combine_de_res$vivo_mlgP > 0),
                      AOM = sum(plot_df$vivo_mlgP > 2 & plot_df$vivo_mlgP > 0 & plot_df$AOM_mlgP>2 & plot_df$AOM_lfc > 0)/sum(combine_de_res$vivo_mlgP > 2 & combine_de_res$vivo_mlgP > 0),
                      APKS = sum(plot_df$vivo_mlgP > 2 & plot_df$vivo_mlgP > 0 & plot_df$APKS_mlgP>2 & plot_df$APKS_lfc > 0)/sum(combine_de_res$vivo_mlgP > 2 & combine_de_res$vivo_mlgP > 0)
)

deg_res <- reshape2::melt(deg_res)
colnames(deg_res) <- c("variable", "Fraction")
g1 <- ggplot(deg_res)+
    geom_col(aes_string(x = "variable", y = "Fraction"), fill = "#006d2c") +
    theme_classic()+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9)) 
ggsave(paste0(mplotPath, test_text, "deg_recap_bar.pdf"), g1, width = 2.2, height=1.4, units = "in", device = "pdf")


