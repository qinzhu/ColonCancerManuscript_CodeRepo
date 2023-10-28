

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
library(ggrastr)

test_text <- "emt_signature_"
library(AUCell)
cur_vis <- clist$`Coculture experiment final set [IFF]`
cur_eset <- eset[,cur_vis@idx]
cur_eset$Cluster <- cur_vis@pmeta$Cluster
emt_genes <- read.delim("~/Documents/CHOP/NingProj/VisCello.cc/data-raw/public_resource/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v7.5.1.tsv", row.names=1)
emt_genes <- strsplit(emt_genes["MAPPED_SYMBOLS",], ",")
names(emt_genes) <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
expressed_g <- which(rowMeans(exprs(cur_eset) > 0) > .01)
expr_hc <- as.matrix(exprs(cur_eset[expressed_g,]))
rownames(expr_hc) <- make.names(fData(cur_eset)$gene_short_name[match(rownames(expr_hc), fData(cur_eset)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_hc)

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(emt_genes, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "coculcutre_0.2cells_AUC.rds"))
cells_AUC <- readRDS(paste0(mstatePath, test_text, "coculcutre_0.2cells_AUC.rds"))
plot_df <- as.data.frame(t(getAUC(cells_AUC)))
plot_df <- cbind(plot_df, pData(cur_eset))

epi_proj <- plot_df[!plot_df$Cluster %in% c(3,4,6) & plot_df$SampleType != "Macrophage",]
stype_color <- c("O" = "#3690c0", "O+M" = "#41ab5d", "T" = "#fd8d3c", "T+M" = "#e7298a")
stype_mapping = c("Organoid" = "O",  "Organoid+Macrophage" = "O+M",  "Tumoroid" = "T", "Tumoroid+Macrophage" = "T+M")
epi_proj$SampleType2 <- stype_mapping[epi_proj$SampleType]

g1 <- ggplot(epi_proj, aes(x = SampleType2, y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) +
    geom_boxplot(aes(fill = SampleType2), outlier.shape=NA) + 
    scale_fill_manual(values = stype_color) +
    theme(text=element_text(family = "Helvetica", size=7),
          legend.text=element_text(size=7),
          axis.text = element_text(size=7),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) + 
    ylim(c(0, 0.25))+
    facet_wrap(~Patient) + 
    theme_bw() 
ggsave(paste0(mplotPath, test_text, "epi_emt_boxplot.pdf"), g1, width = 4.5, height=3, units = "in", device = "pdf")


# T test on emt signature
res_list <- list()
test_group = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
for(pp in c("28", "86", "92")) {
    print(pp)
    pp_proj <- epi_proj[epi_proj$Patient == pp,]
    if(pp %in% c("86", "92")) {
        tres1 <- t.test(pp_proj[[test_group]][pp_proj$SampleType2 == "T+M"], pp_proj[[test_group]][pp_proj$SampleType2 == "T"], alternative = "greater")
        tres2 <- t.test(pp_proj[[test_group]][pp_proj$SampleType2 == "O+M"], pp_proj[[test_group]][pp_proj$SampleType2 == "O"], alternative = "greater")
        res_list[[pp]] <- rbind(
            data.frame(comparison = "T+M_vs_T", p.value = tres1$p.value), 
            data.frame(comparison = "O+M_vs_O", p.value = tres2$p.value) 
        )
    } else {
        tres1 <- t.test(pp_proj[[test_group]][pp_proj$SampleType2 == "T+M"], pp_proj[[test_group]][pp_proj$SampleType2 == "T"], alternative = "greater")
        res_list[[pp]] <- data.frame(comparison = "T+M_vs_T", p.value = tres1$p.value)
    }
}
write.xlsx(res_list, paste0(msheetPath, test_text, "emt_t_res_list.xlsx"),  rowNames = F)
write.xlsx(as.data.frame.matrix(table(epi_proj$SampleType2, epi_proj$Patient)),  paste0(msheetPath, test_text, "emt_cellcount.xlsx"), rowNames=T)

write.xlsx(epi_proj[, c("Patient", "SampleType2", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")], paste0(msheetPath, test_text, "epi_emt_box.xlsx"),  rowNames = T)

# Which (cells) and genes contribute to the signature - perform DE on cells that show high EMT response (top 10 percentile) compared to baseline (no exposure to mac)
eset_epi <- eset[,rownames(epi_proj)]
pData(eset_epi) <- epi_proj
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
de_res <- list()
for(pp in c("28", "86", "92")) { 
    print(pp)
    test_col = "SampleType2"
    for(test_pair in list(c("O+M", "O"), c("T+M", "T"))) {
        eset_de <- eset_epi[, eset_epi$Patient == pp & eset_epi$SampleType2 %in% test_pair]
        if(ncol(eset_de) == 0) next
        id_col <- "gene_id"
        name_col <- "gene_short_name"
        test_method = "sSeq"
        max_count <- 1000
        feature_data <- fData(eset_de)
        
        eset_de$de_group <- eset_de[[test_col]]
        group1_idx <- which(eset_de$de_group == test_pair[1] & eset_de$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION > quantile(eset_de$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,.9))
        group2_idx <-  sample(which(eset_de$de_group == test_pair[2]), min(max_count, length(which(eset_de$de_group == test_pair[2]))))
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
        names(prioritized_genes) <- paste0(pp,"_", test_pair)
        saveRDS(prioritized_genes, paste0(mstatePath,test_text,pp,"_",paste0(test_pair, collapse = "_vs_"),"_prioritized_genes.rds"))
        de_list <- lapply(prioritized_genes, function(x) {
            x %>% dplyr::filter(p_adj <= .01) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj, proportion1, proportion2)
        })
        de_res[[paste0(pp,"_",paste0(test_pair, collapse = "_vs_"))]] <- de_list
        write.xlsx(de_list, paste0(msheetPath,test_text,pp,"_",paste0(test_pair, collapse = "_vs_"),"_de_sig_genes.xlsx"))
    }
}
saveRDS(de_res, paste0(mstatePath,test_text,"emt_de_res.rds"))


de_res <- readRDS(paste0(mstatePath,test_text,"emt_de_res.rds"))
deg_tbls <- lapply(de_res, function(x)x[[1]])
lapply(deg_tbls, nrow)

# Exclude macrophage specific genes from the DEG set due to potential doublet contaminants, 
# macrophage_deg <- read_excel_allsheets(paste0(msheetPath, "Coculture experiment final set [IFF]_macrophage_vs_epi_2022-05-30_de_significant.xlsx"))
# macrophage_deg <- macrophage_deg$macrophage$gene_id
# deg_tbls_mac_excluded <- lapply(deg_tbls, function(x) {
#     x[!x$gene_id %in% macrophage_deg, ]
# })

deg_tbls_prop_filtered <- lapply(deg_tbls, function(x) {
    x[x$proportion1 >= 0.3,, drop=F]
})
lapply(deg_tbls_prop_filtered, nrow)

# deg_tumor <- deg_tbls_prop_filtered[grep("T[+]M", names(deg_tbls_prop_filtered), value = T)]
# deg_organoid <- deg_tbls_prop_filtered[grep("O[+]M", names(deg_tbls_prop_filtered), value = T)]

#deg_tumor_union <- Reduce(union,lapply(deg_tumor, function(x) x$gene_id))

## TO INSERT, Pathway analysis, after running doublet removal and background removal - seems that interfere a lot with the analysis
# Alternatively, require at least 10% epithelial cells express this gene (so it's not likely due to doublet), this may be a much easier solution
gene_id_type = "ENSEMBL"
gene_background <- rownames(eset_epi)
orgdb = "org.Hs.eg.db"
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)
go_res <- lapply(deg_tbls_prop_filtered, function(x) {
    include_g.df <- bitr(x$gene_id, fromType = gene_id_type,
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
saveRDS(go_res, paste0(mstatePath, test_text, "coculture_emt_go_res_prop_filtered.rds"))
go_res_tbl <- lapply(go_res, function(x) {
    if(!is.null(x)) x@result else NA
})
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, "coculture_emt_go_res_prop_filtered.xlsx"))

# GO plot
plot_num = 10
glist <- list()
for(i in 1:length(go_res)) {
    cur_name = names(go_res)[i]
    include_go_simplified <- clusterProfiler::simplify(go_res[[i]], cutoff = .5)
    plot_df<- include_go_simplified@result[1:min(nrow(include_go_simplified@result),plot_num),]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    glist[[i]] <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
        geom_bar(stat = "identity", fill = "#2b8cbe") + 
        #scale_x_discrete(position = 'top') + 
        xlab("GO (biological process)") + 
        coord_flip() +
        theme(
            text=element_text(family = "Helvetica", size=8,color="black"),
            axis.text.x = element_text(family = "Helvetica",lineheight = .8, size=8,color="black"),
            axis.text.y = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
        monocle:::monocle_theme_opts() +
        ggtitle(cur_name)
}

pdf(paste0(mplotPath, test_text, "go_bar.pdf"), width = 5, height = 10)
do.call(grid.arrange, c(glist, ncol=1))
dev.off()

# # plot significance as heatmap
# all_go_terms <- lapply(go_res_tbl, function(x) {
#     x$Description[x$qvalue <= 0.05]
# })
# intersect_go_terms <- Reduce(intersect,  all_go_terms)
# union_go_terms <- Reduce(union,  all_go_terms)
# 
# union_go_terms_sum <- sapply(go_res_tbl, function(x) {
#     union_go_terms %in% x$Description[x$qvalue <= 0.05]
# })
# rownames(union_go_terms_sum) <- union_go_terms
# union_go_terms_sum <- union_go_terms_sum[order(rowSums(union_go_terms_sum), decreasing = T), ]

gene_id_type = "ENSEMBL"
gene_background <- rownames(eset_epi)
orgdb = "org.Hs.eg.db"
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)
kegg_res <- lapply(deg_tbls_prop_filtered, function(x) {
    include_g.df <- bitr(rownames(x), fromType = gene_id_type,
                         toType = c("SYMBOL", "ENTREZID"),
                         OrgDb = orgdb)
    include_g_kegg<- enrichKEGG(gene        = include_g.df$ENTREZID,
                                universe      = bg.df$ENTREZID,
                                organism = "hsa",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
    return(include_g_kegg)
})
saveRDS(kegg_res, paste0(mstatePath, test_text, "coculture_emt_kegg_res_prop_filtered.rds"))
kegg_res <- readRDS(paste0(mstatePath, test_text, "coculture_emt_kegg_res_prop_filtered.rds"))
kegg_res_tbl <- lapply(kegg_res, function(x) {
    x@result
})
write.xlsx(kegg_res_tbl, paste0(msheetPath, test_text, "coculture_emt_kegg_res_prop_filtered.xlsx"))

# KEGG plot
plot_num = 10
glist <- list()
for(i in 1:length(kegg_res)) {
    cur_name = names(kegg_res)[i]
    plot_df<- kegg_res[[i]][1:min(nrow(kegg_res[[i]]),plot_num),]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    glist[[i]] <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
        geom_bar(stat = "identity", fill = "#2b8cbe") + 
        #scale_x_discrete(position = 'top') + 
        xlab("KEGG pathway") + 
        coord_flip() +
        theme(
            text=element_text(family = "Helvetica", size=8,color="black"),
            axis.text.x = element_text(family = "Helvetica",lineheight = .8, size=8,color="black"),
            axis.text.y = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
        monocle:::monocle_theme_opts()  +
        ggtitle(cur_name)
}

pdf(paste0(mplotPath, test_text, "kegg_bar.pdf"), width = 5, height = 10)
do.call(grid.arrange, c(glist, ncol=1))
dev.off()




# Intersect with EMT gene list
emt_deg_tbls <- lapply(deg_tbls, function(x) {
    x[x$gene_name %in% emt_genes$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, ]
})

write.xlsx(emt_deg_tbls, paste0(msheetPath,test_text,"emt_deg_tbls.rds"))

show_genes <- c("VIM", "CD44")
cur_gene <- show_genes[[1]]
#cur_gene <- show_genes[[2]]
epi_proj$gene_expr <- eset@assayData$norm_exprs[which(fData(eset)$gene_short_name == cur_gene),rownames(epi_proj)]
g1 <- ggplot(epi_proj, aes_string(x = "SampleType2", y = "gene_expr", fill = "SampleType2")) +
    geom_boxplot(outlier.colour = NA, alpha = .5)+
    geom_jitter(aes(color = SampleType2), size = 0.3, stroke =0) + 
    scale_fill_manual(values = stype_color)+
    scale_color_manual(values = stype_color)+
    guides(color = F) + 
    theme(text=element_text(family = "Helvetica", size=9),
          axis.text = element_text(family = "Helvetica", size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.margin=margin(0,0,0,0))+
    ggtitle(cur_gene)+
    ylim(c(0, 2.5))+
    xlab(NULL)+
    theme_bw() + 
    facet_wrap(~Patient)

pdf(paste0(mplotPath, test_text,cur_gene, "_gene_expr", ".pdf"),  width = 4, height=2.5)
g1
dev.off()

write.xlsx(epi_proj[,c("Patient", "SampleType2", "gene_expr")], paste0(msheetPath,test_text,"vim_expr.xlsx"), rowNames=T)


# ANOVA on VIM/CD44
for(pp in c("28", "86", "92")) {
    anova_list<-list()
    tukey_list <- list()
    plot_df <- epi_proj[which(epi_proj$Patient == pp),]
    for(show_p in "gene_expr") {
        df <- plot_df[,c(show_p, "SampleType2")]
        colnames(df) <- c("expr", "sample_type")
        lm_res <- lm(expr ~ sample_type, data = df)
        aov_res <- aov(lm_res)
        anova_list[[show_p]] <- data.frame(unclass(summary(aov_res)), check.names = FALSE, stringsAsFactors = FALSE)
        tukey_res <- TukeyHSD(aov_res)
        tukey_res <- as.data.frame(tukey_res$sample_type)
        tukey_res$significance <- ifelse(tukey_res[["p adj"]] <= 0.0001, "****", ifelse(tukey_res[["p adj"]] <= 0.001, "***", ifelse(tukey_res[["p adj"]] <= 0.01, "**", ifelse(tukey_res[["p adj"]] <= 0.05, "*", "ns"))))
        tukey_list[[show_p]] <- tukey_res
    }
    write.xlsx(anova_list, paste0(msheetPath, test_text, cur_gene, "_anova_results_", pp, ".xlsx"),  row.names = T)
    write.xlsx(tukey_list, paste0(msheetPath, test_text, cur_gene,"_tukey_list_", pp, ".xlsx"), row.names = T)
}


# t test on VIM/CD44
res_list <- list()
test_group = "gene_expr"
for(pp in c("28", "86", "92")) {
    print(pp)
    pp_proj <- epi_proj[epi_proj$Patient == pp,]
    if(pp %in% c("86", "92")) {
        tres1 <- t.test(pp_proj[[test_group]][pp_proj$SampleType2 == "T+M"], pp_proj[[test_group]][pp_proj$SampleType2 == "T"], alternative = "greater")
        tres2 <- t.test(pp_proj[[test_group]][pp_proj$SampleType2 == "O+M"], pp_proj[[test_group]][pp_proj$SampleType2 == "O"], alternative = "greater")
        res_list[[pp]] <- rbind(
            data.frame(comparison = "T+M_vs_T", p.value = tres1$p.value), 
            data.frame(comparison = "O+M_vs_O", p.value = tres2$p.value) 
        )
    } else {
        tres1 <- t.test(pp_proj[[test_group]][pp_proj$SampleType2 == "T+M"], pp_proj[[test_group]][pp_proj$SampleType2 == "T"], alternative = "greater")
        res_list[[pp]] <- data.frame(comparison = "T+M_vs_T", p.value = tres1$p.value)
    }
}


##### EMT signature in in vivo data ######

dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

expressed_g <- which(rowMeans(exprs(eset_epi) > 0) > .01)
expr_hc <- as.matrix(exprs(eset_epi[expressed_g,]))
rownames(expr_hc) <- make.names(fData(eset_epi)$gene_short_name[match(rownames(expr_hc), fData(eset_epi)$gene_id)])
cells_rankings_hc <- AUCell_buildRankings(expr_hc)

cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(emt_genes, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
saveRDS(cells_AUC,paste0(mstatePath, test_text, "epi_0.2cells_AUC.rds"))

cells_AUC <- readRDS(paste0(mstatePath, test_text, "epi_0.2cells_AUC.rds"))
plot_df <- as.data.frame(t(getAUC(cells_AUC)))
plot_df <- cbind(plot_df, pData(eset_epi))
stype_color <- c("Colon" = "#1f78b4", "Organoid" = "#33a02c", "Tumor" = "#e31a1c", "Tumoroid" = "#ff7f00")
plot_df <- plot_df[plot_df$SampleType != "LiverMet", ]

stype_color <- c("Colon" = "#1f78b4", "Organoid" = "#33a02c", "Tumor" = "#e31a1c", "Tumoroid" = "#ff7f00")
g1 <- ggplot(plot_df, aes(x = Dataset, y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) +
    geom_boxplot(aes(fill = SampleType), outlier.shape=NA) + 
    ylim(0,0.25) + 
    theme_classic() +
    scale_fill_manual(values = stype_color) +
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "epi_emt_score.pdf"), g1, width = 5, height=5, units = "in", device = "pdf")

# Corrleate fraction of macrophage in tumor with EMT signature
cur_vis <- clist$`All in vivo cells [20210411, IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
compos<-as.data.frame.matrix(table(cur_proj$SampleType, cur_proj$Cell_type))
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

compos<-as.data.frame.matrix(table(cur_proj$Dataset, cur_proj$Cell_type))
compos <- compos[!grepl("Liver", rownames(compos)) & !rownames(compos) %in% c("Tumor_44", "Tumor_49", "Tumor_81", "Tumor_84", "Tumor_87"),]
compos_rel <- compos / (compos$`Epithelial(Normal)` + compos$`Epithelial(Tumor)`)


breaksList <- seq(0,1,by=.1)
pheatmap(compos_rel[,show_ctype], cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2",display_numbers =T, number_format ="%.2f",fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)

# Intersect dataset
dataset_shared<- intersect(names(table(plot_df$Dataset)), rownames(compos_rel))
names(table(plot_df$Dataset))[!names(table(plot_df$Dataset)) %in% dataset_shared]

median_emt_sig <- plot_df %>% group_by_at("Dataset") %>% summarize_at(c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), median)
colnames(compos_rel) <- make.names(colnames(compos_rel))
median_emt_sig <- cbind(median_emt_sig, compos_rel[match(median_emt_sig$Dataset, rownames(compos_rel)),])

median_emt_sig <- median_emt_sig[complete.cases(median_emt_sig),]

ctype_color <- c(
    "Fibroblast" = "#cab2d6", 
    "Myofibroblast" = "#fb9a99",
    "Endothelial"= "#33a02c",
    "Myeloid" = "#6a3d9a",
    "T.cell" ="#a6cee3",
    "B.cell" = "#fdbf6f",
    "Plasma.cell" = "#ff7f00",
    "Mast.cell" = "#ae017e"
)

plot_list <- list()
library(ggpmisc)
for(ctype in names(ctype_color)) {
    plot_df <-  median_emt_sig[,c("Dataset", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ctype)]
    colnames(plot_df) <- c("Dataset", "y","x")
    
    plot_df$log_x <- log(plot_df$x)
    my.formula <- y ~ x
    plot_list[[ctype]] <- ggplot(plot_df, aes(log_x, y)) +
        geom_point(size = 8) + 
        geom_text(aes(label = Dataset), color = "white", size = 1.5) + 
        stat_smooth(method = "lm") +
        stat_poly_eq(formula = my.formula,
                     eq.with.lhs = "italic(hat(y))~`=`~",
                     aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                     parse = TRUE, size = 2)  + 
        theme_bw() +
        ggtitle(ctype) + 
        xlab("log(relative fraction)") + 
        ylab("HALLMARK_EMT")
}

g1<-do.call(grid.arrange,c(plot_list, ncol = 4))
ggsave(paste0(mplotPath, test_text, "compos_cor_emt.pdf"), g1, width = 10, height=5, units = "in", device = "pdf")



