
library(VisCello)
savePath <- "../hcc_final_cello/"
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

library(phateR)

eset <- readRDS("../hcc_final_cello/eset.rds")
clist <- readRDS("../hcc_final_cello/clist.rds")

test_text <- "epi_analysis_20210408_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

cur_proj <- cbind(dep_vis@proj$`UMAP-2D [10PC, IFG]`, pData(eset_epi))
cur_proj <- cur_proj[cur_proj$SampleType != "LiverMet",]

label_data <- cur_proj %>% group_by_at("Dataset") %>% summarize_at(c("UMAP_1", "UMAP_2"), median)

dataset_color <- get_factor_color(sort(unique(cur_proj$Dataset)), "Accent")
#dataset_color <- c(dataset_color[5:length(dataset_color)], dataset_color[1:4])
names(dataset_color) <- sort(unique(cur_proj$Dataset))
library(ggrepel)
library(ggrastr)
g1<-ggplot(cur_proj) +
    geom_point_rast(aes(UMAP_1, UMAP_2, color=Dataset),size=.5, stroke=0) +
    geom_text_repel(data = label_data, aes(UMAP_1, UMAP_2, label=Dataset), size = 3) + 
    scale_color_manual(values = dataset_color) + 
    guides(color=guide_legend(title=NULL, ncol=1, keyheight=.6, override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts() 
ggsave(paste0(mplotPath, test_text, "dataset_umap.pdf"), g1, width = 5.2, height=3.6, units = "in", device = "pdf")


write.xlsx(cur_proj[,c("UMAP_1", "UMAP_2", "Dataset")], paste0(msheetPath,test_text,"dataset_umap.xlsx"), rowNames = T)



# proj_3d <- cbind(epi_PHATE_3d$embedding, pData(eset_epi))
# #proj_3d <- proj_3d[proj_3d$SampleType != "LiverMet",]
# p1 <- plotly::plot_ly(proj_3d,
#                       x = as.formula(paste0("~", "PHATE1")), y = as.formula(paste0("~", "PHATE2")), z = as.formula(paste0("~", "PHATE3")),
#                       text=proj_3d$Dataset,
#                       hoverinfo="text",
#                       marker = list(size = 1), color = ~Dataset, colors = get_numeric_bin_color(sort(unique(proj_3d$Dataset)), "gg_color_hue")) %>%
#     plotly::add_markers(
#         #opacity=alpha_manual[proj$alpha]
#     ) 
# htmlwidgets::saveWidget(p1, paste0(mplotPath,test_text, "phate3d_res.html"))


# Hierarchical clustering
use_sample <- rownames(cur_proj)
eset_epi_ifg <- eset_epi[,use_sample]
# use_ifg <- ifg_select(data = exprs(eset_epi_ifg), cluster = eset_epi_ifg$Dataset, cluster_min_cell_num = 30, min_cluster_expr_fraction = .1, gini_cut_qt = .8, compute_go = F, filePath = mplotPath, fileName = test_text, orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL) 
# eset_epi_ifg <- eset_epi_ifg[use_ifg, ]
expressed_g <- rowMeans(exprs(eset_epi_ifg)>0) > .05
eset_epi_ifg <- eset_epi_ifg[expressed_g, ]
saveRDS(eset_epi_ifg, paste0(mstatePath, test_text, "eset_epi_ifg.rds"))
library(data.table)
# mean expression vector for each dataset
norm_expr_df <- as.data.table(t(as.matrix(eset_epi_ifg@assayData$norm_exprs)))
norm_expr_df$Dataset <- eset_epi_ifg$Dataset
norm_expr_mean <- norm_expr_df[,lapply(.SD,mean,na.rm=TRUE),by=Dataset]

expr_mean_df <- as.data.frame(t(norm_expr_mean[,-1]))
colnames(expr_mean_df) <- norm_expr_mean$Dataset

# library(lsa)

cor_dist <- proxy::dist(t(expr_mean_df),method="correlation")
hc<- hclust(cor_dist, method = "ward.D2")

library(ggtree)
#hc <- hclust(dist(t(expr_mean_df)))
clus <- cutree(hc, 4)
d <- data.frame(label = names(clus), Dataset = names(clus))
p <- ggtree(hc)
g1<-p %<+% d + 
    #layout_dendrogram() + 
    scale_x_reverse()+
    geom_tippoint(size=3, shape=21, aes(fill=Dataset), stroke = 0, show.legend = F) + 
    scale_fill_manual(values = dataset_color) + 
    geom_tiplab(hjust=0, offset = -0.1, show.legend=F, align = T, size=3)  + 
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9))
    # theme(
    #     plot.margin = margin(2, 2, 2, 2, "cm")
    # )
ggsave(paste0(mplotPath, test_text, "dataset_clust_expressedg.pdf"), g1, width = 1, height=3.5, units = "in", device = "pdf")

cor_mtx <- cor(as.matrix(expr_mean_df))
pdf(paste0(mplotPath, test_text, "cor_heatmap_expressedg.pdf"), width = 3.7, height =4)
pheatmap(cor_mtx, cluster_rows = hc, cluster_cols = hc, show_colnames = T, show_rownames = F, color = get_numeric_color("RdBu"), fontsize = 8)
dev.off()

write.xlsx(expr_mean_df, paste0(msheetPath,test_text,"expr_mean_df.xlsx"), rowNames = T)




# Distribution cell to cell distance within each cluster
set.seed(2020)
cell_count_dataset <- table(eset_epi_ifg$Dataset)
use_cell_count <- 500
library(lsa)
cos_list<-lapply(names(cell_count_dataset), function(x) {
    cur_eset <- eset_epi_ifg[, eset_epi_ifg$Dataset == x]
    print(dim(cur_eset))
    cur_expr <- as.matrix(cur_eset@assayData$norm_exprs[, sample(1:ncol(cur_eset), min(use_cell_count, ncol(cur_eset)), replace = F)])
    cos_mtx <- lsa::cosine(cur_expr)
    cos_mtx[lower.tri(cos_mtx, diag = T)] <- NA
    return(cos_mtx)
})
names(cos_list) <- names(cell_count_dataset)
saveRDS(cos_list, paste0(mstatePath, test_text, "cos_list.rds"))

cos_list <- readRDS(paste0(mstatePath, test_text, "cos_list.rds"))
cos_plot <- lapply(1:length(cos_list), function(i) {
    x = cos_list[[i]]
    df<-data.frame(cos = x[!is.na(x)])
    df$dataset <- names(cos_list)[i]
    print(head(df))
    return(df)
})

cos_plot <- do.call(rbind, cos_plot)
cos_plot$sample_type <- sapply(strsplit(cos_plot$dataset, "_"), function(x) x[1])
g1 <- ggplot(cos_plot, aes(x = dataset, y = 1-cos, fill = dataset)) + 
    geom_boxplot(outlier.shape=NA, position = position_dodge(preserve = "single"), lwd=.3, show.legend = F) + 
    scale_fill_manual(values = dataset_color) + 
    #facet_grid(~sample_type, scales = "free_x", space='free') + 
    #coord_flip()+
    ylim(c(0,.7))+
    theme_classic() + 
    theme(text=element_text(family = "Helvetica", size=9),
          legend.text=element_text(size=9),
          axis.text = element_text(size=9),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(mplotPath, test_text, "cosine_het_boxplot.pdf"), g1, width = 4, height=2.3, units = "in", device = "pdf")

write.xlsx(cos_plot, paste0(msheetPath,test_text,"cosine_het_boxplot.xlsx"), rowNames = T)
ds_entries <- 50000
cos_plot_ds <- cos_plot[sample(1:nrow(cos_plot), ds_entries), ]
write.xlsx(cos_plot_ds, paste0(msheetPath,test_text,"cosine_het_boxplot_ds.xlsx"), rowNames = T)


test_pair <- c("Tumor_08", "Tumoroid_08")
# test_pair <- c("Tumor_24", "Tumoroid_24")
# test_pair <- c("Tumor_28", "Tumoroid_28")
group1_idx <- which(cos_plot$dataset == test_pair[1])
group2_idx <-  which(cos_plot$dataset == test_pair[2])
x <- 1-cos_plot$cos
g1_val <- x[group1_idx]
g2_val <- x[group2_idx]
t.test(g1_val, g2_val, alternative = "greater")




use_sample <- table(sapply(strsplit(names(table(cos_plot$dataset)), "_"), function(x) x[2]))
use_sample <- names(which(use_sample >= 2))
cos_plot$patient <- sapply(strsplit(cos_plot$dataset, "_"), function(x) x[2])

compare_pairs <- list(
    c("Colon", "Tumor"),
    c("Colon", "Organoid"), 
    c("Tumor", "Tumoroid"),
    c("Organoid", "Tumoroid")
)

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

sampletype_color <- gg_color_hue(4)
names(sampletype_color) <- unique(cos_plot$sample_type)

glist <- list()
for(i in 1:length(compare_pairs)) {
    cur_pair <- compare_pairs[[i]]
    cos_plot_filtered <- cos_plot[cos_plot$sample_type %in% cur_pair,]
    glist[[i]] <- ggplot(cos_plot_filtered, aes(x = patient, y = cos, fill = sample_type)) + 
        geom_boxplot(outlier.shape=NA, position = position_dodge(preserve = "single",width=1), lwd=.3, show.legend = T) + 
        scale_fill_manual(values = sampletype_color) + 
        #facet_grid(~sample_type, scales = "free_x", space='free') + 
        #coord_flip()+
        theme_classic() + 
        theme(text=element_text(family = "Helvetica", size=9),
              legend.text=element_text(size=9),
              axis.text = element_text(size=9),
              axis.text.x = element_text(angle = 90, hjust = 1),
              legend.margin=margin(0,0,0,0),
              #legend.box.margin=margin(-10,-10,-10,-10),
              plot.margin = unit(c(0,0.4,0,0.1), "cm"))
}
g1<-do.call(grid.arrange,c(glist, ncol = 2))
ggsave(paste0(mplotPath, test_text, "cosine_het_boxplot_grouped.pdf"), g1, width = 9, height=5, units = "in", device = "pdf")





### DE to derive core epithelial signature

### Not used, because not considering background contamination
# test_text <- "tumor_vs_colon_de_"
# dep_vis <- clist$`Epithelial cells [wt OGN, IFF]`
# eset_epi <- eset[, dep_vis@idx]
# bg_count <- 1000
# bg_idx <- sample(which(eset_epi$Dataset %in% c("Colon_07", "Colon_08", "Colon_09")), bg_count, replace = F)
# max_count <- 500
# eset_de <- eset_epi[, c(which(eset_epi$SampleType == "Tumor"), bg_idx)]

### 
set.seed(2020)
test_text <- "tumor_vs_colon_de_"
eset_epi <- eset_epi[,eset_epi$SampleType %in% c("Colon", "Tumor")]
bg_count <- 1000
bg_idx <- sample(which(eset_epi$Dataset %in% c("Colon_07", "Colon_08", "Colon_09")), bg_count, replace = F)

tumor_idx <- which(eset_epi$Cell_type == "Epithelial(Tumor)")

max_count <- 500
eset_de <- eset_epi[,c(tumor_idx, bg_idx)]

eset_de$de_group <- ifelse(eset_de$Cell_type == "Epithelial(Tumor)", eset_de$Dataset, "Background")
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

# 
# deg_list<- readRDS(paste0(mstatePath,"tumor_vs_colon_de_","de_sig_genes.rds"))
# library(UpSetR)
# deg_res <- lapply(deg_list, function(x) {
#     x$gene_name
# })
# deg_res <- deg_res[names(deg_res)!= "Tumor_44"]
# int_ind <- as.list(names(deg_res))
# int_miss1 <- lapply(int_ind, function(x){
#     names(deg_res)[!names(deg_res) %in% x]
# })
# use_intersections <- c(list(names(deg_res)), int_miss1, int_ind)
# pdf(paste0(mplotPath, test_text, "deg_upset.pdf"), width = 4, height=3.5)
# upset(fromList(deg_res), nsets = 10, intersections = use_intersections, order.by = "degree")
# dev.off()
# 
# deg_res2 <- deg_res[c("Tumor_08", "Tumor_24", "Tumor_28")]
# int_ind <- as.list(names(deg_res2))
# int_miss1 <- lapply(int_ind, function(x){
#     names(deg_res2)[!names(deg_res2) %in% x]
# })
# use_intersections <- c(list(names(deg_res2)), int_miss1, int_ind)
# pdf(paste0(mplotPath, test_text, "deg_upset2.pdf"), width = 4, height=3.5)
# upset(fromList(deg_res2), nsets = 10, intersections = use_intersections, order.by = "degree")
# dev.off()
# 
# 


### DE to derive core tumoroid signature
set.seed(2020)
test_text <- "tumoroid_vs_organoid_de_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]
eset_de <- eset_epi[, which(eset_epi$SampleType %in% c("Tumoroid", "Organoid"))]

test_pairs <- list(
    c("Tumoroid_08", "Organoid_08"),
    c("Tumoroid_24", "Organoid_24"),
    c("Tumoroid_28", "Organoid_28")
)
bg_count <- 1000
max_count <- 500

eset_de$de_group <- eset_de$Dataset
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

pdeg_list<- list()
for(test_pair in test_pairs){
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), max_count))
    group2_idx <-  which(pData(eset_de)$de_group == test_pair[2])
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
    pdeg_list[[test_pair[1]]] <- prioritized_genes[[1]]
}

saveRDS(pdeg_list, paste0(mstatePath, test_text, "pdeg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))


library(UpSetR)
deg_res <- lapply(deg_list, function(x) {
    x$gene_name
})
int_ind <- as.list(names(deg_res))
int_miss1 <- lapply(int_ind, function(x){
    names(deg_res)[!names(deg_res) %in% x]
})
use_intersections <- c(list(names(deg_res)), int_miss1, int_ind)
pdf(paste0(mplotPath, test_text, "deg_upset.pdf"), width = 4, height=2.2)
upset(fromList(deg_res), nsets = 10, intersections = use_intersections, order.by = "degree")
dev.off()


# Compare core signatures between tumor and tumoroid
deg_tumor <- readRDS(paste0(mstatePath,"tumor_vs_colon_de_","de_sig_genes.rds")) # UPDATE!
deg_tumor <- lapply(deg_tumor, function(x) x$gene_name)
deg_tumor  <- deg_tumor[names(deg_tumor) != "Tumor_44"]
deg_tumoroid <- readRDS(paste0(mstatePath,"tumoroid_vs_organoid_de_","de_sig_genes.rds"))
deg_tumoroid <- lapply(deg_tumoroid, function(x) x$gene_name)

get_intersect_members <- function (x, ...){
    require(dplyr)
    require(tibble)
    x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
    n <- names(x)
    x %>% rownames_to_column() -> x
    l <- c(...)
    a <- intersect(names(x), l)
    ar <- vector('list',length(n)+1)
    ar[[1]] <- x
    i=2
    for (item in n) {
        if (item %in% a){
            if (class(x[[item]])=='integer'){
                ar[[i]] <- paste(item, '>= 1')
                i <- i + 1
            }
        } else {
            if (class(x[[item]])=='integer'){
                ar[[i]] <- paste(item, '== 0')
                i <- i + 1
            }
        }
    }
    do.call(filter_, ar) %>% column_to_rownames() -> x
    return(x)
}
fromList2 <- function (input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
}
int_ind <- as.list(names(deg_tumor))
int_miss1 <- lapply(int_ind, function(x){
    names(deg_tumor)[!names(deg_tumor) %in% x]
})
use_intersections <- c(list(names(deg_tumor)), int_miss1)
core_tumor <- unique(unlist(lapply(use_intersections, function(x) {
    rownames(get_intersect_members(fromList2(deg_tumor), x))
})))
saveRDS(core_tumor, paste0(mstatePath, "core_tumor_signature.rds"))

int_ind <- as.list(names(deg_tumoroid))
int_miss1 <- lapply(int_ind, function(x){
    names(deg_tumoroid)[!names(deg_tumoroid) %in% x]
})
use_intersections <- c(list(names(deg_tumoroid)), int_miss1)
core_tumoroid <- unique(unlist(lapply(use_intersections, function(x) {
    rownames(get_intersect_members(fromList2(deg_tumoroid), x))
})))

library(VennDiagram)
pk_venn <- list("Tumor" = core_tumor, "Tumoroid" = core_tumoroid)
pdf(paste0(mplotPath, test_text, "_venn_coreg_veg.pdf"), width =3, height =3)
venn.plot <- venn.diagram(pk_venn, NULL, fill=c("Tumor" = "#3288bd", "Tumoroid" = "#d53e4f"), alpha=rep(0.5,length(pk_venn)), cex = 2, cat.fontface=4, category.names=names(pk_venn))
grid.draw(venn.plot)
dev.off()



# Heatmap of expression of core tumor signatures
test_text <- "core_tumor_sig_heatmap_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]
core_tumor <-  readRDS(paste0(mstatePath, "core_tumor_signature.rds"))
eset_plot <- eset_epi[match(core_tumor, fData(eset_epi)$gene_short_name), eset_epi$Dataset %in% c("Tumor_08", "Tumor_24", "Tumor_28", "Tumoroid_08", "Tumoroid_24", "Tumoroid_28")]
max_cell = 100
plot_sample <- unlist(lapply(names(table(eset_plot$Dataset)), function(x) {
    cur_idx <- which(eset_plot$Dataset == x)
    sample(cur_idx, min(max_cell, length(cur_idx)), replace = F)
}))
eset_plot <- eset_plot[, plot_sample]
eset_plot <- eset_plot[rowMeans(exprs(eset_plot) > 0) > .2, ]
plot_meta <- pData(eset_plot)
plot_meta <- plot_meta[, c("SampleType", "Dataset"), drop=F]
plot_meta <- plot_meta[order(plot_meta[["SampleType"]], plot_meta$Dataset),]
non_na_cells_ordered <- rownames(plot_meta)
value <- eset_plot@assayData$norm_exprs[,non_na_cells_ordered]
rownames(value) <- fData(eset_plot)$gene_short_name[match(rownames(value), rownames(eset_plot))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

pdf(paste0(mplotPath, test_text, "core_sig", ".pdf"), width = 6, height=6)
pheatmap::pheatmap(value,
                   color = get_numeric_color("viridis"),
                   cluster_rows = T, cluster_cols = F,
                   #clustering_distance_rows = distfun1, clustering_method = distfun1,
                   show_rownames = T,
                   show_colnames = FALSE, annotation_row = NULL,
                   annotation_col = plot_meta, annotation_names_row = FALSE,
                   annotation_names_col = FALSE, 
                   fontsize = 4)
dev.off()







####################### Update 2020/10/20 ########################

set.seed(2020)
test_text <- "paired_tumoroid_vs_tumor_de_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]
test_pairs <- list(
    c("Tumoroid_08", "Tumor_08"),
    c("Tumoroid_24", "Tumor_24"),
    c("Tumoroid_28", "Tumor_28")
)
eset_de <- eset_epi[, which(eset_epi$Dataset %in% unlist(test_pairs))]

g1_count <- 1000
g2_count <- 1000

eset_de$de_group <- eset_de$Dataset
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
feature_data <- fData(eset_de)

pdeg_list<- list()
for(test_pair in test_pairs){
    print(test_pair)
    group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
    group1_idx <- sample(group1_idx, min(length(group1_idx), g1_count))
    group2_idx <-  which(pData(eset_de)$de_group == test_pair[2])
    group2_idx <- sample(group2_idx, min(length(group2_idx), g2_count))
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
    pdeg_list[[test_pair[1]]] <- prioritized_genes[[1]]
    pdeg_list[[test_pair[2]]] <- prioritized_genes[[2]]
}

saveRDS(pdeg_list, paste0(mstatePath, test_text, "pdeg_list.rds"))

deg_list <- lapply(pdeg_list, function(x){
    x %>% dplyr::filter(significant == TRUE & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(deg_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(deg_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))

lapply(deg_list, nrow)

# Venn diagram
library(VennDiagram)

show_group = "Tumor_"
res_list <- lapply(deg_list[grep(show_group,names(deg_list))], function(x) {
    x$gene_name
})
use_color <- c("#e41a1c", "#377eb8", "#4daf4a")
pdf(paste0(mplotPath, test_text, show_group, "venn_deg.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(res_list, NULL, cex = 2, fill=use_color, cat.fontface=4)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()

show_group = "Tumoroid_"
res_list <- lapply(deg_list[grep(show_group,names(deg_list))], function(x) {
    x$gene_name
})
pdf(paste0(mplotPath, test_text, show_group, "venn_deg.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(res_list, NULL, cex = 2, fill=use_color, cat.fontface=4)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()

show_group = "Tumor_" # show_group = "Tumoroid_"
res_list <- lapply(deg_list[grep(show_group,names(deg_list))], function(x) {
    x$gene_id
})
res_g <- Reduce(intersect, res_list)
gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
gene_background <- rownames(eset_de)
bg.df <- bitr(gene_background, fromType = gene_id_type,
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb)

include_g.df <- bitr(res_g, fromType = gene_id_type,
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = orgdb)
go_res<- enrichGO(gene        = include_g.df$ENTREZID,
                  universe      = bg.df$ENTREZID,
                  OrgDb         = orgdb,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

saveRDS(go_res, paste0(mstatePath, test_text, show_group, "go_res.rds"))
go_res_tbl <- go_res@result
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, show_group, "go_res.xlsx"))

# GO plot
include_go_simplified <- clusterProfiler::simplify(go_res, cutoff = .5)
plot_df<- include_go_simplified@result[1:min(nrow(include_go_simplified@result),15),]
plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
g1 <- ggplot(plot_df, aes(x = Description,y=-log10(qvalue))) +
    geom_bar(stat = "identity", aes(fill = -log10(qvalue))) + 
    #scale_x_discrete(position = 'top') + 
    coord_flip() +
    xlab("") + 
    theme(
        text=element_text(family = "Helvetica", size=8,color="black"),
        axis.text = element_text(family = "Helvetica",lineheight = .8, size=8,color="black")) + 
    ggtitle(paste0("Up-regulated in ", show_group)) + 
    monocle:::monocle_theme_opts() 

pdf(paste0(mplotPath, test_text, show_group, "go_bar.pdf"), width = 6, height = 3)
g1
dev.off()

# Receptors and ligands intersection
PairsLigRec <- read.delim("~/Documents/CHOP/NingProj/HAtlas/preprocess/public_resource/PairsLigRec.txt")
show_group = "Tumor_" 
receptor_list <- lapply(deg_list[grep(show_group,names(deg_list))], function(x) {
    x$gene_name[x$gene_name %in% as.character(PairsLigRec$Receptor.ApprovedSymbol)]
})
ligand_list <- lapply(deg_list[grep(show_group,names(deg_list))], function(x) {
    x$gene_name[x$gene_name %in% as.character(PairsLigRec$Ligand.ApprovedSymbol)]
})
pdf(paste0(mplotPath, test_text, show_group, "venn_receptor.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(receptor_list, NULL, cex = 2, fill=use_color, cat.fontface=4)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()

pdf(paste0(mplotPath, test_text, show_group, "venn_ligand.pdf"), width = 5, height = 5)
venn.plot <- venn.diagram(ligand_list, NULL, cex = 2, fill=use_color, cat.fontface=4)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
grid.draw(venn.plot)
dev.off()


cur_vis <- clist$`Coculture experiment [IFF]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, IFG]`, pData(eset)[cur_vis@idx,])
cur_proj$Cluster <- cur_vis@pmeta$Cluster
eset_epi_cc <- eset[,rownames(cur_proj[cur_proj$Cluster == 1,])]
eset_epi_cc <- eset_epi_cc[,!grepl("Macrophage", eset_epi_cc$Dataset)]

# For each receptor ligand, plot violin
show_group = "Tumor_" # show_group = "Tumoroid_"
res_list <- lapply(deg_list[grep(show_group,names(deg_list))], function(x) {
    x$gene_id[x$gene_name %in% as.character(PairsLigRec$Receptor.ApprovedSymbol)]
})
plot_g <- Reduce(intersect, res_list)
cur_proj <- pData(eset_epi_cc)
g_exprs <- as.data.frame(t(as.matrix(eset_epi_cc@assayData$norm_exprs[match(plot_g,fData(eset_epi_cc)$gene_id), ])))
use_color <- c(
    "Tumoroid_28_24hr" = "#fb9a99",  
    "Tumoroid_28_48hr" = "#e31a1c",
    "CoCulture_28_24hr" = "#b2df8a", 
    "CoCulture_28_48hr" = "#33a02c"
)
glist <- list()
for(g in plot_g) {
    print(g)
    gene_values <- as.matrix(g_exprs[,g,drop=F])
    gname <- fData(eset_epi_cc)$gene_short_name[match(g, fData(eset_epi_cc)$gene_id)]
    ecut <- quantile(gene_values[gene_values > 0], .99)
    #gene_values[gene_values == 0] <- NA
    noise <- rnorm(n = length(x = gene_values[,g]))/1e+05
    gene_values[,g] <- gene_values[,g] + noise
    #ecut <- max(quantile(gene_values,1),1)
    bpGroup = "Dataset"
    colnames(gene_values) <- "expression_level"
    df <- cbind(gene_values, cur_proj)
    glist[[g]]<- ggplot(df, aes_string(x=bpGroup, y="expression_level")) +
        geom_violin(aes_string(fill = bpGroup, color= bpGroup), trim = T, scale = "width") + 
        #geom_jitter(size = .5, aes_string(colour = bpGroup), alpha = .1)+
        scale_color_manual(values = use_color) +
        scale_fill_manual(values = use_color)+
        #ylim(0, ecut) +
        theme_classic() + 
        guides(alpha = F, fill=F,color=F) + 
        xlab("") +
        ylab(gname) +
        theme(text=element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(angle=45, hjust=1, size=8), 
              axis.text.y = element_text(size=8), 
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              legend.margin=margin(0,0,0,0))
}
g1<-do.call(grid.arrange,c(glist, ncol = 6))
ggsave(paste0(mplotPath, test_text, show_group, "gexprviolin.pdf"), g1, width = 10, height=10, units = "in", device = "pdf")








