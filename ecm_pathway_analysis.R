


set.seed(2020)
is_windows = FALSE
if(is_windows) {
    home_path = "C:/Users/qinzh/"
} else {
    home_path = "~"
}
mplotPath <- paste0(home_path, "/Dropbox/ColonManuscript/subplots/")
msheetPath <- paste0(home_path, "/Dropbox/ColonManuscript/sheets/")
mstatePath<- paste0(home_path, "/Dropbox/ColonManuscript/states/")

library(VisCello)
savePath <- "../hcc_final_cello/"
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
sheetPath <- paste0(savePath, "sheets/"); dir.create(sheetPath)
scriptPath <- paste0("../hcc_cello/", "scripts/")

clist <- readRDS("../hcc_final_cello/clist.rds")
eset <- readRDS("../hcc_final_cello/eset.rds")


# For in vivo comparison

test_text <- "ecm_tumor_vs_colon_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

# From de_two_way.R
de_list <- readRDS(paste0(mstatePath, "de_2way/","de2way_tumor_vs_colon_","de_sig_genes.rds"))

non_epi_sig<-readRDS(paste0(mstatePath,"de_non_epi_primary_","de_sig_genes.rds"))
#non_epi_sig_genes <- unique(unlist(lapply(non_epi_sig, function(x) x$gene_id)))
de_list<- lapply(de_list, function(x) {
    x[!x$gene_id %in% non_epi_sig$`Plasma cell`$gene_id, ]
})

gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
gene_background <- rownames(eset_epi)
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
go_res_tbl <- lapply(go_res, function(x) {
    x@result
})
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, "go_res.xlsx"))

# GO plot
glist <- list()
for(i in 1:length(go_res)) {
    include_go_simplified <- clusterProfiler::simplify(go_res[[i]], cutoff = .5)
    plot_df<- include_go_simplified@result[1:10,]
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


pdf(paste0(mplotPath, test_text, "go_bar.pdf"), width = 5, height = 3.5)
do.call(grid.arrange, c(glist, ncol=1))
dev.off()

# For in vitro comparison

test_text <- "ecm_tumoroid_vs_organoid_"

# From de_two_way.R
de_list <- readRDS(paste0(mstatePath, "de_2way/","de2way_tumoroid_vs_organoid_","de_sig_genes.rds"))

gene_id_type = "ENSEMBL"
orgdb = "org.Hs.eg.db"
gene_background <- rownames(eset)
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
go_res_tbl <- lapply(go_res, function(x) {
    x@result
})
write.xlsx(go_res_tbl, paste0(msheetPath, test_text, "go_res.xlsx"))

# GO plot
glist <- list()
for(i in 1:length(go_res)) {
    include_go_simplified <- clusterProfiler::simplify(go_res[[i]], cutoff = .5)
    plot_df<- include_go_simplified@result[1:10,]
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


pdf(paste0(mplotPath, test_text, "go_bar.pdf"), width = 5, height = 3.5)
do.call(grid.arrange, c(glist, ncol=1))
dev.off()


# Perform GO for the other two comparisons






library(matrinetR)
library(igraph)
library(qgraph)

test_text <- "matrinet_tumor_vs_colon_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

# From de_two_way.R
prioritized_genes <- readRDS(paste0(mstatePath, "de_2way/","de2way_tumor_vs_colon_","prioritized_genes.rds"))

# non_epi_sig<-readRDS(paste0(mstatePath,"de_non_epi_primary_","de_sig_genes.rds"))
# #non_epi_sig_genes <- unique(unlist(lapply(non_epi_sig, function(x) x$gene_id)))
# prioritized_genes<- lapply(prioritized_genes, function(x) {
#     x[!x$gene_id %in% non_epi_sig$`Plasma cell`$gene_id, ]
# })

matrisome_genes <- unique(c(matrixDB_edgelist$Gene1, matrixDB_edgelist$Gene2))

expressed_genes = names(which(rowMeans(exprs(eset_epi) >0) > .01))
print(length(expressed_genes))
expressed_gname = fData(eset)$gene_short_name[match(expressed_genes, rownames(fData(eset)))]

valid_genes = intersect(matrisome_genes, expressed_gname)
print(length(valid_genes))
valid_gids = fData(eset)$gene_id[match(valid_genes, fData(eset)$gene_short_name)]

de_tumor = prioritized_genes$Tumor[valid_gids,c('gene_id', 'gene_name', 'p', 'p_adj', 'log2fc', 'proportion1', 'proportion2')]
rownames(de_tumor) = de_tumor$gene_name
# Corrlation between genes
expr_mtx = eset_epi@assayData$norm_exprs[valid_gids, eset_epi$SampleType == "Tumor"]
rownames(expr_mtx) = valid_genes
cor_tumor = cor(t(as.matrix(expr_mtx)))

adj_mtx = matrixDB_adjacency[valid_genes,valid_genes]
use_genes = names(which(rowSums(adj_mtx != 0) != 0))
adj_mtx = adj_mtx[use_genes,use_genes]
diag(adj_mtx) = 0
cor_tumor = cor_tumor[use_genes,use_genes]
cor_tumor[is.na(cor_tumor)]=0
de_tumor = de_tumor[use_genes,]
#Select the genes that are available and valid in both cohorts 
pal1 = get_numeric_color("RdBu")

map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

edge_color1 = apply(cor_tumor, 1, function(x){map2color(x, pal1,limits=c(-.1,.1))})
node_color1 = map2color(de_tumor$log2fc, pal1,limits=c(-2,2))
#node_color1 = map2color(de_tumor$proportion1, pal1,limits=c(0,1))
# qgraph::qgraph(weighted_adjmat$prad[1:10,1:10], edge.color = use_color1[1:10,1:10], curve = -0.2, curveAll = TRUE)
# qgraph::qgraph(weighted_adjmat$prad, edge.color = use_color1, curve = -0.2, curveAll = TRUE)
pdf(paste0(mplotPath, test_text, "qgraph.pdf"), width = 4.5, height = 4)
qgraph::qgraph(adj_mtx, vsize = 2, esize = 2, borders = F, edge.color = edge_color1, color = node_color1, layout ="spring", labels = colnames(adj_mtx), label.cex=.2, label.scale = F)
dev.off()

# use_color2 = map2color(matrinet_GTEx$prad_GTEx$edge_df$cor_C, pal1,limits=c(-1,1))
# qgraph::qgraph(matrinet_GTEx$prad_GTEx$edge_df[,c(1,2)], edge.color = use_color2, curve = -0.2, curveAll = TRUE)
# qgraph::qgraph(matrinet_GTEx$prad_GTEx$edge_df[,c("Gene1", "Gene2", "cor_C")], edge.color = use_color2, directed=F, layout = "spring")
# qgraph::qgraph(matrinet_GTEx$prad_GTEx$edge_df[,c(1,2)], edge.color = use_color2, directed=F, layout = "spring")



test_text <- "matrinet_tumoroid_vs_organoid_"
dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

# From de_two_way.R
prioritized_genes <- readRDS(paste0(mstatePath, "de_2way/","de2way_tumoroid_vs_organoid_","prioritized_genes.rds"))

# non_epi_sig<-readRDS(paste0(mstatePath,"de_non_epi_primary_","de_sig_genes.rds"))
# #non_epi_sig_genes <- unique(unlist(lapply(non_epi_sig, function(x) x$gene_id)))
# prioritized_genes<- lapply(prioritized_genes, function(x) {
#     x[!x$gene_id %in% non_epi_sig$`Plasma cell`$gene_id, ]
# })

matrisome_genes <- unique(c(matrixDB_edgelist$Gene1, matrixDB_edgelist$Gene2))

expressed_genes = names(which(rowMeans(exprs(eset_epi) >0) > .01))
print(length(expressed_genes))
expressed_gname = fData(eset)$gene_short_name[match(expressed_genes, rownames(fData(eset)))]

valid_genes = intersect(matrisome_genes, expressed_gname)
print(length(valid_genes))
valid_gids = fData(eset)$gene_id[match(valid_genes, fData(eset)$gene_short_name)]

de_tumoroid = prioritized_genes$Tumoroid[valid_gids,c('gene_id', 'gene_name', 'p', 'p_adj', 'log2fc', 'proportion1', 'proportion2')]
rownames(de_tumoroid) = de_tumoroid$gene_name
# Corrlation between genes
expr_mtx = eset_epi@assayData$norm_exprs[valid_gids, eset_epi$SampleType == "Tumoroid"]
rownames(expr_mtx) = valid_genes
cor_tumoroid = cor(t(as.matrix(expr_mtx)))
cor_tumoroid[is.na(cor_tumoroid)]=0
adj_mtx = matrixDB_adjacency[valid_genes,valid_genes]
use_genes = names(which(rowSums(adj_mtx != 0) != 0))
adj_mtx = adj_mtx[use_genes,use_genes]
diag(adj_mtx) = 0
cor_tumoroid = cor_tumoroid[use_genes,use_genes]
de_tumoroid = de_tumoroid[use_genes,]
#Select the genes that are available and valid in both cohorts 
pal1 = get_numeric_color("RdBu")

map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

edge_color1 = apply(cor_tumoroid, 1, function(x){map2color(x, pal1,limits=c(-.1,.1))})
node_color1 = map2color(de_tumoroid$log2fc, pal1,limits=c(-2,2))

pdf(paste0(mplotPath, test_text, "qgraph.pdf"), width = 4.5, height = 4)
qgraph::qgraph(adj_mtx, vsize = 2, esize = 2, borders = F, edge.color = edge_color1, color = node_color1, layout ="spring", labels = colnames(adj_mtx), label.cex=.2, label.scale = F)
dev.off()




color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
pdf(paste0(mplotPath, "RdBucolorbar_qgraph.pdf"), width = 2, height = 3)
color.bar(pal1, -1)
dev.off()


write.xlsx(as.data.frame(cor_tumor), paste0(msheetPath, test_text, "cor_tumor.xlsx"), rowNames=T)
write.xlsx(as.data.frame(cor_tumoroid), paste0(msheetPath, test_text, "cor_tumoroid.xlsx"), rowNames=T)


# plot correlation level
library(ggpubr)
diag(cor_tumor) =NA
diag(cor_tumoroid) =NA
identical(rownames(cor_tumor), rownames(cor_tumoroid))
identical(colnames(cor_tumor), colnames(cor_tumoroid))
cor_tumor<- cor_tumor * adj_mtx
cor_tumoroid<- cor_tumoroid * adj_mtx
cor_tumor_df <- reshape2::melt(cor_tumor);cor_tumor_df <- cor_tumor_df[complete.cases(cor_tumor_df) & cor_tumor_df$value!=0,]
cor_tumoroid_df <- reshape2::melt(cor_tumoroid);cor_tumor_df <- cor_tumoroid_df[complete.cases(cor_tumoroid_df) & cor_tumoroid_df$value!=0,]
cor_cbn <- data.frame(cor = c(cor_tumor_df$value, cor_tumoroid_df$value), type = c(rep("Tumor", nrow(cor_tumor_df)), rep("Tumoroid", nrow(cor_tumoroid_df))))
stype_color <- c("Tumor" = "#e31a1c", "Tumoroid" = "#ff7f00")
g1 <-  ggplot(cor_cbn, aes_string(x = "type", y = "cor", fill = "type", color= "type")) +
    geom_boxplot(alpha = .7, outlier.shape = NA) +
    stat_compare_means(aes_string(group = "type"), method = "t.test", label = "p.signif", size = 4, hjust = -1, vjust = .5, method.args = list(alternative = "less")) +
    ylim(c(-0.1,0.1))+
    scale_color_manual(values = stype_color)+ 
    scale_fill_manual(values = stype_color)+ 
    theme_bw()+
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.4, "cm"),
          axis.text = element_text(size=8),
          #axis.text.x = element_text(angle=45, hjust=1, size=8),
          legend.margin=margin(0,0,0,0))+
    xlab(NULL) + 
    ylab("Correlation")
ggsave(paste0(mplotPath, test_text, "cor_box.pdf"), g1, width = 2.5, height=2, units = "in", device = "pdf")

library(ggrastr)











## Check fibroblast and myofibroblast

test_text <- "mesenchymal_invivo_"

cur_vis <- clist$`Mesenchymal cells in vivo [20221224]`
cur_proj <- cbind(cur_vis@proj$`UMAP-2D [30PC, global subset]`, pData(eset)[cur_vis@idx,])
cur_proj <- cur_proj[cur_proj$UMAP_1 < 8, ]
mctype_color <- c(
    "Fibroblast" = "#cab2d6", 
    "Myofibroblast" = "#fb9a99"
)
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="Cell_type", pal=mctype_color, size = .7, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5, rastr=T) + 
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
ggsave(paste0(mplotPath, test_text, "umap_ctype.pdf"), g1, width = 2.8, height=2, units = "in", device = "pdf")


plot_g <- c("ACTA2", "COL1A2")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(plot_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- plot_g
glist <- list()
for(g in plot_g) {
    print(g)
    cur_proj$gene_expr <- g_exprs[[g]]
    ecut <- quantile(cur_proj$gene_expr[cur_proj$gene_expr > 0], .975)
    if(ecut < 1) ecut = 1
    cur_proj$gene_expr[cur_proj$gene_expr > ecut] = ecut
    glist[[g]]<-plotProj(cur_proj, dim_col = c(1,2), group.by="gene_expr", pal="viridis", size = .5, legend.title = g, layer_zero = T, rastr=F) + guides(color = guide_colorbar(barwidth = 10, barheight = 1, title = g, raster = T)) +
        guides(color = guide_colorbar(
            barwidth = 4, barheight = .5, title = g,
            title.theme = element_text(size = 8),
            label.theme = element_text(size = 8))) + 
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        theme(text=element_text(family = "Helvetica", size=8),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.margin=margin(0,0,0,0))+
        monocle:::monocle_theme_opts()
}
#g1<-do.call(grid.arrange,c(glist, ncol = 3))
g1 <- glist[[1]]
ggsave(paste0(mplotPath, test_text, "gexpr_umap.pdf"), g1, width = 1.9, height=2.5, units = "in", device = "pdf")




show_g <- c("ACTA2", "CTHRC1", "BGN", "INHBA", "PDPN")
g_exprs <- as.data.frame(t(as.matrix(eset@assayData$norm_exprs[match(show_g,fData(eset)$gene_short_name), match(rownames(cur_proj), colnames(eset))])))
colnames(g_exprs) <- show_g
save_proj <- cbind(cur_proj[,c("UMAP_1", "UMAP_2", "Patient", "SampleType", "Cell_type")], g_exprs)
write.xlsx(save_proj, paste0(msheetPath,test_text,"save_proj.xlsx"), rowNames=T)


stype_color <- c(
    "Colon" = "#1f78b4",   
    "Tumor" = "#e31a1c",
    "Liver" = "#b2df8a", 
    "LiverMet" = "#6a3d9a"
)
cur_proj$SampleType <- factor(cur_proj$SampleType, levels = names(stype_color))
g1<-plotProj(cur_proj, dim_col = c(1,2), group.by="SampleType", pal=stype_color, size = .7, na.col = "lightgrey", legend=T,legend.size=2, legend.text.size = 8, legend.position = "right", legend.title = NULL, keyheight=.5) + 
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
ggsave(paste0(mplotPath, test_text, "umap_sampletype.pdf"), g1, width = 3.4, height=1.6, units = "in", device = "pdf")

g1<-ggplot(cur_proj, aes_string("UMAP_1", "UMAP_2")) + 
    geom_point(aes(color = SampleType),size = .7,stroke=0) +
    scale_color_manual(values = stype_color)+
    scale_fill_manual(values = stype_color)+
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    facet_wrap(~SampleType) + 
    theme_bw() +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          strip.text.x = element_text(size = 9),
          legend.margin=margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) 
ggsave(paste0(mplotPath, test_text, "sampletype_distribution.pdf"), g1, width = 3.3, height=2.4, units = "in", device = "pdf")


cell_count<-as.data.frame.matrix(table(cur_proj$SampleType, cur_proj$Cell_type))
cell_frac <- as.data.frame(cell_count / rowSums(cell_count))
cell_frac$Sample_type <- rownames(cell_count)
plot_df <- reshape2::melt(cell_frac)
colnames(plot_df) <- c("Sample_type", "Cell_type", "Count")
g1 = ggplot(plot_df, aes(x="", y = Count, fill = Cell_type)) + 
    geom_bar(stat = "identity", width=1) +
    coord_polar("y", start=0) + 
    facet_wrap(~Sample_type) +
    scale_fill_manual(values = mctype_color)+
    theme_bw() +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) 
ggsave(paste0(mplotPath, test_text, "sampletype_pie.pdf"), g1, width = 3.5, height=2.4, units = "in", device = "pdf")




####### Compare fibroblast from tumor vs normal colon ##########
test_text <- "fibroblast_tumor_vs_colon_"
test_col = "SampleType"
test_pair <- c("Tumor", "Colon")
cur_eset = eset[,rownames(cur_proj)]
eset_de <- cur_eset[, cur_eset$Cell_type == "Fibroblast"]

# test_text <- "myofibroblast_LiverMet_vs_Liver_"
# test_col = "SampleType"
# test_pair <- c("LiverMet", "Liver")
# cur_eset = eset[,rownames(cur_proj)]
# eset_de <- cur_eset[, cur_eset$Cell_type == "Fibroblast"]


# Derive and plot DEGs
source(paste0(scriptPath, "cellrangerrkit_fun.R"))
id_col <- "gene_id"
name_col <- "gene_short_name"
test_method = "sSeq"
max_count <- 300
feature_data <- fData(eset_de)

eset_de$de_group <- eset_de[[test_col]]
group1_idx <- which(pData(eset_de)$de_group == test_pair[1])
group1_idx <- sample(group1_idx, min(max_count, length(group1_idx)))
group2_idx <- which(pData(eset_de)$de_group == test_pair[2])
group2_idx <- sample(group2_idx, min(max_count, length(group2_idx)))
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

saveRDS(prioritized_genes, paste0(mstatePath,test_text,"prioritized_genes.rds"))
write.xlsx(prioritized_genes, paste0(msheetPath,test_text,"prioritized_genes.xlsx"))
de_list <- lapply(prioritized_genes, function(x) {
    x %>% dplyr::filter(p_adj <= .01 & log2fc > 1) %>% dplyr::select(gene_id, gene_name, common_mean, dispersion, log2fc, p, p_adj)
})
saveRDS(de_list, paste0(mstatePath,test_text,"de_sig_genes.rds"))
write.xlsx(de_list, paste0(msheetPath,test_text,"de_sig_genes.xlsx"))


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
plotg <- plotg[!is.na(plotg)]
value <- eset_de@assayData$norm_exprs[plotg,non_na_cells_ordered]
rownames(value) <- fData(eset_de)$gene_short_name[match(rownames(value), rownames(eset_de))]

value <- t(scale(t(as.matrix(value))))
limits <- c(-2,2)
if(!is.null(limits)) {
    value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]
}

pdf(paste0(mplotPath, test_text, "deg_", max_gene,".pdf"), width = 6, height=4)
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




library(ggpubr)
test_genes <- c("CTHRC1", "BGN", "INHBA", "PDPN")
eset_plot <- eset[,rownames(cur_proj)]
eset_plot <- eset_plot[,eset_plot$Cell_type == "Fibroblast" & eset_plot$SampleType %in% c("Tumor", "Colon")]
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
    glist[[g]]<- ggplot(df, aes_string(x="SampleType", y="expression_level")) +
        geom_violin(aes_string(fill = "SampleType", color= "SampleType"), trim = T, scale = "width", alpha = .7) +
        geom_jitter(size = .2, stroke=0, aes_string(colour = "SampleType"), alpha = 1, position=position_jitterdodge())+
        scale_color_manual(values = stype_color) +
        scale_fill_manual(values = stype_color)+
        #ylim(0, ecut) +
        theme_classic() +
        ggtitle(g)+
        guides(alpha = F, fill=F,color=F) +
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

g1<-do.call(grid.arrange,c(glist, ncol = 2))
ggsave(paste0(mplotPath, test_text, "violin_gexpr_violin.pdf"), g1, width = 2, height=3, units = "in", device = "pdf")


