



# Preprocess TCGA data
library(TCGAbiolinks)
library(MultiAssayExperiment)
library(maftools)
library(dplyr)
library(ComplexHeatmap)
library(openxlsx)
library(VisCello)

test_text <- "survival_analysis_"
mplotPath <- "~/Dropbox/ColonManuscript/subplots/"
msheetPath <- "~/Dropbox/ColonManuscript/sheets/"
mstatePath<- "~/Dropbox/ColonManuscript/states/"

scriptPath <- paste0("../hcc_cello/scripts/")
source(paste0(scriptPath, "read_excel.R"))
source(paste0(scriptPath, "iff_selection.R"))
source(paste0(scriptPath, "sigma_main.R"))
source(paste0(scriptPath, "plotFunction.R"))
source(paste0(scriptPath, "compute_dimR.R"))

coad.clinical <- read.delim("~/Documents/CHOP/NingProj/TCGA/TCGA-COAD.GDC_phenotype.tsv")
coad.exp <- GDCquery(project = "TCGA-COAD", 
                           data.category = "Transcriptome Profiling", 
                     workflow.type = "HTSeq - Counts",
                           data.type = "Gene Expression Quantification")
GDCdownload(coad.exp)
coad_se <- GDCprepare(query = coad.exp, summarizedExperiment = T)
colData(coad_se)$MSI_state <- coad.clinical$microsatellite_instability[match(colData(coad_se)$sample, coad.clinical$submitter_id.samples)]
saveRDS(coad_se, paste0(mstatePath, test_text, "coad_se.rds"))

read.clinical <- read.delim("~/Documents/CHOP/NingProj/TCGA/TCGA-READ.GDC_phenotype.tsv.gz")
read.exp <- GDCquery(project = "TCGA-READ", 
                     data.category = "Transcriptome Profiling", 
                     workflow.type = "HTSeq - Counts",
                     data.type = "Gene Expression Quantification")
GDCdownload(read.exp)
read_se <- GDCprepare(query = read.exp, summarizedExperiment = T)
colData(read_se)$MSI_state <- read.clinical$microsatellite_instability[match(colData(read_se)$sample, read.clinical$submitter_id.samples)]
saveRDS(read_se, paste0(mstatePath, test_text, "read_se.rds"))

coad_se <- readRDS(paste0(mstatePath, "survival_analysis_", "coad_se.rds"))
read_se <- readRDS(paste0(mstatePath, "survival_analysis_", "read_se.rds"))



combineSE <- function(x, y, ..., nomatch = NA, use.mcols = FALSE) {
    
    args <- unname(list(x, y, ...))
    
    # Check that each object has unique colnames
    colnames <- unlist(lapply(args, colnames))
    if (anyDuplicated(colnames)) {
        stop(paste0("Cannot combine ", class(args[[1]]), " objects with ", 
                    "duplicate 'colnames'"))
    }
    
    # Check that each object has the same assays
    an <- lapply(args, assayNames)
    if (any(sapply(an, function(x, y) any(is.na(match(x, y))), 
                   y = an[[1]]))) {
        stop(paste0("All ", class(args[[1]]), " objects must have ", 
                    "identical 'assayNames'"))
    }
    
    # Combine rowRanges or NAMES slot
    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        # NOTE: rowRanges(x) includes elementMetadata(x) since 
        #       elementMetadata(x) = elementMetadata(rowRanges(x)).
        #       This is the reason for the use.mcols if-else block.
        if (use.mcols) {
            all_rowRanges <- do.call(c, lapply(args, rowRanges))
        } else {
            all_rowRanges <- do.call(c, lapply(args, function(x) {
                rr <- rowRanges(x)
                elementMetadata(rr) <- NULL
                rr
            }))
        }
        rowRanges <- unique(all_rowRanges)
        nr <- length(rowRanges)
    } else {
        if (any(vapply(args, function(x) is.null(x@NAMES), logical(1)))) {
            stop("Cannot combine ", class(args[[1]]), " objects with ", 
                 "'NAMES' set to NULL")
        }
        all_NAMES <- do.call(c, lapply(args, function(x) x@NAMES))
        NAMES <- unique(all_NAMES)
        nr <- length(NAMES)
    }
    
    # Combine colData
    colData <- do.call(rbind, lapply(args, colData))
    
    # Combine elementMetadata
    if (use.mcols) {
        if (is(args[[1L]], "RangedSummarizedExperiment")) {
            elementMetadata <- mcols(rowRanges)
        } else {
            # IDEA: Create DataFrame with all_NAMES/all_rowRanges in one 
            #       column and elementMetadata in others, unique-ify, and 
            #       check that the number of unique rows equals nr.
            # WARNING: This will be slow for large SummarizedExperiment0 
            #          objects
            elementMetadata <- unique(
                cbind(DataFrame(all_NAMES), 
                      do.call(rbind, lapply(args, elementMetadata))))
            # Now drop the all_NAMES column (assumes it is the first column)
            elementMetadata <- elementMetadata[, -c(1L), drop = FALSE]
            # Sanity check
            if (nrow(elementMetadata) != nr) {
                stop(paste0("'elementMetadata' must match across ", 
                            class(args[[1]]), " objects"))
            }
        }
    } else {
        elementMetadata <- DataFrame()
        elementMetadata@nrows <- nr
    }
    
    # Create assays of the correct dimension (fill with 'nomatch')
    # First, create the empty combined assay using the appropriate 
    # storage.mode (guessed from the storage.mode of the assay in the 
    # first sample).
    nomatch <- lapply(seq_along(an[[1]]), function(i) {
        storage.mode(nomatch) <- storage.mode(assay(args[[1]], 
                                                    i, 
                                                    withDimnames = FALSE))
        nomatch
    })
    assays <- lapply(nomatch, function(nm) {
        matrix(nm, nrow = nr, ncol = length(colnames))
    })
    names(assays) <- an[[1]]
    
    # NOTE: I suspect that there are faster and more efficient ways to 
    # combine the assays, perhaps at the C-level.
    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        for (j in seq_along(args)) {
            ol <- findOverlaps(args[[j]], rowRanges, type = "equal")
            for (i in seq_along(assays)) {
                assays[[i]][subjectHits(ol),
                            match(colnames(args[[j]]), colnames)] <- 
                    assay(args[[j]], i, withDimnames = FALSE)
            }
        }
    } else {
        for (j in seq_along(args)) {
            ol <- match(args[[j]]@NAMES, NAMES)
            for (i in seq_along(assays)) {
                assays[[i]][ol, match(colnames(args[[j]]), colnames)] <- 
                    assay(args[[j]], i, withDimnames = FALSE)
            }
        }
    }
    assays <- Assays(assays)
    
    # Combine metadata
    metadata <- do.call(c, lapply(args, metadata))
    
    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        # No need to replace elementMetadata slot since it is part of 
        # rowRanges.
        BiocGenerics:::replaceSlots(args[[1L]], 
                                    rowRanges = rowRanges,
                                    colData = colData, 
                                    assays = assays,
                                    metadata = metadata)
    } else {
        BiocGenerics:::replaceSlots(args[[1L]],
                                    NAMES = NAMES,
                                    colData = colData,
                                    assays = assays,
                                    metadata = metadata,
                                    elementMetadata = elementMetadata)
    }
}



combined_se <- combineSE(coad_se, read_se)
rowData(combined_se) <- rowData(coad_se)
combined_se <- combined_se[,combined_se$MSI_state != "YES"]
saveRDS(combined_se, paste0(mstatePath, test_text, "combined_se.rds"))

#coad_se <- TCGAanalyze_Preprocessing(coad_se)



# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
# save(dataBRCA, geneInfo , file = "dataGeneExpression.rda")

# normalization of genes
#raw_count <- assay(combined_se, "HTSeq - Counts")
combined_se <- readRDS(paste0(mstatePath, test_text, "combined_se.rds"))
dataPrep <- TCGAanalyze_Preprocessing(object = combined_se,
                                      cor.cut = 0.5,
                                      datatype = "HTSeq - Counts")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

saveRDS(dataNorm, paste0(mstatePath, test_text, "dataNorm.rds"))
# boxplot(dataPrep, outline = FALSE)
# 
# boxplot(dataNorm, outline = FALSE)
# 
# dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                   method = "quantile", 
#                                   qnt.cut =  0.25)   
# samplesDown <- colnames(dataPrep)
# 
# dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                   typesample = "TP")
# 
# dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                   typesample = "NT")
# 
# dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP],
#                             mat2 = dataFilt[,dataSmNT],
#                             Cond1type = "Normal",
#                             Cond2type = "Tumor",
#                             fdr.cut = 0.01 ,
#                             logFC.cut = 1,
#                             method = "glmLRT")
# saveRDS(dataDEGs, paste0(resultPath, test_text, "dataDEG_glmLRT.rds"))
# 
# dataDEGs <- readRDS(paste0(resultPath, test_text, "dataDEG_glmLRT.rds"))
# dataDEGs <- dataDEGs %>% tibble::rownames_to_column("gene_id") %>% tibble::add_column(gene_name = rowData(combined_se)$external_gene_name[match(rownames(dataDEGs), rowData(combined_se)$ensembl_gene_id)], .after = 1) %>% dplyr::arrange(FDR)
# deg_res <- list("Normal" = dataDEGs[dataDEGs$logFC > 0, ], "Tumor" = dataDEGs[dataDEGs$logFC < 0, ])
# write.xlsx(deg_res, paste0(resultPath, test_text,"deg_res.xlsx"))


# Plot expression of selected genes
dataNorm <- readRDS(paste0(mstatePath, test_text, "dataNorm.rds"))

plot_g <- c("IL1B", "C1QC","SPP1", "CD44")
plot_id <- rownames(rowData(combined_se))[match(plot_g, rowData(combined_se)$external_gene_name)]
ncol <- 2

colData(combined_se)$Tumor_or_Normal <- ifelse(colData(combined_se)$shortLetterCode == "NT", "Normal", "Tumor")

plot_group <- "Tumor_or_Normal"
plot_meta <- as.data.frame(colData(combined_se))[,plot_group,drop=F]
plot_meta$Tumor_or_Normal <- as.character(plot_meta$Tumor_or_Normal)
g_exprs <- as.data.frame(log10(t(as.matrix(dataNorm[plot_id, match(rownames(plot_meta), colnames(dataNorm))]))+1))
colnames(g_exprs) <- plot_g

# save as violin plots
glist <- list()
for(g in plot_g) {
    print(g)
    gene_values <- as.matrix(g_exprs[,g,drop=F])
    ecut <- quantile(gene_values[gene_values > 0], .99)
    noise <- rnorm(n = length(x = gene_values[,g]))/1e+05
    gene_values[,g] <- gene_values[,g] + noise
    colnames(gene_values) <- "expression_level"
    df <- cbind(gene_values, plot_meta)
    df$Tumor_or_Normal <- ifelse(df$Tumor_or_Normal == "Tumor", "T", "N")
    glist[[g]]<- ggplot(df, aes_string(x=plot_group, y="expression_level")) +
        geom_violin(aes_string(fill = plot_group, color= plot_group), trim = T, scale = "width") + 
        geom_jitter(size = .1, alpha = .5) + 
        stat_compare_means(aes_string(group = plot_group), method = "t.test", label = "p.signif", size = 4, vjust = .5) + 
        #scale_color_manual(values = use_color) +
        #scale_fill_manual(values = use_color)+
        #ylim(0, ecut) +
        theme_classic() + 
        guides(alpha = F, fill=F,color=F) + 
        xlab("Cell type") +
        ylab(g) +
        theme(text=element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(angle=90, hjust=1, size=8), 
              axis.text.y = element_text(size=8), 
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              legend.margin=margin(0,0,0,0))
}
g1<-do.call(grid.arrange,c(glist, ncol = 2))
ggsave(paste0(mplotPath, test_text, "selected_expression.pdf"), g1, width = 3, height=3, units = "in", device = "pdf")



# Survival analysis

# Survival plot of individual genes
plot_g <- c("IL1B", "C1QC","SPP1", "CD44")
plot_meta <- colData(combined_se)[rownames(g_exprs),]
glist <- list()
for(g in plot_g) {
    print(g)
    gene_values <- as.matrix(g_exprs[,g,drop=F])
    hicut <- quantile(gene_values, .55)
    locut <- quantile(gene_values, .45)
    plot_meta[[g]] <- as.character(ifelse(gene_values >= hicut, "high", ifelse(gene_values <= locut, "low", "medium")))
    plot_df <- plot_meta[plot_meta[[g]] != "medium", ]
    TCGAanalyze_survival(plot_df,
                         g,
                         main = "TCGA Set\n COAD",
                         filename = paste0(mplotPath,test_text, g,"_survival.pdf"),
                         height = 5, width=4)
}

plot_meta$CD44_SPP1 <- as.character(ifelse(plot_meta$CD44 == "high" & plot_meta$SPP1 == "high", "CD44_hi_SPP1_hi",
                              ifelse(plot_meta$CD44 == "high" & plot_meta$SPP1 == "low", "CD44_hi_SPP1_lo",
                                     ifelse(plot_meta$CD44 == "low" & plot_meta$SPP1 == "high", "CD44_lo_SPP1_hi",
                                            ifelse(plot_meta$CD44 == "low" & plot_meta$SPP1 == "low", "CD44_lo_SPP1_lo", "else")))))

# TCGAanalyze_survival(plot_meta[plot_meta$CD44_SPP1 != "else",],
#                      "CD44_SPP1",
#                      main = "TCGA Set\n COAD & READ",
#                      filename = paste0(mplotPath,test_text,"CD44_SPP1_survival.pdf"),
#                      height = 8, width=6)

plot_meta$C1QC_SPP1 <- as.character(ifelse(plot_meta$C1QC == "high" & plot_meta$SPP1 == "low", "C1QC_hi_SPP1_lo",
                                                         ifelse(plot_meta$C1QC == "low" & plot_meta$SPP1 == "high", "C1QC_lo_SPP1_hi", "else")))

plot_meta$s <- as.numeric(grepl("dead|deceased", plot_meta$vital_status, ignore.case = TRUE))

test_col = "C1QC_SPP1" 
# test_col = "CD44_SPP1"
plot_df <-plot_meta[plot_meta[[test_col]] != "else",]
require("survival")
library("survminer")
notDead <- is.na(plot_df$days_to_death)
if (any(notDead == TRUE)) {
    plot_df[notDead, "days_to_death"] <- plot_df[notDead, "days_to_last_follow_up"]
}
plot_df$type <- as.factor(plot_df[, test_col])
plot_df$s <- grepl("dead|deceased", plot_df$vital_status, ignore.case = TRUE)
f.m <- formula(survival::Surv(as.numeric(plot_df$days_to_death), event = plot_df$s) ~ plot_df$type)
fit <- survfit(as.formula(paste0("Surv(days_to_death, s) ~ ", test_col)), data = plot_df)
g1<-ggsurvplot(fit,pval=TRUE)
g1 <- g1$plot + theme(legend.position = "right")
ggsave(paste0(mplotPath, test_text, test_col, "_surv.pdf"), g1, width = 7, height=4, units = "in", device = "pdf")




# Survival based on signature score
de_list <- read_excel_allsheets(paste0("~/Dropbox/ColonManuscript/sheets/", "macrophages in vivo [20210411]_IL1B_vs_SPP1_vs_C1QC_2021-04-13_de_significant.xlsx"))
max_gnum = 100
sig_deg <- lapply(de_list, function(x) {
    x$gene_short_name[1:min(max_gnum, nrow(x))]
})

library(AUCell)
expr_norm <- dataNorm
rownames(expr_norm) <- rowData(combined_se)$external_gene_name[match(rownames(expr_norm), rownames(rowData(combined_se)))]
cells_rankings_hc <- AUCell_buildRankings(expr_norm)
cur_thresh <- 0.2
cells_AUC <- AUCell_calcAUC(sig_deg, cells_rankings_hc, aucMaxRank=nrow(cells_rankings_hc)*cur_thresh)
matrixAUC <- getAUC(cells_AUC)

plot_meta <- cbind(colData(combined_se)[colnames(matrixAUC),], t(matrixAUC))

plot_g <- rownames(matrixAUC)
# save as violin plots
glist <- list()
library(ggpubr)
plot_group <- "Tumor_or_Normal"
for(g in plot_g) {
    print(g)
    gene_values <- as.matrix(t(matrixAUC)[,g,drop=F])
    colnames(gene_values) <- "state_signature"
    df <- as.data.frame(cbind(gene_values, plot_meta))
    glist[[g]]<- ggplot(df, aes_string(x=plot_group, y="state_signature")) +
        geom_violin(aes_string(fill = plot_group, color= plot_group), trim = T, scale = "width") + 
        geom_jitter(size = .1, alpha = .5) + 
        stat_compare_means(aes_string(group = plot_group), method = "t.test", label = "p.signif", size = 4, vjust = .5) + 
        theme_classic() + 
        guides(alpha = F, fill=F,color=F) + 
        xlab("Cell type") +
        ylab(g) +
        theme(text=element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(angle=90, hjust=1, size=8), 
              axis.text.y = element_text(size=8), 
              legend.text=element_text(size=8),
              axis.text = element_text(size=8),
              legend.margin=margin(0,0,0,0))
}
g1<-do.call(grid.arrange,c(glist, ncol = 3))
ggsave(paste0(mplotPath, test_text, "mac_state_activity.pdf"), g1, width = 4, height=2, units = "in", device = "pdf")


glist <- list()
for(g in plot_g) {
    print(g)
    gene_values <- as.matrix(t(matrixAUC)[,g,drop=F])
    hicut <- quantile(gene_values, .55)
    locut <- quantile(gene_values, .45)
    plot_meta[[paste0(g, "_statelevel")]] <- as.character(ifelse(gene_values >= hicut, "high", ifelse(gene_values <= locut, "low", "medium")))
    plot_df <- plot_meta[plot_meta[[paste0(g, "_statelevel")]] != "medium", ]
    TCGAanalyze_survival(plot_df,
                         paste0(g, "_statelevel"),
                         main = "TCGA Set\n COAD",
                         filename = paste0(mplotPath,test_text, g,"_statelevel_survival1.pdf"),
                         height = 5.5, width=4.5)
}

g <- "C1QC"; 
#g = "SPP1"; 
#g = "IL1B"
test_col = paste0(g, "_statelevel")
plot_df <-plot_meta[plot_meta[[paste0(g, "_statelevel")]] != "medium", ]
require("survival")
library("survminer")
notDead <- is.na(plot_df$days_to_death)
if (any(notDead == TRUE)) {
    plot_df[notDead, "days_to_death"] <- plot_df[notDead, "days_to_last_follow_up"]
}
plot_df$type <- as.factor(plot_df[, test_col])
plot_df$s <- grepl("dead|deceased", plot_df$vital_status, ignore.case = TRUE)
f.m <- formula(survival::Surv(as.numeric(plot_df$days_to_death), event = plot_df$s) ~ plot_df$type)
fit <- survfit(as.formula(paste0("Surv(days_to_death, s) ~ ", test_col)), data = plot_df)
g1<-ggsurvplot(fit,pval=TRUE)
g1 <- g1$plot + theme(legend.position = "right")
ggsave(paste0(mplotPath, test_text, test_col, "_surv.pdf"), g1, width = 6, height=3.5, units = "in", device = "pdf")


