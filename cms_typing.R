


library(VisCello)
library(CMScaller)

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

test_text <- "cms_typing_"
library(CMScaller)
# Make pseudo bulk for tumor
eset_tumor <- eset[, eset$SampleType %in% c("Tumor", "LiverMet") & eset$Patient!= "40"]
datasets <- names(table(eset_tumor$Dataset))
bulk_expr <- sapply(datasets, function(dset) {
    print(dset)
    cur_eset <- eset_tumor[, eset_tumor$Dataset == dset]
    rowSums(exprs(cur_eset))
})

eset_bulk <- new("ExpressionSet",
            assayData = assayDataNew( "environment", exprs=bulk_expr))
fData(eset_bulk) <- fData(eset)
# Convert to entrez id
fData(eset_bulk)$entrez <- fData(crcTCGAsubset)$entrez[match(fData(eset_bulk)$gene_short_name, fData(crcTCGAsubset)$symbol)]
eset_bulk <- eset_bulk[!is.na(fData(eset_bulk)$entrez) & !duplicated(fData(eset_bulk)$entrez),]
rownames(eset_bulk) <- fData(eset_bulk)$entrez
pdf(paste0(mplotPath, test_text, "res_plot.pdf"), width = 3.5, height = 3.5)
res <- CMScaller(eset_bulk, RNAseq=TRUE, doPlot=TRUE)
dev.off()
head(res)

write.csv(res[order(res$prediction, res$p.value),], paste0(msheetPath, test_text, "res.csv"))






