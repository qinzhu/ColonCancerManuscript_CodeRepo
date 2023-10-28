




library(VisCello)
library(igraph)
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

eset <- readRDS("../hcc_final_cello/eset.rds")
clist <- readRDS("../hcc_final_cello/clist.rds")


test_text <- "cellcycle_epi_human_"

dep_vis <- clist$`Epithelial cells [Final, wt OGN, IFF]`
eset_epi <- eset[, dep_vis@idx]

cc_res <-as.data.frame.matrix(table(eset_epi$Dataset, eset_epi$Cell_cycle_phase))
write.xlsx(cc_res, paste0(msheetPath,test_text,"dataset_barplot.xlsx"), rowNames = T)
cc_res <- cc_res[!grepl("LiverMet", rownames(cc_res)),]
cc_res <- cc_res/rowSums(cc_res)
cc_res$dataset <- rownames(cc_res)
cc_res <- reshape2::melt(cc_res)
cc_color <- c("G1" = "#66bd63", "S" = "#fee08b", "G2M" = "#d73027")
cc_res$variable <- factor(cc_res$variable, levels = names(cc_color))
g1<-ggplot(cc_res, aes(x = dataset, y = value, fill = variable)) + 
    geom_bar(stat="identity") + 
    theme_classic() + 
    scale_fill_manual(values = cc_color) + 
    theme(legend.position="top", 
          axis.text.x=element_text(angle = 90,hjust = 1),
          text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))+
    xlab("Cell type") +
    ylab("Fraction")

ggsave(paste0(mplotPath, test_text, "dataset_barplot.pdf"), g1, width = 3.5, height=3, units = "in", device = "pdf")







