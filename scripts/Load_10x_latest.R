
library(monocle)
library(Matrix)

raw_to_monocle <- function(dirPath, minUMI = 500, minGenes = 200){
    message("loading matrix")
    file_list<-list.files(dirPath, full.names = T)
    barcode.path <- grep("barcodes.tsv", file_list, value = T)
    features.path <- grep("features.tsv", file_list, value = T)
    matrix.path <- grep("matrix.mtx", file_list, value = T)
    raw_matrix <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(raw_matrix) = barcode.names$V1
    rownames(raw_matrix) = feature.names$V1
    message(dim(raw_matrix)[1])
    message(dim(raw_matrix)[2])

    message("filtering matrix")
    ###Filter 1: ~500 UMI and 200 genes detected
    Total_mRNAs <- Matrix::colSums(raw_matrix)
    num_genes_expressed <- Matrix::colSums(raw_matrix > 0)
    raw_matrix=raw_matrix[,Total_mRNAs>=minUMI & num_genes_expressed>=minGenes]
    message(dim(raw_matrix)[1])
    message(dim(raw_matrix)[2])
 #   points(pData(raw_matrix)$Total_mRNAs,pData(raw_matrix)$num_genes_expressed,log='xy',col=2)
    fmeta <- feature.names[,c(1,2)]
    colnames(fmeta) <- c("gene_id", "gene_short_name")
    rownames(fmeta) <- fmeta$gene_id
    pmeta <- data.frame(barcode = colnames(raw_matrix))
    rownames(pmeta) <- pmeta$barcode
    #Convert to monocle object
    cds = newCellDataSet(Matrix(raw_matrix, sparse=T), 
                         phenoData =     new("AnnotatedDataFrame", data = pmeta),
                         featureData =   new("AnnotatedDataFrame", data = fmeta))
    return(cds)
}



combine_monocle_datasets <- function(dirPaths, minUMI = 500, minGenes = 200){
  for(i in 1:length(dirPaths)){
    dirPath <- dirPaths[i]
    message(dirPath)
    dirName <- names(dirPaths)[i]
    raw_matrix <- raw_to_monocle(dirPath, minUMI, minGenes)
    #   points(pData(raw_matrix)$Total_mRNAs,pData(raw_matrix)$num_genes_expressed,log='xy',col=2)
    if(i==1){
      CombinedMatrix = exprs(raw_matrix)
      CombinedPh =pData(raw_matrix)
      rownames(CombinedPh)<-paste(rownames(CombinedPh),dirName,sep="-")
      colnames(CombinedMatrix)<-rownames(CombinedPh)
      CombinedPh[,"barcode"] =  paste(CombinedPh[,"barcode"],dirName,sep="-")
      CombinedPh$Dataset <- dirName
      CombinedF = fData(raw_matrix)
    } else {
      ThisPh = pData(raw_matrix)
      rownames(ThisPh)<-paste(rownames(ThisPh),dirName,sep="-")
      ThisPh[,"barcode"] =  paste(ThisPh[,"barcode"],dirName,sep="-")
      ThisPh$Dataset <- dirName
      CombinedPh = rbind(CombinedPh,ThisPh)
      ThisMatrix = exprs(raw_matrix)
      colnames(ThisMatrix)<-rownames(ThisPh)
      if(!identical(rownames(ThisMatrix), rownames(CombinedMatrix))){
        warning("Gene names does not match!!")
        shared_g <- intersect(rownames(ThisMatrix), rownames(CombinedMatrix))
        CombinedMatrix = cbind(CombinedMatrix[shared_g,],ThisMatrix[shared_g,])
        CombinedF = CombinedF[shared_g,]
      } else {
        CombinedMatrix = cbind(CombinedMatrix,ThisMatrix)
      }
    }        
  }
  #Convert to monocle object
  CombinedMatrix<-Matrix(CombinedMatrix, sparse = T)
  sum(colnames(CombinedMatrix)!=rownames(CombinedPh))
  cds = newCellDataSet(CombinedMatrix,
                       phenoData =     new("AnnotatedDataFrame", data = CombinedPh),
                       featureData =   new("AnnotatedDataFrame", data = CombinedF))
  return(cds)
}



