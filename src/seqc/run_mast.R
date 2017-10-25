suppressMessages(library(MAST))
suppressPackageStartupMessages({library(data.table)})
options(mc.cores = 1) # gives me error messages when I use > 1

loadData <- function(input_data) {
  df <- (read.csv(input_data, row.names=NULL))
}

extractConditions <- function(df) {
  # extract conditions (sg) from column names of the df
  sg <- factor(unlist(df[1]))
  return(sg)
}

annotateDF <- function(df, sg) {
  df[1] <- NULL
  df <- t(df)  
  names(df) <- sg
  return(df)
}

runMAST <- function(df, sg) {
  # extract columns and row information
  # add a cell number column to avoid duplicate row names
  wellKey <- seq_len(dim(df)[2])
  wellKey <- lapply(wellKey, toString)
  condition <- as.numeric(unlist(as.list(sg)))
  cdata <- data.frame(cbind(wellKey=wellKey, condition=condition))
  fdata <- data.frame(primerid=row.names(df))
  
  # create the sca object. Note that we do filtering before
  # we create the test matrix, so no additional filtering of cells is added here
  exprsArray <- as.matrix(df)
  dimnames(exprsArray)[[2]] <- cdata$wellKey
  sca <- FromMatrix(exprsArray, cdata, fdata)
  
  # calculate cellular detection rate
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  colData(sca)$cond <- as.numeric(unlist(as.list(sg)))
  
  # carry out DE analysis
  zlmCond <- zlm.SingleCellAssay(~cond + cngeneson, sca)
  #res <- lrTest(zlmCond, CoefficientHypothesis("cond"))

  #only test the cluster coefficient.
  summaryCond <- summary(zlmCond, doLRT=TRUE)
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='cond' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='cond' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
  
  fcHurdle <- fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- fcHurdle[(fdr<=0.05) & (abs(coef)>=log2(1.25)) ]
  setorder(fcHurdleSig, fdr)
  
  return(fcHurdleSig)
}

saveResult <- function(result, filename) {
  resultDf <- as.data.frame(result)
  colnames(resultDf)[1] = 'gene'
  colnames(resultDf)[2] = 'p'
  colnames(resultDf)[3] = 'logFC'
  colnames(resultDf)[6] = 'p.fdr.adj'
  resultDf <- resultDf[,c('gene','p','p.fdr.adj','logFC')]
  write.table(resultDf, file = filename, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
}

testMAST <- function(input_filename, save_filename) {
  df <- loadData(input_filename)
  sg <- extractConditions(df)
  df <- annotateDF(df, sg)
  result <- runMAST(df, sg)
  saveResult(result, save_filename)
}

# args should be:
# 1. input_filename
# 2. output_filename

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

testMAST(args[1], args[2])
