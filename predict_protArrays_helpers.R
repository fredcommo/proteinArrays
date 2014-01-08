# Predict - proteinArrays Rcode


#################################
# Helper functions
readFile <- function(File = file.path(dir, file), sampleNames=NULL){
  l <- readLines(File, n=50)
  rawData <- read.csv(File, skip = grep("Block", l)-1, sep = '\t')
  if(!is.null(sampleNames))
    for(i in unique(rawData$Block))
      rawData$Block[rawData$Block == i] <- sampleNames[i]
  return(rawData)
}
getData <- function(File = file.path(dir, file), sampleNames=NULL){
  rawData <- readFile(File, sampleNames)
  Blocks <- split(rawData, rawData$Block)
  protNames <- as.character(Blocks[[1]]$Name[-grep('POS|NEG', Blocks[[1]]$Name)])
  
  Data <- lapply(1:length(Blocks), function(i){
    current <- Blocks[[i]]
    blank <- mean(current[,grep("^SNR", colnames(current))][grep("NEG", current$Name)], na.rm = TRUE)
    values <- current[,grep("^SNR", colnames(current))][-grep('POS|NEG', current$Name)]#-blank
#     blank <- mean(current$SNR.635[grep("NEG", current$Name)], na.rm = TRUE)
#     values <- current$SNR.635[-grep('POS|NEG', current$Name)]#-blank
#     blank <- mean(current$F635.Median...B635[grep("NEG", current$Name)], na.rm = TRUE)
#     values <- current$F635.Median...B635[-grep('POS|NEG', current$Name)]-blank
    values[values <= 0] <- 1 #+ abs(rnorm(sum(values<=0), 0, 1e-12))
    values
    })
  Data <- do.call(cbind, Data)
  Data <- .normData(Data)
  colnames(Data) <- names(Blocks)
  Data <- cbind.data.frame(protein = protNames, Data)
  return(Data)
}
.normData <- function(Data){
  qData <- normalize.quantiles(Data)
  qData[qData<=0] <- 1
  return(log2(qData))
}
zeroCorrection <- function(x){
  Sx <- sapply(seq(1, length(x), by = 3), function(i) sd(x[i:(i+2)]))
  Sx <- mean(Sx)
  if(any(x == 0))
    x[x == 0] <- abs(rnorm(sum(x==0), 0, Sx))
  return(x)
}
computeLRatios <- function(Data, Ref, Sample, B){
  prots <- unique(as.character(Data$protein))
  X <- zeroCorrection(Data[ ,colnames(Data) == Ref])
  Y <- zeroCorrection(Data[ ,colnames(Data) == Sample])
  
  output <- lapply(prots, function(prot){
    cat(prot, '\n')
    x <- X[which(Data$protein == prot)]
    y <- Y[which(Data$protein == prot)]
    statTest <- Resamp.T.v4(x, y, B)
    l2r <- statTest$Meany - statTest$Meanx
    cbind(protein = prot, statTest[1:6], Log2Ratio = l2r, statTest[7:8])
  })
  output <- as.data.frame(do.call(rbind, output))
  colnames(output)[grep("Meanx", colnames(output))] <- sprintf("estim_%s", Ref)
  colnames(output)[grep("Meany", colnames(output))] <- sprintf("estim_%s", Sample)
  r <- output[,grep(sprintf("estim_%s$", Ref), colnames(output))]
  s <- output[,grep(sprintf("estim_%s$", Sample), colnames(output))]
  padj <- p.adjust(output$pValue, method = "BH")
  return(cbind.data.frame(output, adjpValue = padj, M = s-r, A = 1/2*(s+r)))
}
plotSignal <- function(rawData, idx, today=format(Sys.Date(), format="%Y-%m-%d")){
  if(!exists("op")) op <- par(no.readonly = TRUE)
  par(mar = c(9, 6, 4, 2), cex.main = 2, cex.sub = 1, cex.lab = 1.5, bty = "n")
  
  prots <- as.factor(rawData$protein)
  sampleNames <- colnames(rawData)[idx]
  col1 <- rgb(0,0,1,.2)
  col2 <- rgb(1,0,0,.2)
  plot(prots, rawData[,idx[1]], ylim = range(rawData[,-1])*1.2, las = 3, type="n",
       xlab = "", ylab = "", main = sprintf("%s Vs. %s", sampleNames[1], sampleNames[2]))
  abline(v=1:nlevels(prots), lty=3, col="grey65")
  plot(prots, rawData[,idx[1]], add=TRUE, col = col1, axes=FALSE)
  plot(prots, rawData[,idx[2]], add=TRUE, col = col2, axes=FALSE)
  sub <- paste('Array_Raybiotech Human RTK Phosphorilation Antibody Array G-series 1_(Ref:AAH-PRTK-G1) -', today)
  title(xlab="Proteins", mgp = c(6, 1, 0))
  title(ylab=expression(Log[2](signal)), mgp = c(4, 1, 0))
  title(sub=sub, mgp = c(6.5, 1, 0))
  legend("topleft", legend = sampleNames, fill = c(col1, col2), cex = 1.25, bty = "n")
  par(op)
}
plotLRatios <- function(ratio, useP = c("p", "adjp"), thresh = .05, today=format(Sys.Date(), format="%Y-%m-%d"),...){
  if(!exists("op")) op <- par(no.readonly = TRUE)
  signifCol = 'steelblue3'; loCol = "red3"; nullCol = 'grey'
  useP <- match.arg(useP)
  switch(useP,
         p = {pIdx <- grep("^pValue$", colnames(ratio))
              leg <- substitute(paste("p<=", thresh), list(thresh=thresh))},
         adjp = {pIdx <- grep("^adjpValue$", colnames(ratio))
                 leg <- substitute(paste("adjp<=", thresh), list(thresh=thresh))}
         )
  p <- ratio[,pIdx]
  x <- ratio$Log2Ratio
  Cols <- ifelse(p <= thresh, signifCol, nullCol)
  
  par(mar = c(9, 8, 4, 2), cex.main = 2, cex.sub = 1, cex.lab = 1.5)
  yAxis <- seq(floor(min(x))-1, ceiling(max(x))+1, by = 2)
  xAxis <- seq(0.75, nrow(ratio)*1.1925, len = nrow(ratio))
  barplot(x, axes = FALSE, col = Cols, border = 'grey30',
          ylim = range(yAxis), xlab = "", ylab = "",...)
  axis(side = 1, at = xAxis, labels = ratio$protein, las = 3, cex.axis = 1.25)
  axis(side = 2, at = yAxis, labels = yAxis, las = 1, cex.axis = 1.5)
  sub <- paste('Array_Raybiotech Human RTK Phosphorilation Antibody Array G-series 1_(Ref:AAH-PRTK-G1) -', today)
  title(xlab="Proteins", mgp = c(6, 1, 0))
  title(ylab="Log2Ratio", mgp = c(4, 1, 0))
  title(sub=sub, mgp = c(6.5, 1, 0))
  yMin <- ifelse(x > 0, 0, x)
  segments(x0 = xAxis, y0 = rep(min(yAxis)*2, length(xAxis)), y1 = yMin, col = 'grey80', lty = 2)
  legend('topleft', legend = leg, cex = 1.5, fill = signifCol, bty = 'n')
  par(op)
}
vPlot <- function(ratio, useP = "adjp", thresh=.05, main = paste(Sample, 'Vs.', Ref),...){
  if(!exists("op")) op <- par(no.readonly = TRUE)
  par(mar = c(6, 5, 4, 2), cex.main = 2, cex.sub = 1, cex.lab = 1.5, cex.axis = 1.25, las = 1)
  useP <- match.arg(useP)
  signifCol = 'steelblue3'; loCol = "red3"; nullCol = 'grey'
  switch(useP,
         p = {pIdx <- grep("^pValue$", colnames(ratio))
              leg <- substitute(paste("p<=", thresh), list(thresh=thresh))},
         adjp = {pIdx <- grep("^adjpValue$", colnames(ratio))
                 leg <- substitute(paste("adjp<=", thresh), list(thresh=thresh))}
  )
  x <- ratio$Log2Ratio
  y <- -log10(ratio[,pIdx])
  Cols <- ifelse(y >= -log10(thresh), signifCol, nullCol)
  plot(x, y, type = "n",
       xlab = expression(Log[2](Ratio)), ylab = expression(-Log[10](p-value)), main = main,
       sub = paste('Array_Raybiotech Human RTK Phosphorilation Antibody Array G-series 1_(Ref:AAH-PRTK-G1) -', today), ...)
  abline(v = 0, lty = 2, lwd = 2)
  points(x, y, pch = 19, col = Cols, cex = 1.5,...)
  legend("bottomleft", legend = leg, pch = 19, col = signifCol, cex = 1.25)
  par(op)
}
HeatMap <- function(X, Factors = NULL, Scale = "none", ...){
  rowCols = rep("lightgrey", nrow(X))
  if(!is.null(Factors)){
    Factors <- as.factor(Factors)
#    rowIds <- Factors
    rowCols <- rgb(seq(.2, .8, len=nlevels(Factors)),
                 seq(.2, .8, len=nlevels(Factors)),
                 seq(.2, .8, len=nlevels(Factors)),
                 .25)
    rowCols <- rowCols[Factors]
    }
  rowIds <- rownames(X)
  heatmap.4(as.matrix(X),
            scale = Scale,
            col = colorpanel(100, "darkblue", "grey95", "darkorange1"),
            RowSideColors = rowCols,
            sepcolor = "grey75", sepwidth = c(0.01, 0.01),
            colsep = 1:ncol(X), rowsep = 1:nrow(X),
            cexCol = 1.25, cexRow = 2, labRow = rowIds, labRowPos = 2,
          density.info = "none", trace = "none", keysize = 1)
}
plotCor <- function(X, Ids, subIds = NULL,...){
  if(!exists("op")) op <- par(no.readonly = TRUE)
  par(mar = c(8, 5, 4, 8), cex.main = 2, cex.sub = 1, cex.lab = 1.5, cex.axis = 1.25, las = 1)
  Ids <- as.factor(Ids)
  rowCols <- rgb(seq(.2, .8, nlevels(Ids)),
                 seq(.2, .8, nlevels(Ids)),
                 seq(.2, .8, nlevels(Ids)),
                 .25)
  rowIds <- colIds <- Ids
  if(!is.null(subIds))
    rowIds <- colIds <- sprintf("%s\n%s", Ids, subIds)
  else rowIds <- colIds <- colnames(X)
#   rownames(X) <- gsub("_Log2Ratio", "", rownames(X))
#   colnames(X) <- gsub("_Log2Ratio", "", colnames(X))
  heatmap.4(as.matrix(X), #Rowv = FALSE, Colv = FALSE, dendrogram = "none",
            scale = "none",
            col = colorpanel(100, "darkblue", "grey95", "darkorange1"),
            # RowSideColors = rowCols[Ids], ColSideColors = rowCols[Ids],
            sepcolor = "grey75", sepwidth = c(0.01, 0.01),
            colsep = 1:ncol(X), rowsep = 1:nrow(X),
            labCol = colIds, cexCol = 2 - log10(ncol(X)/4),
            labRow = rowIds, cexRow = 2 - log10(nrow(X)/4), 
#            labRowPos = 2, labColPos = 3,
            density.info = "none", trace = "none", keysize = .75)
  par(op)
}
.reverse <- function(X, by = c("row", "col")){
  by <- match.arg(by)
  switch(by,
         row = {new <- lapply(nrow(X):1, function(i) X[i,])
                new <- do.call(rbind, new)
                rownames(new) <- rownames(X)[nrow(X):1]},
         col = {new <- lapply(ncol(X):1, function(i) X[,i])
                new <- do.call(cbind, new)
                colnames(new) <- colnames(X)[ncol(X):1]})
  return(new)
}
Resamp.T.v4 <- function(x, y, B = 10000, trimmed = 0){
  
  # compute bootstrap t-values
  
  na.replace <- function(x, trimmed){
    qx <- quantile(x, probs=c(trimmed, 1-trimmed), na.rm = T)
    x[which(x<qx[1] | x>qx[2])] <- NA
    return(x)
  }
  stat.test <- function(x, y){as.numeric(t.test(y, x)$statistic)}
  
  if(any(is.na(x)) | any(is.na(y))){
    NAs <- unique(c(which(is.na(x)), which(is.na(y))))
    x <- x[-NAs]
    y <- y[,-NAs]
    cat(length(NAs), "case(s) suppressed for missingness\n")
  }
  t <- bootX <- bootY <- S2X <- S2Y <- t.boot <- c()
  nx <- length(x)
  ny <- length(y)
  Boot <- lapply(1:B, function(b){
    if(b%%(B/10) == 0) cat("B:", b, "\t")
    xb <- sample(x, nx, replace = TRUE)# + rnorm(nx, 0, sd(x))
    yb <- sample(y, ny, replace = TRUE)# + rnorm(ny, 0, sd(y))
    if(trimmed>0){
      xb <- apply(xb, 1, na.replace, trimmed = trimmed)
      xb <- t(xb)
      yb <- apply(yb, 1, na.replace, trimmed = trimmed)
      yb <- t(yb)
    }
    t <- try(stat.test(xb, yb), silent = TRUE)
    if(class(t)[[1]]=="try-error") t = NA
    c("Meanx" = mean(xb), "S2x" = var(xb), "Meany" = mean(yb), "S2y" = var(yb), "stat" = t)
  })
  cat("\n")
  Boot <- as.data.frame(do.call(rbind, Boot))
  out <- as.data.frame(t(apply(Boot, 2, mean, na.rm = TRUE)))
  p <- pt(abs(out$stat), df = (nx+ny-2), lower.tail = FALSE)*2
  Xbias <- out$Meanx - mean(x, na.rm = TRUE)
  Ybias <- out$Meany - mean(y, na.rm = TRUE)
  return(cbind(out[1:2], "xbias" = Xbias, out[3:4], "ybias" = Ybias, out[5], pValue = p))
}

# End helper functions
#################################
