library(ggplot2)
plotPCA_xiao = function(x,condition, ntop=500, returnData=FALSE)
{
  ## x is counts matrix
  #specific the shape here: http://www.sthda.com/english/wiki/ggplot2-point-shapes
  #0: square, 1 circle, 2 tringle, 15, 16,17 is solid
  # shapes <- c()
  
  library(genefilter)
  # calculate the variance for each gene
  rv <- rowVars(x)
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(x[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  colData <- data.frame(row.names=colnames(x),condition=condition)
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],condition=condition,names=colnames(x))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  p=ggplot(data=d, aes(x=PC1, y=PC2)) + geom_point(size=2,aes(color=condition),shape=21) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +theme_bw()
    #scale_shape_manual(values= c(letters, 0:9)) +
    #scale_color_manual(values=c("#203864","#2f5597","#8faadc","#b4c7e7","#2f5597"))
 #
    #geom_text(data=d,aes(PC1,PC2,label=names))
  return(p)
}

