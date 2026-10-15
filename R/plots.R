#' @title ascat.plotSegmentedData.RRBS
#' @description Plot segmentated BAF LogR 
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#' 
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#' @author Peter Van Loo
#' @noRd
ascat.plotSegmentedData.RRBS <- function (ASCATobj, lim_logR=2) 
{
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    Select_nonNAs = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) = rownames(ASCATobj$Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    png(filename = paste(ASCATobj$samples[arraynr], ".ASPCF.png", 
                         sep = ""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, 
        cex.main = 3, cex.axis = 2)
    r = ASCATobj$Tumor_LogR_segmented[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), 
                                      arraynr]
    beta = ASCATobj$Tumor_BAF_segmented[[arraynr]][, , drop = FALSE]
    plot(c(1, length(r)), c(-lim_logR ,lim_logR), type = "n", xaxt = "n", 
         main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], 
                      ", LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), 
                               arraynr], col = "red", pch = 10, cex = 0.20)
    points(r, col = "blue")
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len + chrk_tot_len_prev)/2
      text(tpos, lim_logR-0.5, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    plot(c(1, length(beta)), c(0, 1), type = "n", xaxt = "n", 
         main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], 
                      ", BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), 
                              arraynr], col = "red", pch = 10, cex = 0.20)
    points(beta, col = "blue")
    points(1 - beta, col = "blue")
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len + chrk_tot_len_prev)/2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    dev.off()
  }
}

#' @title ascat.plotRawData
#' @description Plot BAF LogR
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#' @param pch type of data points in plot
#' @param cex size of data points in plot
#' @param lim_logR y-axis limits on logR plot
#' 
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#' @author Peter Van Loo
#' @keywords internal
ascat.plotRawData.flags = function(ASCATobj, pch, cex, lim_logR) {
  return(1)
  print.noquote("Plotting tumor data")
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    png(filename = paste(ASCATobj$samples[i],".tumour.png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
    plot(c(1,dim(ASCATobj$Tumor_LogR)[1]), c(-lim_logR ,lim_logR ), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[,i],col="red")
    #points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,2,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    
    plot(c(1,dim(ASCATobj$Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[,i],col="red", pch = pch, cex = cex)
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }
  
  if(!is.null(ASCATobj$Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(ASCATobj$Germline_LogR)[2]) {
      png(filename = paste(ASCATobj$samples[i],".germline.png",sep=""), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
      plot(c(1,dim(ASCATobj$Germline_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,2,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(ASCATobj$Germline_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_BAF[,i],col="red", pch = pch, cex = cex)
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
}
