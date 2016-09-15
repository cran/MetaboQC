#' Representate the compounds area (normalized or not) as a function of their injection order to study trends.
#' 
#' Export graphs for each compound included in LCdata matrix in which the area of the specified compound is represented vs the injection order.
#' @param LCdata Matrix of data obtained (mainly by LC-MS) that included four data columns ("Compound Name","Order","QC","Day") and then one coulm for each compound or entity detected (normalized or not).
#' @param g Number of compounds for which the graph should be obtained
#' @param NameDataSet A name for the data set that is going to be used for the pdf file name. It must be given in quotes
#' @examples
#' \dontrun{
#' graphQC(LCdata,3,"datasetName")
#' }
#' @export
#' @import graphics
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @return Multiple graphs of the compounds area (normalized or not) vs the injection order.

graphQC<-function(LCdata,g,NameDataSet){
  a=1;
  PreCompoundsList<-colnames(LCdata)
  InfoColumns<-c("Sample","Order","QC","Day")
  CompoundsList<-setdiff(PreCompoundsList,InfoColumns)
  if(g>(length(CompoundsList))){g<-length(CompoundsList)
  } else {g<-g}
  pdf(file=paste0('xy-Plot for ',g,' compounds from ',NameDataSet,'.pdf'),height=8,width=15)
  while (a<(g+1)) {
    plot(LCdata$Order,LCdata[,CompoundsList[a]],xlab="Order of injection",ylab=CompoundsList[a],main=paste0(CompoundsList[a],' vs injection order'),col=ifelse(LCdata$QC==1,"red",ifelse(LCdata$Day%%2==0,"blue","green")),pch=(ifelse(LCdata$Day%%2==0,16,17)))
    a<-a+1
  }
  dev.off()
}
