#' Generate values for metabolites normalization
#' 
#' According to the area of QC along each day, this function generates values for each sample injected along the day that are going to be used for data normalization.
#' @param LCdata Matrix of data obtained (mainly by LC-MS) that included four data columns ("Compound Name","Order","QC","Day") and then one coulm for each compound or entity detected.
#' @examples
#' \dontrun{
#' correctedLCdata<-QCcorrectionLOESS(LCdata)
#' }
#' @export
#' @import plyr
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats coefficients
#' @importFrom stats poly
#' @return A data set similar to LCdata matrix but with duplicated columns for each compound or entity with the area needed to normalize each of them.

QCcorrectionLOESS<-function(LCdata) {
  globalenv()
  LCdata<-LCdata
  LCdataframe<-as.data.frame(LCdata)
  LCdataQC<-subset(LCdata,LCdata$QC==1)
  PreCompoundsList<-colnames(LCdata)
  InfoColumns<-c("Sample","Order","QC","Day")
  CompoundsList<-setdiff(PreCompoundsList,InfoColumns)
  LCdatanew<-LCdata

  m=1;
    while(m<(length(CompoundsList)+1)){
      lmloess<-loess(LCdata[,CompoundsList[m]]~(LCdata[,"Order"]),span=0.3,degree=2);
      LCdatanew[,paste0(CompoundsList[m],'QC')]<-predict(lmloess,LCdata[,"Order"]);
      m<-m+1
    }
  PreNormLCdata<-LCdatanew

  p=1;
  while(p<(length(CompoundsList)+1)){
    PreNormLCdata[,paste0(CompoundsList[p],'_Norm')]<-PreNormLCdata[,(CompoundsList[p])]/PreNormLCdata[,paste0(CompoundsList[p],'QC')]
    p<-p+1
  }
    NormLCdata1<-PreNormLCdata[,grepl("Norm",names(PreNormLCdata))]
    NormLCdata2<-PreNormLCdata[,InfoColumns]
    NormLCdata<-cbind(NormLCdata2,NormLCdata1)
  return(NormLCdata)
}