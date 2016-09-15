#' Generate values for metabolites normalization
#' 
#' According to the area of QC along each day, this function generates values for each sample injected along the day that are going to be used for data normalization.
#' @param LCdata Matrix of data obtained (mainly by LC-MS) that included four data columns ("Compound Name","Order","QC","Day") and then one coulm for each compound or entity detected.
#' @examples
#' \dontrun{
#' correctedLCdata<-QCcorrectionMultiLOESS(LCdata)
#' }
#' @export
#' @import plyr
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats coefficients
#' @importFrom stats poly
#' @importFrom stats loess
#' @importFrom stats predict
#' @return A data set similar to LCdata matrix but with duplicated columns for each compound or entity with the area needed to normalize each of them.

QCcorrectionMultiLOESS<-function(LCdata) {
  pos=1;
  envir=as.environment(pos);
  LCdata<-LCdata
  LCdataframe<-as.data.frame(LCdata)
  LCdataQC<-subset(LCdata,LCdata$QC==1)
  nDAY<-ddply(LCdataframe,.("QC","Day"),summarise,nSamples=sum(!is.na(LCdataframe$QC)))
  nDAYQC<-nDAY[nDAY$QC==1,]
  nDAYQC<-na.omit(nDAYQC)
  Day1<-list(nDAYQC$Day)
  Day2<-as.data.frame(Day1)
  Day<-dim(Day2)[1]
  a=1;
  while (a<(Day+1)){
    assign(paste0('Day_',Day2[a,]),subset(LCdata,(LCdata$Day==(Day2[a,]))),envir=envir);
    assign(paste0('QCsDay_',Day2[a,]),subset(LCdataQC,(LCdataQC$Day==(Day2[a,]))),envir=envir);
    a<-a+1
  }
  matrix.names<-paste0('Day_',Day2[1,]:Day2[Day,])
  list<-lapply(matrix.names,get, envir=envir)
  matrixQC.names<-paste0('QCsDay_',Day2[1,]:Day2[Day,])
  listQC<-lapply(matrixQC.names,get)
  PreCompoundsList<-colnames(LCdata)
  InfoColumns<-c("Sample","Order","QC","Day")
  CompoundsList<-setdiff(PreCompoundsList,InfoColumns)
  n=1;
  while(n<(Day+1)){
    bb<-as.data.frame(listQC[n],envir=envir)
    cc<-as.data.frame(list[n],envir=envir)
    m=1;
    while(m<(length(CompoundsList)+1)){
      lmloess<-loess(bb[,CompoundsList[m]]~(bb[,"Order"]),span=0.3,degree=2);
      cc[,paste0(CompoundsList[m],'QC')]<-predict(lmloess,LCdata[,"Order"]);
      m<-m+1
    }
    if(exists("PreNormLCdata")){PreNormLCdata<-rbind(PreNormLCdata,cc)
    } else {PreNormLCdata<-cc}
    n<-n+1
  }
    
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