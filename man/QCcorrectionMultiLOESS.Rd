% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QCcorrectionMultiLOESS.R
\name{QCcorrectionMultiLOESS}
\alias{QCcorrectionMultiLOESS}
\title{Generate values for metabolites normalization}
\usage{
QCcorrectionMultiLOESS(LCdata)
}
\arguments{
\item{LCdata}{Matrix of data obtained (mainly by LC-MS) that included four data columns ("Compound Name","Order","QC","Day") and then one coulm for each compound or entity detected.}
}
\value{
A data set similar to LCdata matrix but with duplicated columns for each compound or entity with the area needed to normalize each of them.
}
\description{
According to the area of QC along each day, this function generates values for each sample injected along the day that are going to be used for data normalization.
}
\examples{
\dontrun{
correctedLCdata<-QCcorrectionMultiLOESS(LCdata)
}
}

