% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plantGlycoMSFunctions.R
\name{Read.MS2Data}
\alias{Read.MS2Data}
\title{A function to import MS2 binary data generated with msaccess (ProteoWizard).}
\usage{
Read.MS2Data(IDPdb, dir = "convertedData/MS2Data_Chym1_ELUTE")
}
\arguments{
\item{IDPdb}{a data.frame containing peptide spectrum matches}

\item{dir}{name of the directory containing the MS2 binary data}
}
\description{
this function imports MS2 binary data
}
\examples{
Read.MS2Data()
}
\keyword{data}
\keyword{import}
