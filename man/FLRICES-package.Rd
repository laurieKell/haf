\name{FLRICES-package}
\alias{FLRICES-package}
\alias{FLRICES}
\docType{package}
\title{
FLR-based Implementation of ICES Harvest Control Rules
}
\description{
A comprehensive R package implementing ICES (International Council for the Exploration of the Sea) harvest control rules within the FLR (Fisheries Library for R) framework. The package provides standardized implementations of various ICES advice rules including hockey-stick HCRs, RFB rules for data-limited stocks, and Nephrops-specific rules.
}
\details{
\tabular{ll}{
Package: \tab FLRICES\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2024-12-19\cr
License: \tab GPL-3\cr
}

The package includes the following main functions:

\strong{Hockey-Stick HCR (Category 1 & 2 stocks):}
\itemize{
  \item \code{\link{hcrICESAR}} - Standard ICES hockey-stick harvest control rule
}

\strong{RFB Rule (Category 3 stocks):}
\itemize{
  \item \code{\link{rfbParams}} - Set up RFB control parameters
  \item \code{\link{rfbRule}} - Apply RFB harvest control rule
  \item \code{\link{simulateRFB}} - Full RFB simulation example
}

\strong{Survey-Based HCR:}
\itemize{
  \item \code{\link{hcrSurvey}} - Generic survey-based HCR for stocks managed through survey indices
  \item \code{\link{hcrNephrops}} - Nephrops-specific HCR (wrapper around hcrSurvey)
}

\strong{Visualization:}
\itemize{
  \item \code{\link{hcrPlot}} - Standard HCR plot
  \item \code{\link{hcrPlot2}} - Enhanced HCR plot with management zones
  \item \code{\link{plotSurveyHCR}} - Survey-based HCR visualization
  \item \code{\link{plotNephropsHCR}} - Nephrops HCR visualization
}
}
\author{
Your Name <your.email@example.com>
}
\references{
ICES (2025). Technical Guidelines for Data-limited Stock Assessments. ICES CM 2025/ACOM:68

ICES (2023). ICES Technical Guidelines for Category 1 and 2 stocks. ICES CM 2023/ACOM:XX
}
\keyword{package} 