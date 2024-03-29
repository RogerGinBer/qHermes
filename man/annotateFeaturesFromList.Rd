% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateFeatures-family.R
\name{annotateFeaturesFromList}
\alias{annotateFeaturesFromList}
\title{annotateFeaturesFromList}
\usage{
annotateFeaturesFromList(
  XCMSnExp,
  RHermesExp,
  SOI_id,
  MS2Exp = NA,
  RTtol = 10,
  filter = TRUE,
  quantifyMissingSOI = FALSE
)
}
\arguments{
\item{XCMSnExp}{XCMS object that contains a list of features (ie. after
running \code{groupChromPeaks()}).}

\item{RHermesExp}{RHermes object that contains a SOI (Scans of Interest)
list.}

\item{SOI_id}{Index of the SOI that will be used as reference for the matching.}

\item{MS2Exp}{Numeric. If provided, index of the MS2Exp object inside the
RHermesExp object from which the MS2 scans will be imported and added to
the feature table.}

\item{RTtol}{Numeric. Retention time shift tolerance in seconds.}

\item{filter}{Logical. Whether to remove those features that can't be matched
to a SOI.}

\item{quantifyMissingSOI}{Logical. Whether to quantify via EIC (Extracted Ion
Chromatogram) those SOIs that were not matched with a feature.}
}
\description{
Finds a correspondance between XCMS features and RHermes SOIs
in order to filter the former and obtain a cleaner peak matrix
}
