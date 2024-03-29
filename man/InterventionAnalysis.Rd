% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMTEssence.R
\name{InterventionAnalysis}
\alias{InterventionAnalysis}
\alias{intervention}
\title{my title}
\usage{
intervention(data, metric, pairs, pre_period, post_period, label = NULL,
  plots = F)
}
\arguments{
\item{data}{a data frame}

\item{metric}{conversion metric}

\item{pairs}{a list of market pairs}

\item{pre_period}{vector of two dates representing start and end of pre period}

\item{post_period}{vector of two dates representing start and end of post period}

\item{label}{label to match file path used for finding pairs}

\item{plots}{Boolean for plot inclusion}
}
\value{
list of valid tests
}
\description{
my title
}
