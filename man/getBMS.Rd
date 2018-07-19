\name{getBMS}
\alias{getBMS}
\title{
getBMS
}
\description{
Given data and prior information, this function samples all
     possible model combinations via MC3 or enumeration and returns
     aggregate results. This procedure is performed for all nodes}
\usage{
getBMS(dat, Gmat, predictiveDensity = TRUE)
}
\arguments{
  \item{dat}{
a list of data that includes response vectors and design matrices for all nodes
}
  \item{Gmat}{
an adjacency matrix for prior causal protein structure
}
  \item{predictiveDensity}{
a logical value indicating whether predictive density is computed or not
}
}
\value{
\item{outlast}{an output list from bms functions for all nodes}
\item{pdlist}{an output list from predictive.density function for all nodes }
}
