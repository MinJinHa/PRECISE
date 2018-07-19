\name{getPosteriors}
\alias{getPosteriors}
\title{
getPosteriors
}
\description{
This function calculates posterior graph using median probability model from the results of getBMS
}
\usage{
getPosteriors(outlist, nodes)
}
\arguments{
  \item{outlist}{
object from getBMS
}
  \item{nodes}{
node names
}
}

\value{
\item{G}{p x p adjacency matrix for protein posterior network}
\item{intGlist}{a p length list for edges from upstream data to each of the proteins}
}
\note{
There are three different node types for upstream platforms, miRNA, ME (gene expression explained by methylation) and NME (gene expression not explained by methylation)
}
