\name{getPRECISE}
\alias{getPRECISE}
\title{
getPRECISE
}
\description{
This function is to calculate PRECISE network and score for all patients
}
\usage{
getPRECISE(outlist, pdlist, nodes, delta)
}
\arguments{
  \item{outlist}{
the outlast object from getBMS
}
  \item{pdlist}{
the pdlist object from getBMS
}
  \item{nodes}{
node names}
  \item{delta}{
a cutoff for patient-specific activation status for a protein
}
}
\details{
To determine the activation status of a protein, we compute the posterior probabilities of the protein to lie in the delta-interval around zero (neutral), to be greater than delta, (activated) or less than -delta (suppressed).
}
\value{
 \item{net}{n x p matrix for node statuses, activated(1), suppressed(-1), or neutral(0)}
 \item{net.status}{length n vector for network statuses, activated (1), suppressed (-1) or neutral (0)}
\item{score.mat}{n x 3 matrix for network scores for activated (first column), suppressed (second column), and neutral (third column)}
\item{samplename}{length n vector for sample names. all other objects are in this order}
}
