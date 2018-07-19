\name{getPPD}
\alias{getPPD}
\title{
getPPD
}
\description{
This function is to calculate posterior predictive distribution (delta interval) of each of the proteins for each patient
}
\usage{
getPPD(pd, delta)
}
\arguments{
  \item{pd}{
the pdlist output from getBMS
}
  \item{delta}{
a cutoff for patient-specific activation status for a protein}
}
\details{
To determine the activation status of a protein, we compute the posterior probabilities of the protein to lie in the delta-interval around zero (neutral), to be greater than delta, (activated) or less than -delta (suppressed).
}
\value{
 \item{results}{a matrix includes all the patient-specific results, samplename, fitted value, sd, positive prob. (activated), negative prob. (suppressed), neutral prob., and activation status }
\item{pd}{an object from pred.density }
}
