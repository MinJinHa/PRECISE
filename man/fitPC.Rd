\name{fitPC}
\alias{fitPC}
\title{fitPC}
\description{
Statbility selection for CPDAG using PC algorithm
}
\usage{
fitPC(dat, alpha, stable = FALSE, alpha.array, B, labels, verbose = F)
}
\arguments{
  \item{dat}{
a n x p data matrix where n is the number of samples and p is the number of nodes with colnames as the node names 
}
  \item{alpha}{
a scalar for the tuning parameter for PC algorithm
}
  \item{stable}{
a logical value. If TRUE, selection probabilities are calculated from subsampling 
}
  \item{alpha.array}{
a vector for the set of alpha values for the stability selection
}
  \item{B}{
a scalar for the number of subsampling in the stability selection
}
  \item{labels}{
p x 1 vector for node names
}
  \item{verbose}{
a logical value for details
}
}
\value{
 \item{fit}{pcAlgo object (pcalg package) for CPDAG}
  \item{pcA}{an no. edges x length(alpha.array) x B dimensional array that includes the all edges across all tuning parameters in alpha.array and B subsamples. 1 indicates -> (directed edge from column 1 to column 2 in the rownames object), 2 indicates <- (directed edge from column 2 to column 1 in the rownames object), and 3 indicates - (undirected edge between column 1 and column 2). }
\item{pcSel}{3 x no.edges x length(alpha.array) list that includes 3 types of no.edges x length(alpha.array) matrices of selection probabilities for each of the directions (->, <-, -).}
\item{pcmaxSel}{no.edges x 3 matrix for maximum selection probabilities for each of the directions, -> (first column), <- (second column), - (third column)}
\item{rownames}{no. edges x 2 matrix for pairs of nodes (first node in the first column and second node in the second column) for all possible edges. All the outputs are in the order of the rows of this matrix. 
}
}