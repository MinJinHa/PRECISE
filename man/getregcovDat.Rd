\name{getregcovDat}
\alias{getregcovDat}
\title{
getregcovDat
}
\description{
This function is to obtain a RPPA response vectors and covariate matrices from mRNA, miRNA, methylation data, and protein regulations from the estimated graph (G)
}
\usage{
getregcovDat(Gmat, RPPAdat, mRNA = NULL, miRNA = NULL, Methylation = NULL, anno.miRNA = NULL, anno.methyl = NULL, integrative = T)
}
\arguments{
  \item{Gmat}{
 p x p adjacency matrix (row->column) which includes weights for edges
}
  \item{RPPAdat}{
n x p protein expression matrix
}
  \item{mRNA}{
a list including mRNA Data (p_mRNA x n matrix) and Des (p_mRNA x 2 matrix including GeneSymbol and EntrezID)
}
  \item{miRNA}{
a list including miRNA Data (p_miRNA x n matrix) and Des (p_mRNA x 1 matrix including miRNA names)
}
  \item{Methylation}{
a list including Methylation Data (p_methyl x n matrix) and Des (p_methyl x 4 matrix including REF, GeneSymbol, ChromosomeID and Coordinate ID). Only the cg numbers are used.
}
  \item{anno.miRNA}{
annotation matrix for miRNA including gene name in the first column and miRNA names in the second column
}
  \item{anno.methyl}{
annotation matrix for Methylation including gene name in the first column and Methylation names in the second column
}
  \item{integrative}{
TRUE if use miRNA, Methylation, and mRNA
}
}
\value{
 \item{ylist}{a list including scaled protein response vectors for all nodes}
\item{Xlist}{a list including scaled covariate matrices for all nodes}
}
\note{
Note: all samples should be matched. The rows of RPPAdata and columns of mRNA$Data, miRNA$Data, Methylation$Data should be matched.
}
