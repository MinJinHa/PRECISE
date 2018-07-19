\name{grepVecAnno}
\alias{grepVecAnno}
\title{
grepVecAnno
}
\description{
grep a vector x from a vector y using annotation file
}
\usage{
grepVecAnno(x, y, anno)
}
\arguments{
  \item{x}{
a vector
}
  \item{y}{
a vector
}
  \item{anno}{
an annotation matrix which links x and y. The first and second column includes elements for x and y, respectively.
}
}
\value{
 a list which includes positions of x in y
}
\note{
anno file is not necessarily one to one, multiple positions are possible for an element of x}
