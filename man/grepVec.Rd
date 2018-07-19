\name{grepVec}
\alias{grepVec}
\title{
grepVec
}
\description{
grep a vector x from a vector y
}
\usage{
grepVec(x, y, exact = F, Union = T, ignore.case = F)
}
\arguments{
  \item{x}{
a vector
}
  \item{y}{
a vector
}
  \item{exact}{
a logical value. If true, find exact matches. Otherwise grep the pattern
}
  \item{Union}{
a logical value. If true, find union of matches for y. Otherwise, find matches for each element of y.
}
  \item{ignore.case}{
When exact=F and ignore.case=T, case is ignored. Otherwise, case sensitive matching.
}
}
\value{
position of x in y
}
