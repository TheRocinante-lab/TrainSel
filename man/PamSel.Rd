\name{PamSel}
\alias{PamSel}
\title{
PamSel
}
\description{
Selection using partitioning around medoids.
}
\usage{
PamSel(M=NULL,
distancemethod=c("euclidean",  "maximum", "manhattan", "canberra", "binary", "minkowski"),
ntoselect=NULL)
}
\arguments{
  \item{M}{Data matrix.}

\item{distancemethod}{one of "euclidean",  "maximum", "manhattan",
"canberra", "binary", "minkowski".}
\item{ntoselect}{number of samples to select.
}
}
\details{
Selection by partitioning around medoids.
}
\value{vector of integers that correspond to the selected rows of M.}

\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez, Javier Fernandez-Gonzalez
}

\examples{
\dontrun{
data(WheatData)
solpam<-PamSel(M=Wheat.M, distancemethod="euclidean", ntoselect=20)
solpam
}

}
