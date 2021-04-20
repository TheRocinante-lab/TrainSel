\name{MakeTrainSelData}
\alias{MakeTrainSelData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a data of class 'TrainSel_Data'.
}
\description{
Create a data of class 'TrainSel_Data'. This is a data preprocessor that creates
'TrainSel_Data' for use as input in 'TrainSel'  for the 'data' parameter for the defaults mixed model based selection criteria (CDMEAN).
}
\usage{
MakeTrainSelData(M = NULL,
          K = NULL,
          R=NULL,
          Vk = NULL,
          Ve = NULL,
          lambda = NULL,
          X = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{numerically coded (samples x features) matrix.
}
  \item{K}{positive definite relationship matrix of (samples x samples)
}
  \item{R}{positive definite relationship matrix of the errors of samples. Default value is NULL, in which case R matrix is assumed as an identity matrix.
}
  \item{Vk}{Default is NULL, only used for selecting several sets of samples (Mixed model II): Vk is covariance matrix for different sets of samples.
}
  \item{Ve}{Default is NULL, only used for selecting several sets of samples (Mixed model II): , Ve is the covariance matrix for the errors of sets samples.
}
  \item{lambda}{Defaults is 1. Ratio of Ve to Vk. For mixed model I this is the only variance parameter to specify. For mixed model II specifying it has the effect of multiplying the specified Ve this quantity.
}
  \item{X}{Default is NULL, fixed effects design matrix in mixed effects model. If not specified, mixed model is assumed to have no fixed effect. If specified the number of rows of the matrix match to the total number of samples specified in the TrainSel function by the ntoselect parameter. Mostly useful if the rows of the X matrix are different within each block and an ordered samples are sought.
}

}
\details{
Either use a numerically coded (samples x features) matrix ('M') or a  positive definite relationship matrix of (samples x samples) ('K'). If providing both make sure these matrices include the same samples in the same order.}
\value{
A TrainSel_Data class object. This is actually a list object so you can create the data and then modify the contents of this as needed.
}

\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez
}
\examples{
\dontrun{
data(WheatData) # load data
{
#prepare data for mixed model based criteria (Mixed model I)
DataTS <- MakeTrainSelData(M = Wheat.M)
DataTS <- MakeTrainSelData(K = Wheat.K)
DataTS <- MakeTrainSelData(M = Wheat.M,K = Wheat.K)
X<-matrix(rnorm(100), ncol=2) # X has dimesions 50 x 2.
# this means we will be selecting 50 samples in total.
DataTS <- MakeTrainSelData(M = Wheat.M,K = Wheat.K, X=X)
DataTS <- MakeTrainSelData(M = Wheat.M,K = Wheat.K, X=X, lambda=2)
#lambda=2 means Ve is twice as Vk.
}
{
#prepare data for mixed model based criteria (Mixed model II)
#Two sets of samples will be selected.
Vk=matrix(.5, 2,2)+.5*diag(2) # Random effects in sample 1
# and sample 2 have correlation of 0.5.
#Error variances of sample 1 and sample 2
#are uncorrelated with sample 2 twice the error
#variance than sample 1.
Ve=diag(c(1,2))
DataTS <- MakeTrainSelData(M = Wheat.M, Vk=Vk, Ve=Ve)
DataTS <- MakeTrainSelData(M = Wheat.K, Vk=Vk, Ve=Ve)
# X has dimesions 50 x 2. i.e, the number of selected
# samples in sample 1 and sample 2 sums to 50.
# we will be selecting 50 samples in total.
DataTS <- MakeTrainSelData(M = Wheat.K, Vk=Vk, Ve=Ve, X=X)
DataTS <- MakeTrainSelData(M = Wheat.K, Vk=Vk, Ve=Ve, X=X, lambda=2)
#lambda =2 has the eefect of multiplying Ve by 2.
}

}
}

