\name{TrainSel}
\alias{TrainSel}
\title{
TrainSel
}
\description{
A function for mixed integer programming. We assume that the maximization of the objective functions is of interest. Objective function can be contructed to accomodate binary ("BOOL"), integers (ordered-unordered, set-multi set) selected from candidate set of integers, and continuous variables restricted by upper and lower bounds.
}
\usage{
TrainSel(Data = NULL,
           Candidates = NULL,
           setsizes = NULL,
           ntotal = NULL,
           settypes = NULL,
           Stat = NULL,
           nStat= 1,
           Target = NULL,
           control = NULL)
}
\arguments{
  \item{Data}{This is either a 'TrainSel_Data' object created by 'MakeTrainSelData' function or it is a list that contains the data to be used by the user defined optimization criteria. It can be null if the user defined function only has no "Data" argument (for instance, when the data is called from the R environment).}
  \item{Candidates}{A list of vectors of integers, for 'DBL' type the corresponding element in the list are the lower and upper bounds of the numeric variable of interest.}
\item{setsizes}{A vector of integers. Number of samples in each of the samples.}
\item{ntotal}{An integer value. The total number of unique integers in 'BestSoln_int' in the in the final solution. Default value is 0 which means no ntotal constraint is applied. When ntotal>0 and number of unique integers in the solution is not equal to ntotal, the Stat value is set to -1e_10*abs("number of unique integers in the solution"-ntotal).}
\item{settypes}{Six options are available:"BOOL" Boelian, "UOS" unordered sample, "UOMS" unordered multi-set, "OS" ordered set, "OMS" ordered multi-set, "DBL" numeric variables.
}
\item{Stat}{Either NULL (in that case CDMin is used) or a user defined function that takes an integer vector or a vector of numeric variables and possibly a Data argument. The first input is used as the vector of integers, second input is the vector of numeric values and the last input if provided is used as Data. For univariate optimization this function returns a numeric value which we want to maximize with respect to the integer vector. Note that integer vector is considered to be the concatenation of the "OS", "UOS", "OMS", "UOMS" components provided in 'settypes', and the numeric vector corresponds to the concanetation of the 'DBL' components. There is also minimal suppoert for muti-objective optimization (see details and examples).
}
\item{nStat}{Default is 1. Set >1 for MO optimization, this is the number of objectives.}
\item{Target}{Default is NULL. An integer vector. When used with CDMin, the CDMin is calculated on the specified Target. When used with a user defined function with  non NULL Data argument, Target is added to the Data list and can be accessed by the function in 'Stat' with "Data[["Target"]]".
}
\item{control}{NULL or a TrainSel_Control object created by 'TrainSelControl'.
}
}
\details{

TrainSel package is composed of a set of programs that work for optimally selecting ordered or unordered samples (set or multi-set) from a set of candidates, or several ordered or unordered samples (set or multi-set) from several sets of candidates based on of one (or more) optimization criterion. There is also some support for solving mixed-integer optimization problems and muti-objective optimization.



1- Built in criteria

The generic form of a Gaussian linear mixed model (MM) can be written as: \deqn{y=X\beta+Zg+We,} where \eqn{y} is the \eqn{n \times 1} response variable, \eqn{X} is the \eqn{n \times p} design matrix of \eqn{p \times 1} the fixed effects \eqn{\beta}, \eqn{Z} is the \eqn{n \times q} design matrices of the \eqn{q \times 1} random effect \eqn{g,}  \eqn{W} is the \eqn{n \times r} design matrices of the \eqn{r\times 1} random effect \eqn{e.} The random effect and the residual are independently distributed and are normally distributed as \eqn{g\sim N_{q \times 1}(0_{q \times 1}, G)} and \eqn{e\sim N_{r \times 1}(0_{r \times 1}, D).}


To estimate the precision of the fixed and random effects we can use  the coefficient matrix \cite{Henderson:1975}. For a given design, the the coefficient matrix is written as \tabular{rccl}{
 |\tab \eqn{X' (WDW')^{-1}X,} \tab \eqn{X'(WDW')^{-1}Z}  \tab | \cr
 |\tab  \eqn{Z' (WDW')^{-1}X,} \tab \eqn{Z'(WDW')^{-1}Z + G^{-1}}  \tab |
}

We call the inverse of the coefficient matrix \eqn{H} and this matrix is assumed to be partitioned as
\tabular{rccl}{
    | \tab \eqn{H_{11},} \tab \eqn{H_{12}} \tab | \cr
    | \tab \eqn{H'_{12},} \tab \eqn{H_{22}} \tab |
}
where \eqn{H_{11}}, \eqn{H_{12}}, and \eqn{H_{22}} are, respectively, \eqn{p\times p,} \eqn{p\times q,} and \eqn{q\times q} sub-matrices that correspond to the blocks of the coefficient matrix. Calculation of the \eqn{H} matrix does not involve the response vector, it only depends on the design and the covariance parameters \eqn{D} and \eqn{G.} Therefore, if we assume we have some prior knowledge about these matrices we can use \eqn{H} matrix to derive optimization criteria. The coefficient of determination matrix for the random effect can be expressed as \deqn{CD=(G-H_{22})/G} where the division  is an elementwise division. Each diagonal element of CD is between 0 and 1. CDMin design criterion is obtained by taking the minimum of certain diagonal elements of this matrix.

"CDMin" is written for two forms of MMs and assumes that the focus is on estimating the random effects. We will refer to the specific MMs covered under 'TrainSel' as MM-I and MM-II. Lets assume that the matrices \eqn{K} and \eqn{R} are kernel matrices (covariance matrices with some standardization over the diagonals). We also assume that the design matrix Z is a permutation matrix (with exactly one 1 at each row with remaining elements zero).

MM-I is parametrized as follows: \deqn{y=X\beta+Z(g+e),} where \eqn{y} is the \eqn{n \times 1} response variable, \eqn{X} is the \eqn{n \times p} design matrix of \eqn{p \times 1} the fixed effects \eqn{\beta}, \eqn{Z} is the \eqn{n \times q} design matrices of the \eqn{q \times 1} random effects \eqn{g}  and \eqn{e.} The random effects \eqn{g} and \eqn{e} are independently distributed and are normally distributed as \eqn{g\sim N_{q \times 1}(0_{q \times 1}, \sigma^2_k K)} and \eqn{e\sim N_{q \times 1}(0_{q \times 1}, \sigma^2_e R)} correspondingly.

MM-II is parametrized as follows: \deqn{y=X\beta+Z(g+e),} where \eqn{y} is the \eqn{n \times 1} response variable, \eqn{X} is the \eqn{n \times p} design matrix of \eqn{p \times 1} the fixed effects \eqn{\beta}, \eqn{Z} is the \eqn{n \times qd} design matrices of the \eqn{qd \times 1} random effects \eqn{g}  and \eqn{e} \eqn{ qd\times 1} residual vector. The random effect and the residual are independently distributed and are normally distributed as \eqn{g\sim N_{qd \times 1}(0_{qd \times 1}, K \otimes V_k)} and \eqn{e\sim N_{qd \times 1}(0_{qd \times 1}, R \otimes V_e)} correspondingly. \eqn{V_k} and \eqn{V_e} are \eqn{d times d} covariance matrices.

We use MM-I when we know the the relationship or resemblance of the random effects that are expressed in the relationship matrices \eqn{K} and \eqn{R.} The \eqn{K} matrix can be obtained given the features that describe the levels of the random effects \eqn{g} and the usual assumption is that \eqn{e} are independent and identically distributed. The design problem we seek to solve is to find a good design matrix of the random effects. Two cases arise: If the design matrix of the fixed effects \eqn{X} are also specified so that the order of the random effects matters we are in the ordered sample case. If the fixed effects matrix is not specified or even if it is specified the ordering of the random effects does not matter then we are in the unordered sample case. In either case, we are seeking one set of samples from one set of candidates.

The difference MM-II from MM-I is that in MM-II we want to select a sample from a set of candidates in several related conditions. The relatedness of the \eqn{d} different conditions are expressed in the covariance matrices \eqn{V_k} and \eqn{V_e}.  As with MM-I the existance of a fixed effects design matrix might require an ordered sample.


Once the CD matrix for a design is calculated,  a scalar measure is calculated from this matrix is used as the optimization criterion. If no 'Target' is specified, mean of the diagonals of this matrix is taken as the scalar optimization criterion. If the 'Target' is specified, then the mean is taken over the diagonals that correspond to the 'Target' is used as the optimization criterion.

All of the parameters 'Candidates' and 'Target' are integer vectors of subset of the random effect levels. If in doubt take a look at the the 'labels' in TrainSel_Data" object and see examples for more details and usage.

2- User defined criteria

Most flexible use of TrainSel is by using user defined criteria. This requires the user to  define function. This function follows has two or  three inputs, 'soln_int'  and 'soln_dbl' only or 'soln_int', 'soln_dbl' and 'Data'. 'soln_int' is an integer vector, soln_dbl is a vector of real numbers. Data is a list object. The user defined function should return a scalar value when provided with admisable 'soln_int', 'soln_dbl' and 'Data'. TrainSel can then be used to look for 'soln_int' or 'soln_dbl' that maximize the value of the returned values from this function.


See examples and the package vignettes for details and usage.
}
\value{A list of TrainSelOut class. Best solution found is in "BestSol_Int", "BestSol_Dbl" and best objective value is in "BestVal". The value of the objective function over the iterations is in "maxvec". "convergence" is an integer with value 0 or 1, 0 means "control$niterations" reached before convergence (which is defined as no improvement in the last "conrol$minitbefstop"). 1 means the convergence is achieved (the magnitude of the improvement in the last miniterbefstop iterations is less than 'tolconv').}
\references{
Deniz Akdemir
}
\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez
}

\examples{
\dontrun{

###changing control so that the algorithm runs
##in a short time. You need more than 5 iterations
##for convergence!
control<-TrainSelControl()
control$niterations=5
##############################
####1-Select an unordered set optimizing the
####CDMin criteria.
####We use the 'WheatData' for this example.
####There are 200 genotypes
####each genotype has 4670 markers (Wheat.M),
####therefore we have a good idea how these
####200 genotypes are related. We are planning
####to select a subset of size 50 out of
####these 200 genotypes as a training
####population using the CDMin criterion.

data(WheatData)
control<-TrainSelControl()
control$niterations=10 # increase this for convergence.
TSDAT <- MakeTrainSelData(M = Wheat.M)
sol1<-TrainSel(Data=TSDAT,Candidates=list(1:200),
setsizes=c(50), settypes="UOS", Stat=NULL,
control=control)
sol1
sol1$BestVal
sol1$BestSol_int
rownames(Wheat.M)[sol1$BestSol_int]
plot(sol1$maxvec)


##############################
####2- Select an ordered set optimizing the
####CDMin criteria.
####Same data as 1. We are planning
####to select a subset of size 50 out of
####these 200 genotypes as a training
####population using the CDMin criterion.
####A design matrix for the fixed effects
####X (50 x p) is specified in the data
####statement and there is heterogeneity
####among the rows of X, i.e., order is
####important.
X<-as.matrix(cbind(1, expand.grid(1:5,1:10)));head(X);
TSDAT <- MakeTrainSelData(M = Wheat.M,X=X)
sol2<-TrainSel(Data=TSDAT,Candidates=list(1:200),
setsizes=c(50),settypes="OS",
Stat=NULL,control=control)
print(sol2)
sol2$BestSol_int
plot(sol2$maxvec)



##############################
####3- Select 2 unordered sets optimizing the
####CDMin criteria.
####Same data as 1. We are planning
####to select two subsets of size 20 and 30 out of
####these 200 genotypes  as a training
####population using the CDMin criterion.
####These samples will be evaluated under
####conditions which are known to be realted by the
####covariances \eqn{V_k} and \eqn{V_e}.

Ve=diag(2) # unrelated residuals among the conditions

Vk=.5*diag(2)+matrix(.5,2,2) # the condition 1 and
                            #2 are 50\% correlated

TSDAT <- MakeTrainSelData(M = Wheat.M,Vk=Vk,Ve=Ve)
sol3<-TrainSel(Data=TSDAT,Candidates=list(1:nrow(Wheat.M),(1:nrow(Wheat.M))+nrow(Wheat.M)),
setsizes=c(20,30), settypes=c("UOS", "UOS"),
Stat=NULL,control=control)
sol3$BestSol_int ##
# note that the integers 201:400 refer to
# individuals 1:200 in second environment
plot(sol3$maxvec)

##############################
####4- Select D-Optimal subset of 50 samples
####from 200 samples in the 'WheatData' dataset
####using 10 principal components.
Doptfunc<-function(soln,Data){
determinant(crossprod(Data[[1]][soln,]),logarithm = TRUE)$modulus
}
PCData<-svd(scale(Wheat.M, scale=FALSE,center=TRUE))$u[,1:10]
sol4<-TrainSel(Data=list(PCData), Candidates=list(1:200),
setsizes=c(50), Stat=Doptfunc,settypes="UOS", control=control)
sol4
}

}
