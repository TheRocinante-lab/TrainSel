\name{TrainSelControl}
\alias{TrainSelControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Settings for the optimization algorithm
}
\description{
Settings for the optimization algorithm
}
\usage{
TrainSelControl(npop = 300,
           nelite = 10,
           mutprob = .01,
           mutintensity = 2,
           niterations = 500,
           niterSANN = 200,
           stepSANN = 1e-2,
           minitbefstop = 200,
           tolconv = 1e-7,
           progress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npop}{
Size of the genetic algorithm (GA) population at each cycle}
\item{nelite}{
number of elite solutions selected at each cycle of the GA}
\item{mutprob}{
Chance of mutation for each generated solution}
  \item{mutintensity}{
Number of mutations for each generated solution}
  \item{niterations}{
Maximum Number of GA iterations}
  \item{niterSANN}{
Number of simmulated annealing iterations for each elite solution at each cycle of the GA.}
\item{stepSANN}{Default is 1e-2. A parameter for the cooling of the SANN, for smaller values the system is cooled slower.}
\item{minitbefstop}{If the GA does not improve the objective function more than 'tolconv' in the last 'minitbefstop' iterations the algorithm will stop.}
\item{tolconv}{Default value is 1e-7. Covergence tolerance parameter.}
\item{progress}{Default value is  TRUE. Shows the progress of the algorithm if TRUE.}
}
\details{
Constructor function for the control parameters.
}
\value{
An object of class TrainSel_Control. This is aculally a list sos the output can be manuplated by manuplating the elements of this list.
}
\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez
}
\examples{
\dontrun{
control<-TrainSelControl()
control
control<-TrainSelControl(npop=30, nelite=5)
control
control$npop<-100
control$niterSANN=1000
}
}
