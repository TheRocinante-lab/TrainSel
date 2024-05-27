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
           nelite = 8,
           nEliteSaved = 0,
           mutprob = .01,
           mutintensity = 2,
           niterations = 500,
           niterSANN = 0,
           stepSANN = 1,
           minitbefstop = 200,
           GABurnIn = 0,
           maxiterSANN = 1e+07,
           minitbefSANN = 0,
           SANNcooldown = 0,
           tolconv = 1e-7,
           progress=TRUE,
           nislands=1,
           npopIslands=0,
           neliteIslands = 0,
           nEliteSavedIslands = 0,
           niterIslands=0,
           niterSANNislands = 0,
           minitbefstopIslands = 20,
           verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npop}{
Size of the genetic algorithm (GA) population at each cycle}
\item{nelite}{
number of elite solutions selected at each cycle of the GA}
\item{nEliteSaved}{
number of elite solutions unaffected by SANN}
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
\item{GABurnIn}{Number of iterations of GA alone before starting GA+SANN. Only used if niterSANN > 0}
\item{maxiterSANN}{Maximum number of consecutive cycles of GA for which SANN will be performed. If 'maxiterSANN' is reached, SANN will be disabled until conditions for restarting it are met.}
\item{minitbefSANN}{If the GA does not improve the objective function more than 'tolconv' in the last 'minitbefSANN' iterations, start GA+SANN.}
\item{SANNcooldown}{If SANN is disabled after reaching 'maxiterSANN', it will not be enabled again until 'SANNcooldown' cycles of GA are performed.}
\item{tolconv}{Default value is 1e-7. Covergence tolerance parameter.}
\item{progress}{Default value is  TRUE. Shows the progress of the algorithm if TRUE.}
\item{nislands}{Default value is  1.  Values >1 will first run the optimization algorithm independently on nislands islands. The best solutions from these island populations will be used as an initial population for another round of optimization algoritm to output a final solution. Similar strategies can be implemented manually using the 'InitSols' argument of the 'TrainSel' function.}
\item{npopIslands}{Value of 'npop' used for the island populations.}
\item{neliteIslands}{Value of 'nelite' used for the island populations.}
\item{nEliteSavedIslands}{Value of 'nEliteSaved' used for the island populations.}
\item{niterIslands}{Value of 'niterations' used for the island populations.}
\item{niterSANNislands}{Value of 'niterSANN' used for the island populations.}
\item{minitbefstopIslands}{Value of 'minitbefstop' used for the island populations.}
\item{verbose}{Defaults to TRUE.}
}
\details{
Constructor function for the control parameters.
}
\value{
An object of class TrainSel_Control. This is actually a list and the output can be manuplated by manuplating the elements of this list.
}
\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez, Javier Fernandez-Gonzalez
}
\examples{
\dontrun{
control<-TrainSelControl()
control
control<-TrainSelControl(npop=30, nelite=5)
control
control$npop<-100
control$niterSANN=1
}
}
