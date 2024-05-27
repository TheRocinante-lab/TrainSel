\name{TimeEstimation}
\alias{TimeEstimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computational Time Estimation
}
\description{
Estimation of the time needed to run TrainSel optimization for a user defined optimization criterion and a given set of parameters and data.
}
\usage{
TimeEstimation(Data=NULL,
               Candidates=NULL,
               setsizes=NULL,
               settypes=NULL,
               control=NULL,
               Stat=NULL,
               n_average = 10,
               verbose = TRUE,
               nStat = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{List that contains the data to be used by the user defined optimization criteria. It can be null if the user defined function has no "Data" argument.
}
  \item{Candidates}{A list of vectors of integers. For 'DBL' type, the corresponding element in the list is a vector with the lower and upper bounds of the numeric variable of interest.)
}
  \item{setsizes}{A vector of integers. Number of samples from each of the elements in 'Candidates'.
}
  \item{settypes}{Six options are available:"BOOL" Boelian, "UOS" unordered sample, "UOMS" unordered multi-set, "OS" ordered set, "OMS" ordered multi-set, "DBL" numeric variables.
}
  \item{control}{NULL or a TrainSel_Control object created by 'TrainSelControl' or 'SetControlDefault'.
}
  \item{Stat}{User defined optimization criterion. For more details, see documentation on 'TrainSel' function.
}
  \item{n_average}{Number of times 'Stat' is run during the estimation of its computational time. The average of the 'n_average' runs will be used for the estimation of the total time needed.
}
  \item{verbose}{Defaults to true, resulting in additional details being printed to sdterr.
}
  \item{nStat}{Default is 1. Set >1 for multi-objective optimization, this is the number of objectives.
}
}
\details{
The error of the estimation is smaller for slower optimization criteria. The time calculated by this function assumes that convergence is not reached. If Stat is extremely fast (faster than a few miliseconds every time it is run), this function will underestimate the actual computational time.
}
\value{
A list with the following elements:

time_Stat_vector. It contains the 'n_average' values for the time (in seconds) needed to run 'Stat'.

eval_total. Number of times that 'Stat' is computed during optimization assuming no convergence.

time_total_seconds. Estimation of the time in seconds needed to complete the optimization process.
}

\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez, Javier Fernandez-Gonzalez
}
\examples{
\dontrun{
#Toy optimization criterion with a mostly controlled compuational time
Stat <- function(soln, Data) {
  sleepTime <- Data[["sleepTime"]]
  Sys.sleep(sleepTime)
  return(sum(soln))
}

Data <- list(sleepTime = 1)

#SetControlDefault allows to easily set parameters from tested configurations
#type ?SetControlDefault for more details
control <- SetControlDefault(size="small",
                             complexity = "low_complexity",
                             verbose=TRUE)

Candidates <- list(1:100)
setsizes <- c(50)
settypes <- c("UOS")

#Estimation of time in seconds needed for optimization
TimeResults <- TimeEstimation(Data=Data,
                              Candidates=list(1:100),
                              setsizes=c(50),
                              settypes=c("UOS"),
                              control=control,
                              Stat=Stat)
print(TimeResults)
}
}

