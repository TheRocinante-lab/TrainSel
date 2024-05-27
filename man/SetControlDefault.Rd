\name{SetControlDefault}
\alias{SetControlDefault}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Settings for the optimization algorithm
}
\description{
Creates an object of the class 'TrainSel_Control'. It allows to easily get a configuration well suited for the user needs.
}
\usage{
SetControlDefault(size="medium",
                  complexity = "low_complexity",
                  verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{It has four options: "demo", "small", "medium", "large". Larger values are needed if there is a high number of combinations in the search space. Some factors that increase it are a higher number of candidates, allowing for repetition, ordered optimization, etc. The only setting that works without a valid license key is "demo".
}
  \item{complexity}{it has two possible values: "low_complexity" and "high_complexity". The former only employs a GA for fast convergence, but it presents a high risk of getting stuck in a local maximum. The latter uses a mixture of islands, GA and SANN to reduce the risk of converging at a suboptimal solution, with the disadvantage of a slower learning rate.
}
  \item{verbose}{Defaults to true, resulting in additional details being printed to sdterr.
}

}
\details{
The best combination of "size" and "complexity" can be known from previous experience with the optimization problem. When working with unfamiliar problems, we suggest the following steps:

1) Identify if computational time is a limiting factor. The time needed to perform optimization can be found using "TimeEstimation" function.

2) If optimization is fast, setting "size" to "large" and "complexity" to "high_complexity" will ensure the best results for single-objective optimization. For multi-objective optimization (nStat > 1), "complexity" should always be set to "low_complexity".

3) If optimization is slow and you are limited by computational time, set "complexity" to "low_complexity" and "size" as large as possible for the present time limitations.

4) If convergence is reached early, you can decrease "size" and increase "complexity" next time optimization is run for a similar optimization problem.

5) If the algorithm is far from convergence ("maxvec" in the output of "TrainSel" was increasing relatively fast when reaching the limit of iterations), "size" has to be increased. If "size" was already set to "large" and convergence is not reached, it is recomended to manually increase "niterations" in the output of SetControlDefault.
}



\value{
An object of class TrainSel_Control. This is actually a list and the output can be manuplated by manuplating the elements of this list.
}

\author{
Deniz Akdemir, Simon Rio, Julio Isidro Sanchez, Javier Fernandez-Gonzalez
}
\examples{
\dontrun{
control <- SetControlDefault(size="small",
                             complexity = "low_complexity",
                             verbose=TRUE)
control
control$npop<-300
control$niterations<-300
control
}
}

