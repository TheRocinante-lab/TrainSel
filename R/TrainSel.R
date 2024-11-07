########################TrainSelData
MakeTrainSelData <-
  function(M = NULL,
           ##Features Matrix Samples
           K = NULL,
           ##Relationship Matrix Samples
           R = NULL,
           ##Relationship Matrix errors of Samples
           Vk = NULL,
           ##Relationship Matrix Blocks
           Ve = NULL,
           ##Relationship Matrix errors of Blocks
           lambda = NULL,
           ##ratio of Ve to Vk
           X = NULL) {

    if ((is.null(M) & is.null(K))) {
      stop("At least one of M (features) or K (similarity matrix) has to be present.")
    }
    if ((is.null(M) & !is.null(K))) {
      M<-svd(K)$u
      rownames(M)<-rownames(K)
    }

    namesinM<-rownames(M)
    if (is.null(namesinM)){
      rownames(M)<-1:nrow(M)
      namesinM<-rownames(M)
    }
    ####Mixed model stats
    if (is.null(K)) {
      K = cov(t(M))
      K =nearPDc(K / mean(diag(K)))
    }
    if (is.null(rownames(K))){rownames(K)<-colnames(K)<-namesinM}
    if (is.null(R)) {
      message("Assuming that the residual covariance in the mixed model is identity.")
      R = diag(nrow(K))
      rownames(R) <- colnames(R) <- namesinM
    }
    if ((!is.null(Ve) & !is.null(Vk))){
      if (is.null(lambda)) {
        lambda = 1
      }
      if (is.null(rownames(Vk))){rownames(Vk)<-1:nrow(Vk)}
      if (is.null(rownames(Ve))){rownames(Ve)<-rownames(Vk)}


      bigK <- kronecker(Vk, K, make.dimnames = TRUE)
      bigR <- kronecker(lambda * Ve, R, make.dimnames = TRUE)
      labels<-data.frame(intlabel=1:nrow(bigK),names=rownames(bigK))
    } else {
      if (is.null(lambda)) {
        lambda = 1
      }
      bigK <- K
      bigR <- lambda * R
      labels<-data.frame(intlabel=1:nrow(bigK),names=rownames(bigK))

    }

    if (is.null(X)){
      return(
        list(
          G = bigK,
          R = bigR,
          lambda = lambda,
          labels=labels,
          Nind=nrow(K),
          class = "TrainSel_Data"
        ))
    }

    if (!is.null(X)){
      return(
        list(
          X=X,
          G = bigK,
          R = bigR,
          lambda = lambda,
          labels=labels,
          Nind=nrow(K),
          class = "TrainSel_Data"
        ))
    }
  }













##########################TrainSel
TrainSelInner <-
  function(Data = NULL,
           Candidates = NULL,
           setsizes = NULL,
           ntotal = NULL,
           settypes = NULL,
           Stat = NULL,
           nStat= 1,
           Target = NULL,
           control = NULL,
           InitSol=NULL,
           Username=NULL,
           Password=NULL,
           Verbose=TRUE) {

    if (is.null(Username) | is.null(Password)) {
      FullVersion <- FALSE
      if (Verbose) {
        print("Using demo version.")
        write("Limited to selecting 100 integer and 3 double solutions with hyperparameters in SetControlDefault(size=\"demo\")", stdout())
        print("If you are interested in upgrading to full version, please contact us at \x6A.\151s\x69d\162o\x40u\160m\x2Ee\163")
      }
      Username<-0
      Password<-0
    } else {
      FullVersion <- o_8ac3596136014076863521e808b79394(Username, Password, Verbose)
      if (!FullVersion) {Sys.sleep(5)}
    }




    ########Check types
    ###Data can only be NULL or list
    if ((!is.null(Data) & !is.list(Data))){
      stop("Data is either NULL or a list")
    }


    if (is.null(ntotal)) {
      ntotal <- 0
    }

    if (!is.function(Stat)){
      if (is.null(Stat))
      {
        CD=TRUE
        Stat=function(){}
      }  else {
        stop('Stat is either NULL (meaning CDMEAN is used) or it is not a function')
      }
    } else {
      CD=FALSE
    }
    if (is.null(Data) & CD){
      stop("no data for mixed model based method")
    }

    if (is.null(Target)) {
      Target <- as.integer(c())
    }


    if (is.null(InitSol)){
      InitSol<-list(solnIntMat=matrix(as.integer(c()), ncol=0, nrow=0),solnDBLMat=matrix(as.numeric(c()), ncol=0, nrow=0) )
    } else {
      solnIntMat<-InitSol[["solnIntMat"]]
      solnDBLMat<-InitSol[["solnDBLMat"]]
      if (is.null(solnIntMat)){
        solnIntMat=matrix(as.integer(c()), ncol=0, nrow=0)
      }
      if (is.null(solnDBLMat)){
        solnDBLMat=matrix(as.numeric(c()), ncol=0, nrow=0)
      }
      InitSol<-list(solnIntMat=solnIntMat,solnDBLMat=solnDBLMat)
    }

    if (nStat==1){
      out<-TrainSelC(Data=Data, CANDIDATES =Candidates,setsizes =setsizes,settypes=settypes,Stat = Stat,CD=CD,Target=Target,control=control, ntotal=ntotal, InitSol=InitSol, o_4c3b474d6e250cd9804178f26306d565=Username, o_5c85d553f3828f5855fe83513707eda8=Password, o_d5ead238fb380c7d6fa344cc58cb043a=FullVersion)
      class(out)<-"TrainSelOut"
    } else {
      out<-TrainSelCMOO(Data=Data, CANDIDATES =Candidates,setsizes =setsizes,settypes=settypes,Stat = Stat,nstat=nStat, control=control, InitSol=InitSol, o_4c3b474d6e250cd9804178f26306d565=Username, o_5c85d553f3828f5855fe83513707eda8=Password, o_d5ead238fb380c7d6fa344cc58cb043a=FullVersion)
      class(out)<-"TrainSelOut"
    }
    return(out)
  }









TrainSel<-
  function(Data = NULL,
           Candidates = NULL,
           setsizes = NULL,
           ntotal = NULL,
           settypes = NULL,
           Stat = NULL,
           nStat= 1,
           Target = NULL,
           control = NULL,
           InitSol=NULL,
           Username=NULL,
           Password=NULL,
           Verbose=TRUE
           ) {


    ####  control is null or  and object not of type TrainSel_control
    if ((!is.null(control) & (class(control)[1]!="TrainSel_Control"))){
      stop("control can be NULL or a TrainSel_Control object.")
    }

    if (is.null(control)){
      control=TrainSelControl()
    }

    nislands<-control[["nislands"]]


    nDBL <- sum(settypes == "DBL")
    if (nislands==1){

        out<-TrainSelInner(Data = Data,
                           Candidates = Candidates,
                           setsizes = setsizes,
                           ntotal = ntotal ,
                           settypes = settypes,
                           Stat = Stat,
                           nStat= nStat,
                           Target = Target,
                           control = control,
                           InitSol=InitSol,
                           Username=Username,
                           Password=Password,
                           Verbose=Verbose)




    } else {
      controlInner<-control
      #controlInner$progress=FALSE
      #Parameters for the islands
      controlInner$niterations<-control$niterIslands
      controlInner$minitbefstop<-control$minitbefstopIslands
      controlInner$niterSANN<-control$niterSANNislands
      controlInner$nEliteSaved<-control$nEliteSavedIslands
      controlInner$nelite<-control$neliteIslands
      controlInner$npop <- control$npopIslands
      #controlInner$dynamicNelite<-FALSE


        solList<-lapply(1:nislands, function(x){write(paste0("\nIsland ", x, ":"), stderr());
                                                              TrainSelInner(Data = Data,
                                                              Candidates = Candidates,
                                                              setsizes = setsizes,
                                                              ntotal = ntotal ,
                                                              settypes = settypes,
                                                              Stat = Stat,
                                                              nStat= nStat,
                                                              Target = Target,
                                                              control = controlInner,
                                                              Username=Username,
                                                              Password=Password,
                                                              Verbose=Verbose)})




      write("\nFinished Islands. Starting Final Step:", stderr())


      solnIntMat<-NULL
      solnDBLMat<-NULL
      condint<-length(intersect(settypes,c("BOOL","UOS","UOMS","OS","OMS")))>0
      conddbl<-length(intersect(settypes,c("DBL")))>0
      if (condint){
        solnIntMat<-lapply(solList, function(x){x$BestSol_int})
        solnIntMat<-Reduce("cbind", solnIntMat)
      }
      if (conddbl){
        solnDBLMat<-lapply(solList, function(x){x$BestSol_DBL})
        solnDBLMat<-Reduce("cbind", solnDBLMat)
      }

      if (condint>0 & conddbl>0){
        InitSol<-list(solnIntMat=solnIntMat, solnDBLMat=solnDBLMat)
      } else if (condint==0 & conddbl>0) {
        InitSol<-list(solnDBLMat=solnDBLMat)
      } else if (condint>0 & conddbl==0) {
        InitSol<-list(solnIntMat=solnIntMat)
      } else {
        InitSol=NULL
      }



        out<-TrainSelInner(Data = Data,
                           Candidates = Candidates,
                           setsizes = setsizes,
                           ntotal = ntotal ,
                           settypes = settypes,
                           Stat = Stat,
                           nStat= nStat,
                           Target = Target,
                           control = control,
                           InitSol=InitSol,
                           Username=Username,
                           Password=Password,
                           Verbose=Verbose)


    }
    return(out)
  }





TrainSelControl <-  function(npop = 100,
                             npopIslands=0,
                             nelite = 8,
                             neliteIslands = 0,
                             mutprob = .01,
                             mutintensity = 2,
                             niterations = 100,
                             niterSANN = 0,
                             niterSANNislands = 0,
                             stepSANN = 1,
                             minitbefstop = 20,
                             minitbefstopIslands = 20,
                             tolconv = 1e-7,
                             progress=TRUE,
                             nislands=1,
                             niterIslands=0,
                             nEliteSaved = 0,
                             nEliteSavedIslands = 0,
                             GABurnIn = 0,
                             maxiterSANN = 1e+07,
                             minitbefSANN = 0,
                             SANNcooldown = 0,
                             verbose = TRUE)  {


  if (verbose) {
    write("*************************************************************************", stderr())
    write("We recommend using SetControlDefault() and TimeEstimation() functions", stderr())
    write("to tune the hyperparameters in the TrainSel_Control object", stderr())
    write("*************************************************************************", stderr())
  }

  return(structure(
    list(
      npop = npop,
      npopIslands = npopIslands,
      nelite = nelite,
      neliteIslands=neliteIslands,
      mutprob = mutprob,
      mutintensity = mutintensity,
      niterations = niterations,
      niterSANN = niterSANN,
      niterSANNislands = niterSANNislands,
      stepSANN = stepSANN,
      minitbefstop = minitbefstop,
      minitbefstopIslands = minitbefstopIslands,
      tolconv = tolconv,
      progress=progress,
      nislands=nislands,
      niterIslands=niterIslands,
      nEliteSaved=nEliteSaved,
      nEliteSavedIslands=nEliteSavedIslands,
      GABurnIn=GABurnIn,
      maxiterSANN=maxiterSANN,
      minitbefSANN = minitbefSANN,
      SANNcooldown=SANNcooldown
    ),
    class = c("TrainSel_Control")
  )
  )


}




PamSel<-function(M=NULL,distancemethod=c("euclidean",  "maximum", "manhattan", "canberra", "binary", "minkowski"),ntoselect=NULL){
  stopifnot(distancemethod%in%c("euclidean",  "maximum", "manhattan", "canberra", "binary", "minkowski"))
  D=as.matrix(dist(M, method = distancemethod))
  rownames(D)<-colnames(D)<-rownames(M)
  pamx <- cluster::pam(D, k=ntoselect)
  namesselected <- rownames(pamx$medoids)
  return(which(rownames(M)%in%namesselected))
}






#Generate control object with default parameters for several scenarios
SetControlDefault <- function(size="medium",
                              complexity = "low_complexity",
                              verbose=TRUE) {

  CNTRL<-TrainSelControl(verbose = FALSE)

  if (size == "small") {

    if (complexity == "low_complexity") {
      CNTRL$npop = 200
      CNTRL$nelite = 5
      CNTRL$niterations = 200
      CNTRL$niterSANN = 0
      CNTRL$minitbefstop = 100
      CNTRL$nislands=1


    } else if (complexity == "high_complexity") {



      CNTRL$npop = 200
      CNTRL$nelite = 5
      CNTRL$niterations = 200
      CNTRL$niterSANN = 1
      CNTRL$stepSANN = 1
      CNTRL$nEliteSaved = 3
      CNTRL$maxiterSANN = 1
      CNTRL$SANNcooldown = 1
      CNTRL$minitbefstop = 100
      CNTRL$nislands=5
      CNTRL$npopIslands=100
      CNTRL$neliteIslands = 5
      CNTRL$niterIslands=50
      CNTRL$minitbefstopIslands = 20

    } else {
      stop("Invalid value for \"complexity\". It should be \"low_complexity\" or \"high_complexity\"")
    }



  } else if (size == "medium") {



    if (complexity == "low_complexity") {

      CNTRL$npop = 300
      CNTRL$nelite = 8
      CNTRL$niterations = 500
      CNTRL$niterSANN = 0
      CNTRL$minitbefstop = 200
      CNTRL$nislands=1



    } else if (complexity == "high_complexity") {




      CNTRL$npop = 300
      CNTRL$nelite = 8
      CNTRL$niterations = 500
      CNTRL$niterSANN = 1
      CNTRL$stepSANN = 1
      CNTRL$nEliteSaved = 4
      CNTRL$maxiterSANN = 1
      CNTRL$SANNcooldown = 1
      CNTRL$minitbefstop = 200
      CNTRL$nislands=8
      CNTRL$npopIslands=100
      CNTRL$neliteIslands = 5
      CNTRL$niterIslands=100
      CNTRL$minitbefstopIslands = 20




    } else {
      stop("Invalid value for \"complexity\". It should be \"low_complexity\" or \"high_complexity\"")
    }


  } else if (size == "large") {



    if (complexity == "low_complexity") {


      CNTRL$npop = 300
      CNTRL$nelite = 8
      CNTRL$niterations = 1500
      CNTRL$niterSANN = 0
      CNTRL$minitbefstop = 200
      CNTRL$nislands=1



    } else if (complexity == "high_complexity") {


      CNTRL$npop = 300
      CNTRL$nelite = 12
      CNTRL$niterations = 1500
      CNTRL$niterSANN = 1
      CNTRL$stepSANN = 1
      CNTRL$nEliteSaved = 6
      CNTRL$maxiterSANN = 1
      CNTRL$SANNcooldown = 1
      CNTRL$minitbefstop = 200
      CNTRL$nislands=12
      CNTRL$npopIslands=100
      CNTRL$neliteIslands = 5
      CNTRL$niterIslands=200
      CNTRL$minitbefstopIslands = 20






    }
  } else if (size == "demo") {

    if (complexity == "low_complexity") {

      CNTRL$npop = 100
      CNTRL$nelite = 5
      CNTRL$niterations = 100
      CNTRL$niterSANN = 0
      CNTRL$minitbefstop = 50
      CNTRL$nislands=1



    } else if (complexity == "high_complexity") {
      stop("Demo does not suppport \"high_complexity\" settings")
    } else {
      stop("Invalid value for \"complexity\". It should be \"low_complexity\" or \"high_complexity\"")
    }


  } else {
    stop("Invalid value for \"size\". It should be \"small\", \"medium\" or \"large\"")
  }

  if (verbose) {

    write("*************************************************************************", stderr())
    write("Current settings:", stderr())
    write(paste0("- Size: ", size), stderr())
    write(paste0("- Complexity: ", complexity), stderr())

    write("Guidelines:", stderr())
    write("- You may use TimeEstimation function to estimate computational time", stderr())
    write("- Size should be as large as allowed by available computational resources", stderr())
    write("- If you are limited by computational time and don't reach convergence,
  use the low_complexity setting. Otherwise, high_complexity is preferred", stderr())
    write("- low_complexity is always preferred for in multi-objective optimization.", stderr())
    write("*************************************************************************", stderr())
  }
  return(CNTRL)
}





TimeEstimation <- function(Data=NULL,
                           Candidates=NULL,
                           setsizes=NULL,
                           settypes=NULL,
                           control=NULL,
                           Stat=NULL,
                           n_average = 10,
                           verbose = TRUE,
                           nStat = 1) {


  nDBL <- sum(settypes == "DBL")


  mc.cores <- 1

  time_one_run_vector <- c()

  for (i in 1:n_average) {

    soln_int <- c()
    soln_dbl <- c()

    for (i in 1:length(settypes)) {
      type <- settypes[i]

      if (type == "DBL") {
        lower_bound <- Candidates[[i]][1]
        upper_bound <- Candidates[[i]][2]
        cands <- seq(lower_bound, upper_bound, abs(upper_bound-lower_bound)/100)
        soln_dbl <- c(soln_dbl, sample(cands, setsizes[i], replace = TRUE)) #replace = TRUE to be robust to allowing repetition
      } else {
        soln_int <- c(soln_int, sample(Candidates[[i]], setsizes[i], replace = TRUE)) #replace = TRUE to be robust to allowing repetition
      }

    }



    #only int
    if (nDBL == 0) {

      #No user-provided Data
      if (is.null(Data)) {
        start_time <- Sys.time()
        output <- Stat(soln_int)
        end_time <- Sys.time()
        time_one_run_vector <- c(time_one_run_vector, difftime(end_time, start_time, units = "secs"))


        #With user-provided Data
      } else {

        start_time <- Sys.time()
        output <- Stat(soln_int, Data)
        end_time <- Sys.time()
        time_one_run_vector <- c(time_one_run_vector, difftime(end_time, start_time, units = "secs"))

      }

      #only dbl
    } else if (nDBL == length(settypes)) {

      #No user-provided Data
      if (is.null(Data)) {

        start_time <- Sys.time()
        output <- Stat(soln_dbl)
        end_time <- Sys.time()
        time_one_run_vector <- c(time_one_run_vector, difftime(end_time, start_time, units = "secs"))

        #With user-provided Data
      } else {

        start_time <- Sys.time()
        output <- Stat(soln_dbl, Data)
        end_time <- Sys.time()
        time_one_run_vector <- c(time_one_run_vector, difftime(end_time, start_time, units = "secs"))

      }

      #double plus int
    } else {

      #No user-provided Data
      if (is.null(Data)) {

        start_time <- Sys.time()
        output <- Stat(soln_int, soln_dbl, Data)
        end_time <- Sys.time()
        time_one_run_vector <- c(time_one_run_vector, difftime(end_time, start_time, units = "secs"))


        #With user-provided Data
      } else {

        start_time <- Sys.time()
        output <- Stat(soln_int, soln_dbl, Data)
        end_time <- Sys.time()
        time_one_run_vector <- c(time_one_run_vector, difftime(end_time, start_time, units = "secs"))

      }

    }

  }


  time_one_run <- mean(time_one_run_vector)


  if (nStat == 1) {

    eval_islands <- 0

    if (control$nislands > 1) {
      eval_islands_GA <- control$nislands*(ceiling(control$npopIslands/mc.cores) + control$niterIslands*ceiling((control$npopIslands-control$neliteIslands)/mc.cores))


      if ((control$npopIslands-control$neliteIslands)%%mc.cores == 0) {
        idle_cores <- 0
      } else {
        idle_cores <- (mc.cores - (control$npopIslands-control$neliteIslands)%%mc.cores)
      }
      #idle_cores
      efficiency_islands_GA <- ((control$npopIslands-control$neliteIslands))/
        ((control$npopIslands-control$neliteIslands)+idle_cores)


      counter <- 0
      eval_islands_SANN <- 0
      efficiency_islands_SANN <- 1
      if (control$niterSANNislands > 0) {
        while (counter <= control$niterIslands) {
          eval_islands_SANN <- eval_islands_SANN + control$maxiterSANN*control$niterSANNislands*ceiling((control$neliteIslands -  control$nEliteSaved)/mc.cores)
          counter <- counter + control$maxiterSANN + control$SANNcooldown
        }
     # }
      eval_islands_SANN <- eval_islands_SANN*control$nislands


      if ((control$neliteIslands -  control$nEliteSaved)%%mc.cores == 0) {
        idle_cores <- 0
      } else {
        idle_cores <- (mc.cores - (control$neliteIslands -  control$nEliteSaved)%%mc.cores)
      }
      #idle_cores
      efficiency_islands_SANN <- ((control$neliteIslands -  control$nEliteSaved))/
        ((control$neliteIslands -  control$nEliteSaved)+idle_cores)

      }

      eval_islands <- eval_islands_GA + eval_islands_SANN


    } else {
      efficiency_islands_GA <- 1
      eval_islands_GA <- 0
      efficiency_islands_SANN <- 1
      eval_islands_SANN <- 0
    }


    eval_main_GA <- (ceiling(control$npop/mc.cores) + control$niterations*ceiling((control$npop-control$nelite)/mc.cores))

    if ((control$npop-control$nelite)%%mc.cores == 0) {
      idle_cores <- 0
    } else {
      idle_cores <- (mc.cores - (control$npop-control$nelite)%%mc.cores)
    }
    #idle_cores
    efficiency_main_GA <- ((control$npop-control$nelite))/
      ((control$npop-control$nelite)+idle_cores)

    counter <- 0
    eval_main_SANN <- 0
    efficiency_main_SANN <- 1
    if (control$niterSANN > 0) {
      #if (control$maxiterSANN + control$SANNcooldown >= 1) {
        while (counter <= control$niterations) {
          eval_main_SANN <- eval_main_SANN + control$maxiterSANN*control$niterSANN*ceiling((control$nelite -  control$nEliteSaved)/mc.cores)
          counter <- counter + control$maxiterSANN + control$SANNcooldown
        }
      #}
      #eval_main_SANN <- eval_main_SANN


      if ((control$nelite -  control$nEliteSaved)%%mc.cores == 0) {
        idle_cores <- 0
      } else {
        idle_cores <- (mc.cores - (control$nelite -  control$nEliteSaved)%%mc.cores)
      }
      #idle_cores
      efficiency_main_SANN <- ((control$nelite -  control$nEliteSaved))/
        ((control$nelite -  control$nEliteSaved)+idle_cores)
    }

    eval_main <- eval_main_GA + eval_main_SANN
    eval_total <- eval_main + eval_islands


    average_efficiency <- (efficiency_islands_GA*eval_islands_GA + efficiency_islands_SANN*eval_islands_SANN +
                             efficiency_main_GA*eval_main_GA + efficiency_main_SANN*eval_main_SANN)/eval_total



  } else {


    if (control$nislands%%mc.cores == 0) {
      idle_cores <- 0
    } else {
      idle_cores <- (mc.cores - control$nislands%%mc.cores)
    }
    idle_cores
    efficiency_islands <- (control$nislands)/
      (control$nislands+idle_cores)


    efficiency_main <- 1/mc.cores

    eval_islands <- 0

    if (control$nislands > 1) {
      eval_islands_GA <- ceiling(control$nislands/mc.cores)*(control$npopIslands + control$niterIslands*ceiling(control$npopIslands*0.5))

      #counter <- 0
      eval_islands_SANN <- 0
      # if (control$maxiterSANN + control$SANNcooldown >= 1) {
      #   while (counter <= control$niterIslands) {
      #     eval_islands_SANN <- eval_islands_SANN + control$maxiterSANN*control$niterSANNislands*(control$neliteIslands -  control$nEliteSaved)
      #     counter <- counter + control$maxiterSANN + control$SANNcooldown
      #   }
      # }
      #eval_islands_SANN <- eval_islands_SANN*control$nislands

      eval_islands <- eval_islands_GA + eval_islands_SANN

    }


    eval_main_GA <- (control$npop + control$niterations*ceiling(control$npop*0.5))

    #counter <- 0
    eval_main_SANN <- 0


    eval_main <- eval_main_GA + eval_main_SANN
    eval_total <- eval_main + eval_islands

    average_efficiency <- (efficiency_main*eval_main + efficiency_islands*eval_islands)/eval_total

  }

  if (verbose) {

 {
      write("*************************************************************************", stderr())
      write("TimeEstimation:", stderr())
      write(paste0("- Average time needed to run Stat once: ", time_one_run, " seconds"), stderr())
      write(paste0("- Total number of times Stat is run during optimizaiton: ", eval_total), stderr())
      write(paste0("- Total time needed for optimization process: ", round(eval_total*time_one_run, 6), " seconds\n"), stderr())
      write("Important remarks:", stderr())
      write("- Estimation based on the time required to run Stat. Internal TrainSel
  processes are disregarded. It is accurate unless Stat is extremely fast", stderr())
      write("- It is assumed that convergence is not reached", stderr())
      write("*************************************************************************", stderr())

    }
  }

  if (nStat == 1) {
    output_list <- list(time_Stat_vector = time_one_run_vector,
                        eval_total=eval_total,
                        time_total_seconds=eval_total*time_one_run#,
                       # mc_efficiency_islands_GA=efficiency_islands_GA,
                       # mc_efficiency_islands_SANN=efficiency_islands_SANN,
                       # mc_efficiency_main_GA=efficiency_main_GA,
                       # mc_efficiency_main_SANN=efficiency_main_SANN,
                       # mc_average_efficiency=average_efficiency
                       )
  } else {
    output_list <- list(time_Stat_vector = time_one_run_vector,
                        eval_total=eval_total,
                        time_total_seconds=eval_total*time_one_run#,
                       # mc_efficiency_islands_GA=efficiency_islands,
                       # mc_efficiency_islands_SANN=efficiency_islands,
                       # mc_efficiency_main_GA=efficiency_main,
                       # mc_efficiency_main_SANN=efficiency_main,
                       # mc_average_efficiency=average_efficiency
                       )
  }

  return(output_list)

}


