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





##############################################

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
           InitSol=NULL) {






########Check types
###Data can only be NULL or list
if ((!is.null(Data) & !is.list(Data))){
  stop("Data is either NULL or a list")
}



####  control is null or  and object of type Trainsel_control
    if ((!is.null(control) &  (class(control)[1]!="TrainSel_Control"))){
      stop("control can be NULL or an integer vector.")
    }


######################
    #########################
  if (is.null(ntotal)) {
      ntotal <- 0
    }

if (!is.function(Stat)){
  if (is.null(Stat))
    {
      CD=TRUE
      Stat=function(){}
  }  else {
    stop('Stat is either NULL (meaning CDMEAN is used) or it is a function')
  }
} else {
  CD=FALSE
}
  if (is.null(Data) & CD){
    stop("no data for MM method")
  }

    if (is.null(control)){
      control=TrainSelControl()
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
out<-TrainSelC(Data=Data, CANDIDATES =Candidates,setsizes =setsizes,settypes=settypes,Stat = Stat,CD=CD,Target=Target,control=control, ntotal=ntotal, InitSol=InitSol)
class(out)<-"TrainSelOut"
} else {
out<-TrainSelCMOO(Data=Data, CANDIDATES =Candidates,setsizes =setsizes,settypes=settypes,Stat = Stat,nstat=nStat, control=control, InitSol=InitSol)
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
           InitSol=NULL) {

    nislands<-control[["nislands"]]
    mc.cores<-control[["mc.cores"]]

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
                InitSol=InitSol)
} else {
  controlInner<-control
  controlInner$progress=FALSE
  solList<-parallel::mclapply(1:nislands, function(x){TrainSelInner(Data = Data,
                                                               Candidates = Candidates,
                                                               setsizes = setsizes,
                                                               ntotal = ntotal ,
                                                               settypes = settypes,
                                                               Stat = Stat,
                                                               nStat= nStat,
                                                               Target = Target,
                                                               control = controlInner)}, mc.cores=mc.cores)
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
                  InitSol=InitSol)

}
return(out)
}


PamSel<-function(M=NULL,distancemethod=c("euclidean",  "maximum", "manhattan", "canberra", "binary", "minkowski"),ntoselect=NULL){
  stopifnot(distancemethod%in%c("euclidean",  "maximum", "manhattan", "canberra", "binary", "minkowski"))
  D=as.matrix(dist(M, method = distancemethod))
  rownames(D)<-colnames(D)<-rownames(M)
  pamx <- cluster::pam(D, k=ntoselect)
  namesselected <- rownames(pamx$medoids)
  return(which(rownames(M)%in%namesselected))
}

#################################TrainSelControl


TrainSelControl <-
  function(npop = 300,
           nelite = 10,
           mutprob = .01,
           mutintensity = 2,
           niterations = 500,
           niterSANN = 200,
           stepSANN = 1e-2,
           minitbefstop = 200,
           tolconv = 1e-7,
           progress=TRUE,
           nislands=1,
           mc.cores=1)
  {
    structure(
      list(
        npop = npop,
        nelite = nelite,
        mutprob = mutprob,
        mutintensity = mutintensity,
        niterations = niterations,
        niterSANN = niterSANN,
        stepSANN = stepSANN,
        minitbefstop = minitbefstop,
        tolconv = tolconv,
        progress=progress,
        nislands=nislands,
        mc.cores=mc.cores
      ),
      class = c("TrainSel_Control")
    )
  }


