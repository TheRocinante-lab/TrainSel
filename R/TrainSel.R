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
TrainSel <-
  function(Data = NULL,
           Candidates = NULL,
           setsizes = NULL,
           ntotal = NULL,
           settypes = NULL,
           Stat = NULL,
           nStat= 1,
           Target = NULL,
           control = NULL) {






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


if (nStat==1){
out<-TrainSelC(Data=Data, CANDIDATES =Candidates,setsizes =setsizes,settypes=settypes,Stat = Stat,CD=CD,Target=Target,control=control, ntotal=ntotal)
class(out)<-"TrainSelOut"
} else {
out<-TrainSelCMOO(Data=Data, CANDIDATES =Candidates,setsizes =setsizes,settypes=settypes,Stat = Stat,nstat=nStat, control=control)
class(out)<-"TrainSelOut"
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
           progress=TRUE)
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
        progress=progress
      ),
      class = c("TrainSel_Control")
    )
  }


