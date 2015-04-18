compareToMatlab <- function(mat.file,
                            dream.obj){
  library(R.matlab)
  mat <- readMat(mat.file)
  ## CHECK INPUTS
  mat.control <- getMatlabControl(mat)
  cat("\nDid matlab and R use same inputs?")
  cat("\n names in Matlab not in R:")
  cat(setdiff(names(mat.control),names(dream.obj$control)))
  cat("\n names in R not in Matlab: ")
  nn <- setdiff(names(dream.obj$control),names(mat.control))
  nn <- nn[!nn %in% c("Rthres","parallel","burnin.length","maxtime","REPORT")]
  cat(nn)
  common.names <- intersect(names(dream.obj$control),names(mat.control))
  mat.control <- mat.control[common.names]
  r.control <- dream.obj$control[common.names]
  cat("\n identical:",identical(mat.control,r.control))
  cat("\n all.equal:",all.equal(mat.control,r.control))
  eq <- sapply(1:length(r.control),
               function(i) r.control[[i]]==mat.control[[i]])
  names(eq) <- names(r.control)
  cat("\n element-wise:\n")
  print(eq)
  cat("\n")
  ## CHECK OUTPUTS
  ## Compare sequences
  mat.seq <- getMatlabSeq(mat)
  R.seq <- dream.obj$Sequences
  cat("\nDo outputs have same dimensions?\n")
  print(sapply(list(matlab=mat.seq,R=R.seq),function(x) c(
                                     variables=nvar(x),
                                     chains=nchain(x),
                                     iterations=niter(x),
                                     thin=thin(x)
                                     )))
  cut.start=1 + (end(mat.seq) - 1) * (1 - 0.5)
  mat.m <- as.matrix(window(mat.seq,start=cut.start,thin=100))
  R.m <- as.matrix(window(R.seq,start=cut.start,thin=100))
  cat("\n ks.test that samples are from same distribution for each variable\n")
  pvals <- sapply(1:ncol(mat.m),function(i) ks.test(R.m[,i],mat.m[,i])$p.value)
  names(pvals) <- colnames(R.m)
  print(round(pvals,2))
  ## Graphically
  plotMCMCQQ(mat.m,R.m)
}##compareToMatlab

plotMCMCQQ <- function(mat.m,R.m){
  library(reshape)
  library(lattice)
  colnames(mat.m) <- colnames(R.m)
  mm <- rbind(
              cbind(melt(mat.m),source="mat"),
              cbind(melt(R.m),source="R")
              )
  mm2 <- cast(mm,X1+X2~source)
  stopifnot(!any(is.na(mm2)))
  print(xyplot(mat~R|X2,data=mm2,as.table=T,scales='free',
               main="QQ plot of distribution of R and matlab code",
               panel=function(x,y,...){
                 e <- sort(x) ##quantile(x,1:length(x)/length(x))
                 o <- sort(y) ##quantile(y,1:length(y)/length(y))
                 panel.points(e,o)
                 panel.lines(e,e)
               }
               ))
}##plotMCMCQQ

getMatlabSeq <-
  function(mat) {
    NSEQ <- dim(mat$Reduced.Seq)[3]
    NDIM <- dim(mat$Reduced.Seq)[2]-2
    as.mcmc.list(lapply(1:NSEQ,function(i)
                        mcmc(mat$Reduced.Seq[,1:NDIM,i],
                             thin=as.numeric(mat$Extra[dimnames(mat$Extra)[[1]]=="T"])
                             )))
  }

getMatlabControl <- function(mat){
  MCMCPar <- lapply(mat$MCMCPar[,,1],as.vector)
  Extra <- lapply(mat$Extra[,,1],function(x) ifelse(all(dim(x)==c(1,1)),as.vector(x),x))
  list(
       thin.t=Extra$T,
       pCR.Update=Extra$pCR=="Update",
       ndim=MCMCPar$n,
       nseq=MCMCPar$seq,
       ndraw=MCMCPar$ndraw,
       nCR=MCMCPar$nCR,
       gamma=MCMCPar$Gamma,
       DEpairs=MCMCPar$DEpairs,
       steps=MCMCPar$steps,
       eps=MCMCPar$eps,
       outlierTest=MCMCPar$outlierTest,
       boundHandling=tolower(Extra$BoundHandling)
       )
}



writeMatlabDREAMSettings <- function(dream.obj,ModelName,InitPopulation,
                                     matlab.dream.dir,
                                     in.mat.file="in.mat",
                                     out.mat.file="out.mat",
                                     run.m.file="run_once.m"
                                     ){
  library(R.matlab)
  control <- dream.obj$control
  Extra <- list(
                pCR=as.matrix(ifelse(control$pCR.Update,"Update","Fixed")),
                reduced_sample_collection="Yes",
                T=control$thin.t,
                InitPopulation=InitPopulation,
                save_in_memory="No",
                BoundHandling=paste(toupper(substring(control$boundHandling,1,1)),
                  substring(control$boundHandling,2),sep=""),
                DR="No",
                DRscale=10
                )
  MCMCPar <-
    with(control,list(
                      n=as.numeric(ndim),
                      seq=nseq,
                      DEpairs=DEpairs,
                      Gamma=gamma,
                      nCR=nCR,
                      ndraw=ndraw,
                      steps=steps,
                      eps=eps,
                      outlierTest=outlierTest
                      ))
  pars <- eval(dream.obj$call$pars)
  ParRange <- list(
                   minn=matrix(sapply(pars, min),nrow=1),
                   maxn=matrix(sapply(pars, max),nrow=1)
                   )
  Measurement <- list(
                      MeasData=NA,
                      Sigma=NA,
                      N=0
                      )
  func.types <- c("posterior.density","calc.loglik","calc.rmse","logposterior.density","calc.weighted.rmse")
  ## Write and run
  oldwd <- setwd(matlab.dream.dir)
  writeMat(in.mat.file,
           Extra=Extra,
           MCMCPar=MCMCPar,
           ParRange=ParRange,
           Measurement=Measurement,
           ModelName=ModelName,
           option=as.matrix(which(func.types==dream.obj$call$func.type))
           )
  cat(sprintf('
load %s;
[Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option);
save -6 %s
',in.mat.file,out.mat.file),file=run.m.file)
  setwd(oldwd)
  cat(sprintf("
Please ensure that there is a file named '%s.m' in the directory '%s'.
Run '%s' in Matlab and then run readMat('%s') in R to obtain Matlab's output.
",ModelName,matlab.dream.dir,run.m.file,out.mat.file))
} ##writeMatlabDREAMSettings
