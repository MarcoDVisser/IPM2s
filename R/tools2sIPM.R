######################################################################
## Function to build the IPMS from fit estimates
## and useful function for the analysis 
## This script a set of function to construct and iterate
## IPMS.
## First Build: Dec, 2014, Milwaukee, WI, USA
## Update: 02.02.2015, Gamboa, Panama
## Update: 12.02.2015, Gamboa, Panama
## Update: 04.08.2015, Rheden, Netherlands
## Update: 09.03.2016, Nijmegen, Netherlands
## Update: 25.01.2017, Princeton, NJ
######################################################################
######################################################################

##' Fecundity function 
##'
##' @param dbh diameter at breast height 
##' @param S Seed production per year per unit reproductive basal area
##' @param R Reproductive probabilty (of tree of size dbh)
##' @param F Fraction of crown bearing fruit
##' @param It Iteration for numerical sensitivity calculation.
##' a vector of values, corresponding to parameter dbh and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @import MASS
##' @author Marco D. Visser et al. 
##' @export
fecundity<-function(dbh,S,R,F=1,It=1){
S*R*F*(pi*(dbh/2)^2)*It
}

##' CF: Crown fraction bearing fruit
##' 
##' @param dbh diameter at breast height 
##' @param model list containing the best fitting model
##' and a vector of coefficients
##' @param It Iteration for numerical sensitivity calculation.
##' a vector of values, corresponding to parameter dbh and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
CF<-function(dbh,model, It=1){
  tempmodel <- get(as.character(model[[1]]))
  tempmodel(dbh, betas=model[[2]])*It
}


##' RSC: Reproductive size curve
##' 
##' @param dbh diameter at breast height 
##' @param model list containing the best fitting model
##' and a vector of coefficients
##' @param It Iteration for numerical sensitivity calculation.
##' a vector of values, corresponding to parameter dbh and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
RSC<-function(dbh,model, It=1){
  tempmodel <- get(as.character(model[[1]]))
  tempmodel(dbh, betas=model[[2]])*It
}

##' drecruitsize: density distribution of recruits
##' 
##' @param hght seedling height 
##' @param modelparam a list with the model name
##" and fitted parameters model for initial recruit size
##' @param It Iteration for numerical sensitivity calculation. Also see ?ItMat.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1.
##' @author Marco D. Visser et al. 
##' @export
drecruitsize<-function(hght,modelparam,It=1){
  model <- modelparam[[1]]
  param <- modelparam[[2]]
  if(model=="dexp"){
    return(dexp(hght,param*It))
  }

  if(model=="dlnorm"){
    return(dlnorm(hght,param[1]*It,param[2]))
  }

  else {
    pdf <- get(as.character(model))
    return(pdf(hght,param[1],param[2]*It))
  }
}

##' precruitsize: density distribution of recruits
##' 
##' @param hght seedling height 
##' @param modelparam a list with the model name
##" and fitted parameters model for initial recruit size
##' @param It Iteration for numerical sensitivity calculation.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @param Elt only return the parameter relavent for elasticity
##' @author Marco D. Visser et al. 
##' @export
precruitsize<-function(hght,modelparam,It=1,Elt=FALSE){
   model <- modelparam[[1]]
   param <- modelparam[[2]]
   if(Elt==FALSE){
     
  if(model=="dexp"){
    return(pexp(hght,param*It))
  }

  if(model=="dlnorm"){
    return(plnorm(hght,param[1]*It,param[2]))
  }

  else {
    model<-gsub("d", "p",model)
    pdf <- get(as.character(model))
    return(pdf(hght,param[1],param[2])*It)
    
  } } else {

  if(model=="dexp"){
    return(param)
  }

  if(model=="dlnorm"){
    return(param[1])
  }

  else {
    model<-gsub("d", "p",model)
    return(param[2])

  }}
}
 

##' dSeedlinggrowth distribution of heights t + 1 as a function
##' of initial height
##' 
##' @param hghtref seedling initial height
##' @param hghttar seedling target height
##' @param model fitted model for seedling growth and survival
##' @param gridsize relavent IPM gridsize 
##' @param L lower limit of the IPM. If NULL, seedling to 
##' tree transition are calculated and gridsize must be of 
##' the translated gridsizes according to the height transition
##' model!
##' @param It Iteration for numerical sensitivity calculation.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @param Elt only return the parameter relavent for elasticity
##' @author Marco D. Visser et al. 
##' @export
dSeedlinggrowth<-function(hghtref,hghttar,model=4,param,sigma,
                          gridsize,L,It=1,Elt=FALSE){

  if(!is.null(L)){
    upper <- hghttar+gridsize
    lower <- hghttar
    ## any sizeclass extend lower than the IPM?
    lowest <- lower<=L
    midpoint <- hghtref + (gridsize/2)
  } else {
  hghtref <- hghtref
    lower <- gridsize[[1]]
    upper <- gridsize[[2]]
  lowest <- rep(0,length(lower))
   midpoint <- hghtref
 }

  tempmodel <-get(as.character(model))
  mu <- tempmodel(midpoint,betas=param,time=1)*It
  if(Elt==TRUE) {return(mu)}
  Trs <- pnorm(upper,mu,sigma)-pnorm(lower,mu,sigma)
  ## Integration of lower tail
  IntTrs <- pnorm(lower,mu,sigma)*lowest
#print(list(Trs,IntTrs,upper,lower,hghtref))
  return((Trs+IntTrs)) # cm gridsize
}

##' pSeedlinggrowth: cummulative distribution of growth  
##' transitions beyond a given height
##' 
##' @param hghtref seedling initial height
##' @param hghttar seedling target height
##' @param model fitted model for seedling growth
##' @param It Iteration for numerical sensitivity calculation.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
pSeedlinggrowth<-function(hghtref,hghttar,model,It=1){
hghtref <- hghtref # Joe's data was in cm
hghttar <- hghttar # Joe's data was in cm

betas <- fixef(model)
sigma <- sd(residuals(model))
mu <- hghtref+betas[1]+betas[2]*hghtref
return(pnorm(hghttar,mu*It,sigma))
}

##' Seedlingsurvival: survival of Seedlings as a function
##' of initial height
##' 
##' @param hght Seedling initial height
##' @param model list containing the fitted model function
##' and a set of parameters for the model
##' @param gridsize IPM mesh size
##' @param It Iteration for numerical sensitivity calculation.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
Seedlingsurvival<-function(hght,model,gridsize,L,It=1){

  if(!is.null(L)){
    midpoint <- hght + (gridsize/2)
  } else {
   midpoint <- hght
 }

  tempmodel <-get(as.character(model[[1]]))
  tempmodel(midpoint,betas=model[[2]])*It

}

##' sdl2splg: Seedling to sapling transition.
##' height to dbh translation
##' @param hght hght to dbh
##' @param dbh dbh to height, ignored if NULL
##' @param model fitted parameters model for tree-height relations
##' @param sp species 6 letter code in capitals
##' @author Marco D. Visser et al. 
##' @export
sdl2splg <- function(hght,dbh=NULL,model,sp){
 B0 <- model$intercept[model$sp==sp]
 B1 <- model$slope[model$sp==sp]
   
  if(is.null(dbh)){
    return(as.numeric(B0+B1*hght))
  } else {
    return(as.numeric((dbh-B0)/B1))
  }
}

##' dTreegrowth: distribution of dbhs t + 1 as a function
##' of initial dbh
##' 
##' @param dbhref tree initial dbh
##' @param dbhtar tree target dbh
##' @param model fitted model for tree growth and survival
##' @param gridsize relavent IPM gridsize 
##' @param L lower limit of the IPM
##' @param It Iteration for numerical sensitivity calculation.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @param Elt only return the parameter relavent for elasticity
##' @author Marco D. Visser et al. 
##' @export
dTreegrowth<-function(dbhref,dbhtar,model,param,sigma,gridsize,L,U,It=1,
                      Elt=FALSE){

    upper <- dbhtar+(gridsize)
    lower <- dbhtar
    mid <- dbhref+(gridsize/2)
    ## Largest and smallest class id
    lowest <- dbhtar<=L
    largest <- dbhtar>=U
        
    ## Predictions
    tempmodel <- get(as.character(model))
    mu <- tempmodel(mid,time=1,betas=param)*It
    if(Elt==TRUE) {return(mu)}
    Trs <- pnorm(upper,mu,sigma)-pnorm(lower,mu,sigma)
    ## Integration of lower tail
    IntTrs <- pnorm(L,mu,sigma)*lowest 
    ## Integration of the upper tail
    UpIntTrs <- (1-pnorm(U+gridsize,mu,sigma))*largest
    
    ##print(list(Trs,IntTrs,upper,lower,dbhref))
    return((Trs+IntTrs+UpIntTrs)) # cm gridsize
}


##' Treesurvival: survival of trees as a function
##' of initial dbh
##' 
##' @param dbh tree initial dbh
##' @param model fitted survival model
##' with first number the best model
##' from TSmodelnames
##' @param gridsize IPM mesh size for dbh
##' @param It Iteration for numerical sensitivity calculation.
##' A vector of values, corresponding to parameter dbh/hght and of the same
##' length, with with to interate the vital rate function. Defaults to 1,
##' which means no iteration will be done. The function ItMat can be
##' use to build iterations. see ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
Treesurvival<-function(dbh,model,gridsize,It=1){
  midpoint <- dbh + (gridsize/2)
  tempmodel <-get(as.character(model[[1]]))
  (tempmodel(midpoint,betas=model[[2]],onlyfixedeffect=FALSE)*It)^(1/9.8)
}

######################################################################
#IPM construction tools
######################################################################

##' Iteration matrix builder
##'
##' Constructs and iteraction matrix used in numerical
##' sensitivity calculations. Inputs are the vital rate,
##' stage and size that should be iterated.
##'
##' @param vitalrate a character specifying the vital rate
##' that should be iterated can be c("seed production",
##'  "reproduction", "recruit size", "seedling growth",
##'  "seedling survival", "tree growth",
##'  "tree survival")
##' @param delta the iteration factor, defualts to 1.01
##' @param size at which size must the iteration take place?
##' this is expected to be an integer corresponding to
##' a position in seedling or treeclasses.  
##' @param seedlingclasses the (discretized) seedling classes
##' which are input into the IPM function and kernel functions
##' @param treeclasses the (discretized0 tree classes
##' which are input into the IPM function and kernel functions
##' @author Marco D. Visser et al. 
##' @export
ItMat <- function(vitalrate="tree growth",delta=1.01,size,
                    seedlingclasses,treeclasses){

  vrates <- c( "seedling growth",
              "seedling survival", 
               "seed production",
              "reproduction",
              "recruit size",
              "tree growth",
              "tree survival")

  Smat <- matrix(1,ncol=2,nrow=length(seedlingclasses))
  Tmat <- matrix(1,ncol=5,nrow=length(treeclasses))

  if(grep(vitalrate,vrates)<3){
    Smat[size,grep(vitalrate,vrates)] <- delta
  } else{
    Tmat[size,grep(vitalrate,vrates)-2] <- delta
  }

return(list(Smat,Tmat))

}


##' Rkernel: Recruits kernel
##' 
##' Function returns recruitment (N individuals) between 
##' dbhref and hghttar 
##' @param dbhref mother tree diameter
##' @param hghttar recruit hght
##' @param models list of models or parameters for seed 
##' to seedling transition, reproduction, seed production, 
##' and inititail recruit height
##' @param gridsize gridsize of the seedling IPM
##' @param smodels Tree survival models
##' @param tImat Tree iteration matrix, which in this case
##' is expected to be a vector with values for seed production, reproduction,
##' recruitsize  tree growth &  tree survival
##' @author Marco D. Visser et al. 
##' @export
Rkernel<-function(dbhref,hghttar,models,gridsize,treegridsize,smodels,
                    tImat){
  Sv <- Treesurvival(dbh=dbhref,model=smodels,
                       gridsize=treegridsize,It=tImat[5])

  Sv*models[[1]]*fecundity(dbhref+(treegridsize/2),
                           R=RSC(dbhref+(treegridsize/2),models[[2]]
                             ,It=tImat[2]),
                           F=CF(dbhref+(treegridsize/2),models[[2]][[3]]
                             ,It=tImat[2]),
                           S=models[[3]],It=tImat[1])*
  (precruitsize(hghttar+(gridsize),modelparam=models[[4]],It=tImat[3])-
   precruitsize(hghttar,modelparam=models[[4]],It=tImat[3])
   )
}

##' Sdlgskernel: Seedling growth and survival transitions
##' 
##' Function returns transition probabilities between 
##' hghtref and hghttar or dbh target
##' @param hghtref seedling initial hght,reference height
##' @param hghttar target height
##' @param dbhhghtthres height at 1 cm dbh
##' @param models list of models for seedling growth and survival
##' @param gridsize seedling matrix grid size
##' @param L lower limit of the IPM
##' @param sImat A seedling iteration matrix. See ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
sdlGSkernel<-function(hghtref,hghttar,models=NULL, gridsize=10,L=0.01,sImat){
  Sv <- Seedlingsurvival(hght=hghtref,model=models[[1]],gridsize=gridsize,L=L,
                         It=sImat[,2])
Gr <- dSeedlinggrowth(hghtref,hghttar,model=models[[2]][[1]],
                      param=models[[2]][[2]],sigma=models[[2]][[3]],gridsize,L,
                      It=sImat[,1])
return((Sv*Gr))
}

##' treeTranskernel: Seedling to tree transition matrix
##' 
##' Function returns transition probabilities between 
##' a referebce height (hghtref) and dbh target (dbhtar)
##' @param hghtref seedling initial hght,reference height
##' @param dbhtar target dbh
##' @param dbhhghtthres height at 1 cm dbh
##' @param models list of models for dbh height translation, 
##' seedling growth and survival
##' @param gridsize seedling matrix grid size (for survival)
##' @param sp 6 letter species code
##' @param sImat A seedling iteration matrix. See ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
treeTranskernel<-function(hghtref,dbhtar,models=NULL,gridsize.tree,
                          gridsize.sdl,sp=NULL,sImat){

    ## translate dbh to hght
    hghttar <-  sdl2splg(dbh=dbhtar,model=models[[1]],sp=sp)
    Lower <- hghttar
    Upper <- sdl2splg(dbh=dbhtar+gridsize.tree,model=models[[1]],sp=sp)


    Sv <- Seedlingsurvival(hght=hghtref,model=models[[2]],
                       gridsize=gridsize.sdl,L=1,It=sImat[,2])
    Gr <- dSeedlinggrowth(hghtref,hghttar,model=models[[3]][[1]],
                      param=models[[3]][[2]],sigma=models[[3]][[3]],
                      gridsize=list(Lower,Upper),
                      L=NULL,It=sImat[,1])
return((Sv*Gr))
}


##' Tree growth and survival kernel
##' 
##' Function returns transition probabilities between 
##' dbhref and dbhtar
##' @param dbhref tree initial dbh,reference dbh
##' @param dbhtar target dbh
##' @param models list of models for tree growth and survival
##' @param gridsize seedling matrix grid size
##' @param L lower limit of the IP
##' @param tImat A tree iteration matrix. See ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
treeGSkernel<-function(dbhref,dbhtar,models,gridsize,
                       L=10,tImat,U=max(treeclasses)){
Sv <- Treesurvival(dbhref,model=models[[1]],gridsize=gridsize,It=tImat[,5])
Gr <- dTreegrowth(dbhref,dbhtar,model=models[[2]][[1]],
                  param=models[[2]][[2]],sigma=models[[2]][[3]],gridsize,L,U,
                  It=tImat[,4])
return((Sv*Gr))
}




##' Construct IPM
##'
##' Function integrates all previous functions
##' and construct a full IPM
##' @param seedlingclasses vector of seedling sizes
##' @param treeclasses vector of treeseedling sizes
##' @param seedlingclasseswidth seedling grid size
##' @param treeclasseswidth tree grid size
##' @param sizerange min and max sizes
##' @param hghtdbhthreshold seedling to tree threshold size in height
##' @param dbhthreshold seedling to tree threshold size in dbh
##' @param f.sp seed production model parameters
##' @param f.sts seed 2 seedling parameters
##' @param f.repro size dependent reproduction model
##' @param f.sdlhght initial recruit height model
##' @param f.s2t allometric function
##' @param f.sg seedling growth model
##' @param f.ss seedling survival model
##' @param f.tg tree growth model
##' @param f.ss tree survival model
##' @param Imat A iteration matrix used for numerical sensitivity calculations.
##' See ?ItMat.
##' @author Marco D. Visser et al. 
##' @export
constructIPM <- function(seedlingclasses,treeclasses,
                         seedlingclasseswidth,treeclasseswidth,
                         sizerange,hghtdbhthreshold,dbhthreshold,
                         f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                         f.ss,f.tg,f.ts,
                         Imat = Imat) {

  sImat <- Imat[[1]]
  tImat <- Imat[[2]]

  
## Reproduction
######################################################################
# Reproduction kernel model list

  rMods <- list(f.sts$rrate,f.repro,f.sp$fec,
                 list(f.sdlhght$bestmod,
                      na.omit(as.numeric(unlist(
                                f.sdlhght[,grep("p_",names(f.sdlhght))])))
 
                    )
                )

  tMods <- list(f.ts[[2]],f.ts[[1]])
  
                                        # make recruitment kernel
  rIPM <- matrix(rep(treeclasses,each=length(seedlingclasses)),
                 ncol=length(treeclasses),
                 nrow=length(seedlingclasses),byrow=F)
  rIPM <- sapply(seq_along(rIPM[1,]),function(X) 
                 Rkernel(rIPM[,X],seedlingclasses, rMods,
                         gridsize=seedlingclasswidth,
                        treegridsize=treeclasswidth,tMods,tImat=tImat[X,]))

                                        #zero recruitment kernel
                                        #used in various calculations
  zrIPM <- matrix(0, ncol=length(treeclasses),
                 nrow=length(seedlingclasses),byrow=F)

  
  ## Seedlings
######################################################################
  ## Seedling growth and survival model list
  prepsdlgrowmodel <- list(f.sg[[2]],f.sg[[1]],f.sg[[1]][[3]])
  prepsdlsurvmodel <- list(f.ss[[2]],f.ss[[1]])
  
  sMods <- list(prepsdlsurvmodel,prepsdlgrowmodel)


                                        # make seedling kernel
  sIPM <- matrix(rep(seedlingclasses,each=length(seedlingclasses)),
                 ncol=length(seedlingclasses),
                 nrow=length(seedlingclasses),byrow=F)
  sIPM <- apply(sIPM,2,function(X) 
                sdlGSkernel(X,seedlingclasses,
                            sMods,gridsize=seedlingclasswidth,
                            L=sizerange[1],sImat)) # apply returns columns

  ## Seedling to Tree
######################################################################
                                        # Tree transition  model list

  s2tMods <- list(f.s2t,prepsdlsurvmodel,prepsdlgrowmodel)

                                        # make seedling to tree transition kernel
  s2tIPM <- matrix(rep(seedlingclasses,each=length(treeclasses)),
                   ncol=length(seedlingclasses),
                   nrow=length(treeclasses),byrow=F)

  s2tIPM <- sapply(seq_along(s2tIPM[1,]),function(X) 
                  treeTranskernel(s2tIPM[,X],treeclasses,
                                  s2tMods,gridsize.tree=treeclasswidth,
                                  gridsize.sdl=seedlingclasswidth,
                                  sp=SpFitList[i],sImat=
                                  t(matrix(sImat[X,])))) # apply returns columns

  ## Trees
######################################################################
  ## Tree growth and surival list

  prepTrgrowmodel <-  list(f.tg[[2]],na.omit(f.tg[[1]]),f.tg[[1]][[3]])
  prepTrsurvmodel <- list(f.ts[[2]],f.ts[[1]])
  
  tMods <- list(prepTrsurvmodel,prepTrgrowmodel)

                                        # make  tree growth and survival kernel
  tIPM <- matrix(rep(treeclasses,each=length(treeclasses)),
                 ncol=length(treeclasses),
                 nrow=length(treeclasses),byrow=F)

  tIPM <- sapply(seq_along(tIPM[1,]),function(X)  
                 treeGSkernel(tIPM[,X],treeclasses, tMods,
                              gridsize=treeclasswidth,
                             L=dbhthreshold,tImat=
                             t(matrix(tImat[X,])))) # apply returns columns


## Full Tree IPM
######################################################################
                                        #start full IPM construction
  fulldim <- length(treeclasses)+length(seedlingclasses)
  IPM<-matrix(0,ncol=fulldim,nrow=fulldim)

                                        # top right
  IPM[1:length(seedlingclasses),1:length(seedlingclasses)] <- sIPM
                                        # top left
  IPM[1:length(seedlingclasses),(length(seedlingclasses)+1):fulldim] <- rIPM
                                        # bottom right
  IPM[(length(seedlingclasses)+1):fulldim,1:length(seedlingclasses)] <- s2tIPM
                                        # bottom left
  IPM[(length(seedlingclasses)+1):fulldim,
      (length(seedlingclasses)+1):fulldim] <- tIPM

  ## Growth and survival matrix 
  GSM <- IPM
  GSM[1:length(seedlingclasses),(length(seedlingclasses)+1):fulldim] <- zrIPM
  FM <- IPM - GSM

  return(list(IPM=IPM,GSM=GSM,FM=FM,rIPM=rIPM,sIPM=sIPM,s2tIPM=s2tIPM,tIPM=tIPM,
              tMods=tMods,s2tMods=s2tMods,
                sMods=sMods,rMods=rMods))
  
}


require(MASS)

##' SEMatrix: Calculate Sensitivity or Elasticity Matrix
##' 
##' Function returns elasticity or sensitivity
##' @param IPM the IPM. As square matrix
##' @param type either elasticity or sensitivity
##' @export
SEMatrix <- function(IPM,type="elasticity"){
  ## right and left eigen decomposition
  Eig <- eigen(IPM)
  w <- Re(Eig$vectors[,1])
  v <- Re(eigen(t(IPM))$vectors[,1])

  ## stablestable and reproductive values matrices
  Vmat <- apply(IPM,1,function(X) v)
  Umat <- apply(IPM,2,function(X) w)
  
  ## Sensitivity Matrix
  scalarprod <- sum(v*w)

  Smat	<-outer(v,w)/scalarprod
  if(type=="sensitivity"){
    return(Smat)} else {
      return((IPM*Smat)/Re(Eig$values[1]))
    }
}

# parameters - an IPM (with full survival and fecundity complement)
# returns - the sensitivity of every transition for pop growth
sens<-function(A) {
	w<-Re(eigen(A)$vectors[,1]); 
	v<-Re(eigen(t(A))$vectors[,1]);
	vw<-sum(v*w);
	s<-outer(v,w)
	return(s/vw); 
}   

#parameters - an IPM (with full survival and fecundity complement)
# returns - the elasticity of every transition for pop growth
elas<-function(A) {
	s<-sens(A)
	lam<-Re(eigen(A)$values[1]);
	return((s*A)/lam);
}


##' Transpose IPM for plotting
##' Function flips the IPM for correct plotting
##'
##' @param ipm the square matrix IPM
##' @export 
tI <- function(ipm){
  limz <- dim(ipm)
  t(ipm)[,limz[2]:1]
}

##' Calculate Mean Life Expectancy
##' Returns the mean life expectancy for each size
##' @param GSM A growth-survival matrix
##' @export 
LifeExp <- function(GSM){
## get life expectancy
lifeExpect <- colSums(ginv(diag(nrow(GSM)) - GSM))
}

##' Calculate mean passage time to a certain size
##' Returns the mean time from class to each size
##' @param class the row/column number from which to calculate
##' passage time
##' @param GSM A growth-survival matrix
##' @export 
PassTime <- function(class,GSM){
  dims <- dim(GSM)[1]
  Tprime <- GSM
  Tprime[, class] <- 0
  Mprime <- 1 - colSums(GSM)
  Mprime[class] <- 0
  Mprime <- rbind(Mprime, rep(0, dims))
  Mprime[2, class] <- 1
  Bprime <- Mprime %*% ginv(diag(dims) - Tprime)
  Bprime[2, ][Bprime[2, ] == 0] <- 1
  diagBprime <- diag(Bprime[2, ])
  Tc <- diagBprime %*% Tprime %*% ginv(diagBprime)
  eta1 <- ginv(diag(dims) - Tc)
  time.to.absorb <- colSums(eta1)
  time.to.absorb[class:length(time.to.absorb)] <- 0
  return(time.to.absorb)
}

##' Calculate Net Reproductive Rate 
##' The net reproduction rate (NRR) is the average number of daughters
##' that would be born to a female (or a group of females) if she
##' passed through her lifetime conforming to the age-specific fertility
##" and mortality rates of a given year.
##' @param GSM A growth-survival matrix
##' @param FM A fertility matrix
##' @export 
NetR0 <- function(GSM,FM){
  s <- length(diag(GSM))
   N <- try(solve(diag(s) - GSM), silent=TRUE)
   if(class(N)=="try-error"){r<-NA}
   else{
     R <- FM %*% N
     r<-lambda(R)
   }
   r
}

##' Calculate Generation time 
##' Returns the mean time between parent birth and offspring birth
##' @param GSM A growth-survival matrix
##' @param FM A fertility matrix
GenTime <- function(GSM,FM){
  s <- length(diag(GSM))
   N <- try(solve(diag(s) - GSM), silent=TRUE)
  if(class(N)=="try-error"){generation.time<-NA}
   else{
     R <- FM %*% N
     Ro<- lambda(R)
   lambda <- lambda(GSM+FM) 
   generation.time = log(Ro)/log(lambda)   
  }
   generation.time


}


##' Retrive the dominant eigenvalue
##' find the dominant eigenvalue and returns it
##' function shamelessly stolen from the PopBio package
##' by Stubben et al. Which you should really be using
##' instead of this!
##' @param A a square projection matrix
##' @param ... Additional parameters passed to eigen
##' @export 
lambda<-function(A, ...)
{
    ev <- eigen(A,only.values=TRUE, ...)
    # R sorts eigenvalues in decreasing order, according to Mod(values)
    #  ususally dominant eigenvalue is first (ev$values[1]), except for
    #  imprimitive matrices with d eigenvalues of equal modulus
    # this should work for most cases
    lmax <- which.max(Re(ev$values))
    lambda <- Re(ev$values[lmax])
    lambda
}


##' Numerically Iterate size-specific vital rates and retrive
##' the dominant eigenvalue.
##' @param vitalrate which vital rate to iterate? See ?ItMat.
##' @param delta the factor of iteration
##' All parameters are the same as in the constructIPM function.
##' @export 
NumericIterate <- function(seedlingclasses,treeclasses,
                     seedlingclasseswidth,treeclasseswidth,
                     sizerange,hghtdbhthreshold,dbhthreshold,
                     f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                            f.ss,f.tg,f.ts,vitalrate="tree growth",
                            delta=1.01) {

     vrates <- c( "seedling growth",
              "seedling survival", 
               "seed production",
              "reproduction", "recruit size",
              "tree growth",
              "tree survival")

     if(grep(vitalrate,vrates)<3){
       sizes <- seq_along(seedlingclasses)
     } else {
       sizes <- seq_along(treeclasses)
           }

     lambdalist <- vector("list",seq_len(sizes))
     
     for(i in sizes){
       
       ## No elasticity calculations needed for this step
       
       Imat <- ItMat(vitalrate=vitalrate,delta=delta,size=i,
                                      seedlingclasses,treeclasses)
       
       ## Build the IPM
       fullIPM <- constructIPM(seedlingclasses,treeclasses,
                               seedlingclasseswidth,treeclasseswidth,
                               sizerange,hghtdbhthreshold,dbhthreshold,
                               f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                               f.ss,f.tg,f.ts,Imat=Imat)

       ## Imp list, transition kernel, fertility matrix
       lambdalist[[i]] <- lambda(fullIPM[["IPM"]])


       
     }

     
    ## collect models from final ipm abd calculate vital rates
  Rmods <-  fullIPM[["rMods"]]
  SdlMods <- fullIPM[["sMods"]]
  s2tMods <- fullIPM[["s2tMods"]]
  tMods <- fullIPM[["tMods"]]

     ## calculate vital rates
   SvSdl <- Seedlingsurvival(hght=seedlingclasses,model=SdlMods[[1]]
                                   ,gridsize=seedlingclasswidth,L=seedlingclasses[1])
   SvTr <- Treesurvival(dbh=treeclasses,model=tMods[[1]],
                              gridsize=treeclasswidth)

   F <-  fecundity(treeclasses+treegrid,R=RSC(treeclasses+treegrid,
                                          Rmods[[2]]),S=Rmods[[3]])
   R<- RSC(treeclasses+treegrid,Rmods[[2]])
  
     dR <-precruitsize(seedlingclasses+(sdlgrid),modelparam=Rmods[[4]],
                       Elt=TRUE)

   GrSdl <- dSeedlinggrowth(seedlingclasses,seedlingclasses,
                     model=SdlMods[[2]][[1]],
                            param=SdlMods[[2]][[2]],sigma=SdlMods[[2]][[3]],
                            seedlingclasswidth,seedlingclasses[1],
                      It=1,Elt=TRUE)

   GrTr <- dTreegrowth(treeclasses,treeclasses,model=tMods[[2]][[1]],
                     param=tMods[[2]][[2]],sigma=tMods[[2]][[3]],
                     treeclasswidth,L=10,U=9000,
                  It=1,Elt=TRUE)

     vitalrates <- list(GrSdl,SvSdl,F,R,dR,GrTr,SvTr)

                     
     names(vitalrates) <- c( "seedling growth",
              "seedling survival", 
               "seed production",
              "reproduction", "recruit size",
              "tree growth",
                            "tree survival")
                     
return(list(lambdas=unlist(lambdalist),vitalrates=vitalrates[[vitalrate]]))

}

##' Calculate  elasticty values from a set of lambdas returned by NumericIterate
##' 
##' With the original lambda and a list of iterated lamdas, elasticities
##' are calculated
##' @param lambda the original eigen value
##' @param lambdalist a list of numerically iterated lambdas and corresponding
##' vital rates returned by NumericIterate.
##' @param delta the iteration factor
##' @param vrlist list with the original size specific vital rates
##' @export 
NumericElast <- function(lambdalist,lambda,delta=1.01){
  lambdas <- lambdalist$lambdas
  vrates <- lambdalist$vitalrates
  sens = abs(lambdas-lambda)/abs(vrates*abs(delta-1))
  elasticity = sens * (vrates/lambda)
 ## elasticity = abs(lambdas - lambda) / lambda / abs(delta-1)
return(elasticity)
}


##'  lmepredict: function to predict vital rates
##'  conform to the standard approach, but with
##'  more control on some processes.
##' LME predict function
##' @param x design matrix
##' @param liana which liana class?
##' @param mod fit lme4 model and normalization coef (in a list)
##' @param sp which species to use
##' @param beta used to provide fit with legacy code (mod = beta)
##' @param time time correction for binomial models?
##' @param onlyfixedeffects logical: only used the fixed effects in predictions
##' @param switchoffliana numeric: which liana cat (1:4) to switch off?
##' @param returngrowth return growth instead of size next?
##' 
##' @export 
lmepredict <- function(x,mod=NULL,liana=NULL,sp=NULL,betas=NULL,time=1,
                       onlyfixedeffect=FALSE,switchoffliana=NULL,
                       returngrowth=FALSE){


  
    if(onlyfixedeffect==TRUE){
    return(fixedpredict(x,mod=mod,liana=liana,sp=sp,betas=betas,
                        time=time))

    }
  if(is.null(mod)){ mod <- betas}

  xorig <- x
  
  if(is.null(liana)){ liana <- get("GlobalLianaLevel",pos=globalenv())}
  
  if(!is.null(liana)){tmp <- rep(0,4)
                      tmp[liana] <- 1
                      liana <- tmp
                    } else {
                      liana <- rep(0,4)
                    }

  if(!is.null(switchoffliana)){liana[switchoffliana] <- 0}
  
  if(is.null(sp)){ sp <- get("foc.sp",pos=globalenv())}
  
  if(is.list(mod)){
    norms <- mod[[2]]
    mod <- mod[[1]]
    x <- (x-norms[1])/norms[2]
  }
  
  modcall <- mod@call
  
  if(grepl("liana", modcall[[2]][3])){

    ## additive model or not?
    if(grepl("* as.factor", modcall[[2]][3])){

    fe <- fixef(mod)
    re <- ranef(mod)$sp
    spn <- which(rownames(re)==sp)

    ## build model frames
    n <- length(x)
    i.frame <- c(1,liana)

    s.frame <- c(1,liana)
    
    ## build theta vectors
    ## Random slopes, intercepts or both?
    if(grepl("*as.factor", modcall[[2]][3])){
      ## random intercepts + slopes
      int.theta <- sum(as.numeric(c(fe[c(1,3:6)]+re[spn,c(1,3:6)]))*i.frame)
      slope.theta <- sum(as.numeric(c(fe[2]+re[spn,2],fe[-c(1:6)]))*s.frame)
    }
    if(grepl(":as.factor", modcall[[2]][3])){
    ## random slopes
      int.theta <- sum(as.numeric(c(fe[c(1,3:6)]+re[spn,1]))*i.frame)
      ## get global liana level
      gll <- get("GlobalLianaLevel",pos=globalenv())
      ## special frame random slopes (LCC=0:4)
      rs.frame <- ifelse(gll==0,c(1,liana),c(0,liana))
      slope.theta <- sum(as.numeric(c(fe[2]+re[spn,2],fe[-c(1:6)]))*s.frame,
                         as.numeric(re[spn,-c(1,2)])*rs.frame)
    } else {
      ## random intercepts
      int.theta <- sum(as.numeric(c(fe[c(1,3:6)]+re[spn,c(1,3:6)]))*i.frame)
      slope.theta <- sum(as.numeric(c(fe[2]+re[spn,2],fe[-c(1:6)]))*s.frame)
    }
    
      ## Predictions
      pred <- int.theta+slope.theta*x
    
    if(mod@call[4]=="binomial()") {
      pred <- (plogis(pred))^time
    }

  } else {

    fe <- fixef(mod)
    re <- ranef(mod)$sp
    spn <- which(rownames(re)==sp)

    ## build model frames
    n <- length(x)
    i.frame <- c(1,liana)

    s.frame <- c(1)
    
    ## build theta vectors
    int.theta <- sum(as.numeric(c(fe[c(1,3:6)]+re[spn,c(1,3:6)]))*i.frame)
    slope.theta <- sum(as.numeric(c(fe[2]+re[spn,2])*s.frame))
    pred <- int.theta+slope.theta*x
    
    if(mod@call[4]=="binomial()") {
      pred <- (plogis(pred))^time
    }
  }
    
  } else {

    fe <- fixef(mod)
    sizef <- fe[grep("size",names(fe))]
    re <- ranef(mod)$sp
    re <- as.numeric(re[rownames(re)==sp,])
    beta0 <- fe[grep("Intercept",names(fe))]
    pred <- (beta0+re[1])+(sizef+re[2])*x
    if(mod@call[4]=="binomial()") {
      pred <- plogis(pred)
    } 
  }

  if(modcall[[2]][2]=="grw()"){
   pred <- xorig+pred
 }
  if(modcall[[2]][2]=="baagr()"&returngrowth==FALSE){
   x <- xorig
   pred <- (pi*(xorig/2)^2)+pred
   pred <- ifelse(pred<=0,10,pred)
   pred <- sqrt(pred/pi)*2
  }

  return(pred) 
}

##' LME predict function (only fixed effects)
##' @param x design matrix
##' @param liana which liana class?
##' @param mod fit lme4 model and normalization coef (in a list)
##' @param sp which species to use
##' @export
fixedpredict <- function(x,mod=NULL,liana=NULL,sp=NULL,betas=NULL,time=1){


  if(is.null(mod)){ mod <- betas}

  xorig <- x
  
  if(is.null(liana)){ liana <- get("GlobalLianaLevel",pos=globalenv())}
  
  if(!is.null(liana)){tmp <- rep(0,4)
                      tmp[liana] <- 1
                      liana <- tmp
                    } else {
                      liana <- rep(0,4)
                    }
  if(is.null(sp)){ sp <- get("foc.sp",pos=globalenv())}
  
  if(is.list(mod)){
    norms <- mod[[2]]
    mod <- mod[[1]]
    x <- (x-norms[1])/norms[2]
  }
  
  modcall <- mod@call
  
  if(grepl("liana", modcall[[2]][3])){

    ## additive model or not?
    if(grepl("* as.factor", modcall[[2]][3])){

    fe <- fixef(mod)

    ## build model frames
    n <- length(x)
    i.frame <- c(1,liana)

    s.frame <- c(1,liana)
    
    ## build theta vectors
    ## Random slopes, intercepts or both?
    if(grepl("*as.factor", modcall[[2]][3])){
      ## random intercepts + slopes
      int.theta <- sum(as.numeric(fe[c(1,3:6)])*i.frame)
      slope.theta <- sum(as.numeric(c(fe[2],fe[-c(1:6)]))*s.frame)
    }
    if(grepl(":as.factor", modcall[[2]][3])){
    ## random slopes
      int.theta <- sum(as.numeric(c(fe[c(1,3:6)]))*i.frame)
      ## get global liana level
      gll <- get("GlobalLianaLevel",pos=globalenv())
      slope.theta <- sum(as.numeric(c(fe[2],fe[-c(1:6)]))*s.frame)
    
    } else {
      ## random intercepts
      int.theta <- sum(as.numeric(c(fe[c(1,3:6)]))*i.frame)
      slope.theta <- sum(as.numeric(c(fe[2],fe[-c(1:6)]))*s.frame)
    }
    
      ## Predictions
      pred <- int.theta+slope.theta*x
    
    if(mod@call[4]=="binomial()") {
      pred <- (plogis(pred))^time
    }

  } else {

    fe <- fixef(mod)
    
        ## build model frames
    n <- length(x)
    i.frame <- c(1,liana)

    s.frame <- c(1)
    
    ## build theta vectors
    int.theta <- sum(as.numeric(c(fe[c(1,3:6)])*i.frame))
    slope.theta <- sum(as.numeric(fe[2]*s.frame))
    pred <- int.theta+slope.theta*x
    
    if(mod@call[4]=="binomial()") {
      pred <- (plogis(pred))^time
    }
  }
    
  } else {

    fe <- fixef(mod)
    sizef <- fe[grep("size",names(fe))]

    beta0 <- fe[grep("Intercept",names(fe))]
    pred <- (beta0)+(sizef)*x
    if(mod@call[4]=="binomial()") {
      pred <- plogis(pred)
    } 
  }

  if(modcall[[2]][2]=="grw()"){
   pred <- xorig+pred
  }
  if(modcall[[2]][2]=="baagr()"){
   x <- xorig
   pred <- (pi*(xorig/2)^2)+pred
   pred <- ifelse(pred<=0,10,pred)
   pred <- sqrt(pred/pi)*2
  }

  return(pred) 
}



  
##' linearityCheck
##'
##' @param lmeMod mixed model 
##' @param xname name of explanatory variable along which linearity
##'  needs to be checked
##' @import gam
##' @import lme4
##' @author Marco D. Visser et al. 
##' @export
linearityCheck <- function(lmeMod,xname){

    y <- residuals(lmeMod)
    x <- lmeMod@frame[,xname]
    xaxz<-seq(range(x)[1],range(x)[2],length.out=15)
    xaxz2<-seq(range(x)[1],range(x)[2],length.out=100)

    Mod1<-gam(y~lo(x))

    axlabs <- xaxz
    
    plot(x,y,pch=".",main="Residual plot",xaxt="n",
         ylab='residual growth',xlab="size")


    axis(1,at=xaxz,labels=round(axlabs,2))

    predz<-predict(Mod1,newdata=data.frame(x=xaxz2))
    se<-sd(y)

    lines(xaxz2,predz,col='red',lty=3,lwd=4)
    lines(xaxz2,predz+se,col='green',lty=1,lwd=3)
    lines(xaxz2,predz-se,col='green',lty=1,lwd=3)
    abline(h=0,lty=2,lwd=2,col=rgb(0,0,0,alpha=0.6))

    legend("bottom",legend=c("residuals","Moving average", "CI"),
           col=c("black","red","green"),pch=18,lty=c(3,3,1))

}
