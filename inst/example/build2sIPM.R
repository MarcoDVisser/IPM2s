################################################################################
################################################################################
##  Example script for building IPMs for a given species list
## (object SpFitList) from fitted models
######################################################################
## Update 20.01.2017, Princeton, NJ

## Load all fit model objects
require(lme4)
require(IPM2s)

################################################################################
################################################################################
##  Step 1: Simulate datasets
################################################################################
################################################################################
## Species
SpFitList <- c("A","B","C")

################################################################################
## < 1 cm dbh Individuals = Seedlings
################################################################################

## Number of seedlings per species
N <- 1000
Seedlings <- data.frame(
    id=paste0(rep(SpFitList,each=N),rep(1:N,length(SpFitList))),
    height=rep(seq(1,2000,length.out=N),length(SpFitList)),
    sp=rep(SpFitList,each=N))

################################################################################
## 1 cm dbh + Individuals = Trees
################################################################################
cat("\t ##################################################### \n ")
cat("\t \t \t Creating example datasets \n")
cat("\t ##################################################### \n ")

## Number of trees per species
N <- 1000
Trees <- data.frame(
    id=paste0(rep(SpFitList,each=N),rep(1:N,length(SpFitList))),
    dbh=rep(1:N,length(SpFitList)),
    sp=rep(SpFitList,each=N))


## Survival & growth per species
sB0 <- rep(c(1.3,1.6,1.8),each=N) # Seedlings
sB1 <- rep(c(0.001,0.002,0.003),each=N) # Seedlings
B0 <- rep(c(1.7,2.2,3),each=N) # Trees
B1 <- rep(c(0.004,0.003,0.003),each=N) # Trees
B2 <- rep(c(-1.5,-1,-1),each=N) # Trees
    
## Reproduction curve
B3 <- rep(c(-5,-10,-50),each=N) # Trees
B4 <- rep(c(0.02,0.03,0.08),each=N) # Trees
B5 <- rep(c(-2,-1,-1),each=N) # Trees


## Seedling annual survival
Seedlings$S <- rbinom(N*length(SpFitList),1,
                      prob=plogis(sB0+sB1*Seedlings$height))
## Seedling  annual  growth
Seedlings$G <- rnorm(N*length(SpFitList),
                     mean=rev(B0)+rev(B1)*Seedlings$height,sd=1) 

## random liana score
Trees$liana <- sample(0:4,3*N,prob=c(.6,.1,.1,.1,.1),replace=TRUE)

## Tree annual survival
Trees$S <- rbinom(N*length(SpFitList),1,
                  prob=plogis(B0+B1*Trees$dbh+B2*Trees$liana))
## Tree annual growth
Trees$G <- rnorm(N*length(SpFitList),
                 mean=rev(B0)+rev(B1)*Trees$dbh+rev(B2)*Trees$liana,sd=1)

## Tree annual reproduction
Trees$R <- rbinom(N*length(SpFitList),1,
                  prob=plogis(B3+B4*Trees$dbh+B5*Trees$liana))

## Tree annual crown fraction
Trees$C <- plogis(B3+B4*Trees$dbh+B5*Trees$liana)

## Tree dbh height allometry
Trees$H <- SSlogis(Trees$dbh,rev(2e4*B0)+rev(10*B1)*Trees$dbh,400,410)



## All other vital rates are kept constant for simplicity
## (as not to needlessly increase complexity of an example)

cat("\t ##################################################### \n ")
cat("\t \t \t Done creating example datasets \n")
cat("\t ##################################################### \n ")


################################################################################
################################################################################
## Step 2: Fit vital rate models
################################################################################
################################################################################


## Timing starts.. this may take some time (about 3-4 min on most systems)
## Speed up by setting N lower above

StartT <- Sys.time()
cat("\t ##################################################### \n ")
cat( "\t \t \t Fitting vital rate models \n")
cat("\t ##################################################### \n ")

## First Seedlings
## Standardize size (aids in model convergence)
## Growth
Seedlings$size <- scale(Seedlings$height)
temp <- unlist(attributes(Seedlings$size)[2:3])
sdlgrownorm <- list(size=c(temp[1],temp[2]))
sdlgrowmain.saved <- lmer(G~size+(1+size|sp),data=Seedlings)

## Check the linear assumption of the model
## Red lines should be centred around zero
## and show no trend, or curves
linearityCheck(sdlgrowmain.saved,"size")

##in this case sigma is equal among species
sdlsigmas <- array(rep(sigma(sdlgrowmain.saved),
                       length(SpFitList)))
names(sdlsigmas) <- SpFitList

## Survival
sdlsurvmain.saved <- glmer(S~size+(1+size|sp),data=Seedlings,family=binomial)
sdlsurvnorm <- sdlgrownorm

## Note that whenever models warn about convergence problems.
## This can be caused by scale differences between variables.
## The tolerance is optimal for variables on a scale 
## close to -1 to 1. Some models here have
## other scales, and  may converge on the correct
## solution despite the warning.
## User can test this by checking the relative
## tolerance and seeing if it is accectable.
##  User can also update the model to
## test if more iterations solves this.
## The models above and below seem to have found the correct parameters
## (the ones used to create the data), despite the
## warnings.

## Warnings on convergence, update model to run try again
sdlsurvmain.saved <- update(sdlsurvmain.saved)

## still complaining check relative values
relgrad <- with(sdlsurvmain.saved@optinfo$derivs,solve(Hessian,gradient))

## Check if the gradiant is fine 
max(abs(relgrad))<0.001 ## model can be assumed to have converged 

## Check the linear assumption of the model
linearityCheck(sdlsurvmain.saved,"size")

## Now Trees

## Growth
Trees$size <- scale(Trees$dbh) # Normalization again
growmain.saved <- lmer(G~size+as.factor(liana)+(1+size+as.factor(liana)|sp),
                       data=Trees)

## Check the linear assumption of the model
linearityCheck(growmain.saved,"size")

temp <- unlist(attributes(Trees$size)[2:3])
grownorm <- list(size=c(temp[1],temp[2]))

##in this case sigma is equal among species
treesigmas <- array(rep(sigma(sdlgrowmain.saved),
                        length(SpFitList)))
names(treesigmas) <- SpFitList

## Survival             
survmain.saved <- glmer(S~size+as.factor(liana)+(1+size+as.factor(liana)|sp),
                        data=Trees,family=binomial)

## Check the linear assumption of the model
linearityCheck(survmain.saved,"size")

survnorm <- grownorm

## Reproduction             
repmain.saved <- glmer(R~size+as.factor(liana)+(1+size+as.factor(liana)|sp),
                     data=Trees,family=binomial)
repnorm <- grownorm

## Check the linear assumption of the model
linearityCheck(repmain.saved,"size")

## crown fraction
fecmain.saved <- lmer(C~size+as.factor(liana)+(1+size+as.factor(liana)|sp),
                     data=Trees)
fecnorm <- grownorm

## here we see linearity problems (created on purpose, to illustrate)
linearityCheck(fecmain.saved,"size") 

## Allomatric model
s2tmod.saved <- lmer(dbh~H+(1+dbh|sp), data=Trees)

## CONSTANT VITAL RATES (for simplicity)
## seed production
seedprod.saved <- data.frame(sp=SpFitList, fec=rep(1,length(SpFitList)))
## recruitment and height distributions

recruitrate.saved <- data.frame(sp=SpFitList, rrate=rep(.1,length(SpFitList)))

sdlheightresults.saved <- data.frame(sp=SpFitList, v=rep(1.05,length(SpFitList)),
                               lambda=rep(12.5,length(SpFitList)),
                               shape=rep(1.05,length(SpFitList)),
                               scale=rep(1.05,length(SpFitList)))

## Size ranges
Sizeranges <- data.frame(sp=SpFitList, maxdbh=rep(2000,length(SpFitList)))



EndT <- Sys.time()

cat("\t ##################################################### \n ")
cat( "\t \t All models fit. Work started at \n")
cat( "\t \t",format(StartT),"\n")
cat( "\t \t Work ended at \n")
cat( "\t \t",format(EndT),"\n")
cat( "\t \t Which took \n")
cat( "\t \t",format(EndT-StartT),"\n")
cat("\t ##################################################### \n ")
cat("\n")

################################################################################
## Setup bootstrap

Nboots <- 10
Lianaeigenvalues <- array(dim=c(length(SpFitList),5,Nboots),
                          dimnames=list(SpFitList,
                                        c("l0","l1","l2","l3","l4"),
                                        1:Nboots))


## Start bootstrap operations
merModSamplerFE <- function(mod,n=1) {
mvrnorm(n,fixef(mod),vcov(mod))
}


cat("\t ##################################################### \n ")
cat("\t Starting bootstrap with",Nboots, " iterations for ",length(SpFitList), "Species \n ")
cat("\t ##################################################### \n ")




## time the execution
StartT <- Sys.time()

for(j in 1:Nboots){
      
################################################################################
    ## iterate model parameters 
################################################################################
    ## NON MerMod
    ## subset all fit objects to focal species
    seedprod <- seedprod.saved
    recruitrate <- recruitrate.saved
    sdlheightresults <- sdlheightresults.saved

    ## MerMod
    
    ## allometric equations
    s2tmod <- s2tmod.saved
    s2tmod@beta <- merModSamplerFE(s2tmod)
    IPM_s2t <- data.frame(sp=rownames(ranef(s2tmod)$sp),
                          intercept=fixef(s2tmod)[1]+ranef(s2tmod)$sp[,1],
                          slope=fixef(s2tmod)[2]+ranef(s2tmod)$sp[,2])
    
    ## structure list(modelname, list(modelobject, norms))
    ## or  list(modelname, modelobject) if no norms
    repmain <- repmain.saved
    fecmain <- fecmain.saved
    repmain@beta <- merModSamplerFE(repmain)
    fecmain@beta <- merModSamplerFE(fecmain)
    
    sdlgrowmain <- sdlgrowmain.saved
    sdlgrowmain@beta <- merModSamplerFE(sdlgrowmain)
    sdlsurvmain <- sdlsurvmain.saved
    sdlsurvmain@beta <- merModSamplerFE(sdlsurvmain)
    growmain <- growmain.saved
    growmain@beta <- merModSamplerFE(growmain)
    survmain <- survmain.saved
    survmain@beta <- merModSamplerFE(survmain)

    
    ## Start building all IPMs
######################################################################
    DimSdl <- 100
    DimTree <- 400

    

    for(i in 1: length(SpFitList)) {
        
        foc.sp <- SpFitList[i]
        
        cat("\r \t Bootstrap ", j, "Species # ", i,"    ")

        sizerange<-c(0,Sizeranges$maxdbh[Sizeranges$sp==foc.sp]*1.1)
        

        hghtdbhthreshold <- sdl2splg(dbh=10,model=IPM_s2t
                                    ,sp=SpFitList[i])
        dbhthreshold <-  sdl2splg(hght=hghtdbhthreshold,model=IPM_s2t
                                 ,sp=SpFitList[i])
        
        seedlingclasses <- seq(sizerange[1],hghtdbhthreshold,length.out=DimSdl)
                                        # remove final class
        seedlingclasses <- seedlingclasses[-length(seedlingclasses)]
        seedlingclasswidth<-mean(diff(seedlingclasses))
        
        treeclasses<-seq(dbhthreshold,sizerange[2],length.out=DimTree)
        treeclasswidth<-mean(diff(treeclasses))
        
        classes<-c(seedlingclasses,treeclasses)

        ## subset all fit objects to focal species
        f.sp <- seedprod[seedprod$sp==foc.sp,c("fec","sp")]
        f.sts <- recruitrate[recruitrate$sp==foc.sp,c("rrate","sp")]

        ## structure list(modelname, list(modelobject, norms))
        ## or  list(modelname, modelobject) if no norms
        f.repro <- list("lmepredict",repmain,
                    list("lmepredict",list(fecmain,as.numeric(fecnorm$size))))
        
        f.sdlhght <- data.frame("dweibull",
                                sdlheightresults$shape[sdlheightresults$sp==foc.sp],
                                sdlheightresults$scale[sdlheightresults$sp==foc.sp])
        names(f.sdlhght) <- c("bestmod","p_1","p_2")
        
        f.s2t <- IPM_s2t[IPM_s2t$sp==foc.sp,]

        f.sg <- list(list(sdlgrowmain,sdlgrownorm$size,
                          sdlsigmas[names(sdlsigmas)==foc.sp])
                    ,"lmepredict")

        f.ss <- list(list(sdlsurvmain,sdlsurvnorm$size),"lmepredict")
        
        f.tg <- list(list(growmain,grownorm$size,treesigmas[names(treesigmas)==foc.sp])
                    ,"lmepredict")
        f.ts <- list(list(survmain,survnorm$size),"lmepredict")


        ## No elasticity calculations needed for this step
        Imat <- ItMat(vitalrate="tree growth",delta=1,size=1,
                      seedlingclasses,treeclasses)

        ## set the liana infestation level
        GlobalLianaLevel <- NULL

        ## Build the IPM
        fullIPM <- constructIPM(seedlingclasses,treeclasses,
                                seedlingclasseswidth,treeclasseswidth,
                                sizerange,hghtdbhthreshold,dbhthreshold,
                                f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                                f.ss,f.tg,f.ts,Imat=Imat)

        ## set the liana infestation level
        GlobalLianaLevel <- 1
        ## Build the IPM
        Liana1fullIPM <- constructIPM(seedlingclasses,treeclasses,
                                      seedlingclasseswidth,treeclasseswidth,
                                      sizerange,hghtdbhthreshold,dbhthreshold,
                                      f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                                      f.ss,f.tg,f.ts,Imat=Imat)[c("IPM","GSM")]

        ## Build the IPM
        GlobalLianaLevel <- 2
        Liana2fullIPM <- constructIPM(seedlingclasses,treeclasses,
                                      seedlingclasseswidth,treeclasseswidth,
                                      sizerange,hghtdbhthreshold,dbhthreshold,
                                      f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                                      f.ss,f.tg,f.ts,Imat=Imat)[c("IPM","GSM")]
        ## Build the IPM
        GlobalLianaLevel <- 3
        Liana3fullIPM <- constructIPM(seedlingclasses,treeclasses,
                                      seedlingclasseswidth,treeclasseswidth,
                                      sizerange,hghtdbhthreshold,dbhthreshold,
                                      f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                                      f.ss,f.tg,f.ts,Imat=Imat)[c("IPM","GSM")]
        ## Build the IPM
        GlobalLianaLevel <- 4
        Liana4fullIPM <- constructIPM(seedlingclasses,treeclasses,
                                      seedlingclasseswidth,treeclasseswidth,
                                      sizerange,hghtdbhthreshold,dbhthreshold,
                                      f.sp,f.sts,f.repro,f.sdlhght,f.s2t,f.sg,
                                      f.ss,f.tg,f.ts,Imat=Imat)[c("IPM","GSM")]
        
        Lianaeigenvalues [i,,j] <- c(lambda(fullIPM[["IPM"]]),
                                     lambda(Liana1fullIPM[["IPM"]]),
                                  lambda(Liana2fullIPM[["IPM"]]),
                                  lambda(Liana3fullIPM[["IPM"]]),
                                  lambda(Liana4fullIPM[["IPM"]]))
        
    }



   
}
EndT <- Sys.time()

cat("\t ##################################################### \n ")
cat( "\t \t Calculations complete. Work started at \n")
cat( "\t \t",format(StartT),"\n")
cat( "\t \t Work ended at \n")
cat( "\t \t",format(EndT),"\n")
cat( "\t \t Which took \n")
cat( "\t \t",format(EndT-StartT),"\n")
cat("\t ##################################################### \n ")
cat("\n")

