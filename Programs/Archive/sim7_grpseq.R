#### Simulation 2 comparing RMST estimation methods


##### RMSTdiff comparison simulation 2
##### this simulation explores power at different sample sizes and same treatment effect and data generated without covariates

source("Programs/source.R")
#### checking type 1 error

reps <- 1000
ss <- c(150,200,250,300,350,400)
LTFU <- .02
bstar <- 2
accru <- 10
rm.times <- NULL
bd.type <- c(1,1)
int.time <- c(0.3,0.7,1)
alpha <- .05
low.err <- 0.1

Sys.time()



for (k in 1:length(ss)){

     result <- matrix(nrow=reps, ncol=52)
     set.seed(12)
     for(i in 1:reps){
          
          data <- simdat.cov(n=(5*ss[k]),               #total number of observations
                             nevents=ss[k],
                             accru=accru,           # accrual rate per month
                             LTFU=LTFU,            # loss to follow up proportion
                             pos=0.5,             # proportion of biomarker positive subjects
                             trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                             female=0.5,          # probability of being female
                             age.min=55,         # minimum age of subjects
                             age.max=85,         # maximum age of subjects
                             dis1=.1,            # probability of having disease severity 1
                             dis2=.2,            # probability of having disease severity 2
                             dis3=.3,            # probability of having disease severity 3
                             dis4=.2,            # probability of having disease severity 4
                             b0=6,              # the median survival for a male of minimum age with disease status 1
                             bstar=bstar,            # increase in median survival for being treated with B
                             bmark=0,             # coefficient for being biomarker positive
                             b1=-1,              # the coefficient for age in the model to generate survival times
                             b2=3,                 # coefficient for being female
                             b3=0,              # coefficient for having disease status 2
                             b4=0,              # coefficient for having disease status 3
                             b5=0,              # coefficient for having disease status 4
                             b6=0,               # coefficient for having disease status 5
                             b7=0,               # coefficent for interaction between trt B and biomarker positive
                             PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                             shapeB=0.5,            # shape of weibull for trt B if PH=FALSE
                             clin.trial=FALSE
          )
          
          covs <- NULL
          
          
          result[i,] <- sim.trial(data,
                                  n=ss[k],
                                  int.time=int.time,
                                  alpha=alpha,
                                  low.err=low.err,
                                  bound.type=bd.type,
                                  RM.times=rm.times,
                                  covs=covs)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"fem.age.nocoef.LTFU",LTFU,"times",rm.times,".rds",sep="")
     
     saveRDS(result, file=paste("Results","sim7grpseq", filetmp, sep="/"))
     
}


for (k in 1:length(ss)){
     
     result <- matrix(nrow=reps, ncol=60)
     set.seed(12)
     for(i in 1:reps){
          
          data <- simdat.cov(n=(5*ss[k]),               #total number of observations
                             nevents=ss[k],
                             accru=accru,           # accrual rate per month
                             LTFU=LTFU,            # loss to follow up proportion
                             pos=0.5,             # proportion of biomarker positive subjects
                             trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                             female=0.5,          # probability of being female
                             age.min=55,         # minimum age of subjects
                             age.max=85,         # maximum age of subjects
                             dis1=.1,            # probability of having disease severity 1
                             dis2=.2,            # probability of having disease severity 2
                             dis3=.3,            # probability of having disease severity 3
                             dis4=.2,            # probability of having disease severity 4
                             b0=6,              # the median survival for a male of minimum age with disease status 1
                             bstar=bstar,            # increase in median survival for being treated with B
                             bmark=0,             # coefficient for being biomarker positive
                             b1=-1,              # the coefficient for age in the model to generate survival times
                             b2=3,                 # coefficient for being female
                             b3=0,              # coefficient for having disease status 2
                             b4=0,              # coefficient for having disease status 3
                             b5=0,              # coefficient for having disease status 4
                             b6=0,               # coefficient for having disease status 5
                             b7=0,               # coefficent for interaction between trt B and biomarker positive
                             PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                             shapeB=0.5,            # shape of weibull for trt B if PH=FALSE
                             clin.trial=FALSE
          )
          
          covs <- c("age","female")
          
          
          result[i,] <- sim.trial(data,
                                  n=ss[k],
                                  int.time=int.time,
                                  alpha=alpha,
                                  low.err=low.err,
                                  bound.type=bd.type,
                                  RM.times=rm.times,
                                  covs=covs)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"fem.age.fem.age.LTFU",LTFU,"times",rm.times,".rds",sep="")
     
     saveRDS(result, file=paste("Results","sim7grpseq", filetmp, sep="/"))
     
}


