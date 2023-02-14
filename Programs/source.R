#==============================================================================
# FILENAME: sim_source.R
# PROJECT: 	RMST pseudo-observation evaluation
# PURPOSE: source functions to simulate trial results 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 4.2.0
#==============================================================================
#Notes: 





# =============================================================================

# install.packages("survival")
# install.packages("survRM2")
# install.packages("pseudo")
# install.packages("geepack")
# install.packages("flexsurv")
# install.packages("plyr")
# install.packages("data.table")
# install.packages("ldbounds")

library(geepack)
library(survival)
library(survRM2)
library(pseudo)
library(flexsurv)
library(plyr)
library(data.table)
library(ldbounds)


simdat.cov <- function(n,               #total number of observations
                       nevents,         # total sample size
                       accru,           # accrual rate per month
                       LTFU,            # loss to follow up proportion
                       pos,             # proportion of biomarker positive subjects
                       trtB,            # proportion of subjects randomly assigned to treatment B
                       female,          # probability of being female
                       age.min,         # minimum age of subjects
                       age.max,         # maximum age of subjects
                       dis1p,            # probability of having disease severity 1
                       dis2p,            # probability of having disease severity 2
                       dis3p,            # probability of having disease severity 3
                       dis4p,            # probability of having disease severity 4
                       b0,              # the median survival for a male of minimum age with disease status 1, biomarker negative, treated with A
                       bstar,             # coefficient for being treated with treatment B
                       bmark,            # the coefficient for being biomarker positive
                       b1,              # the coefficient for age in the model to generate survival times
                       b2,              # coefficient for being female
                       b3,              # coefficient for having disease status 2
                       b4,              # coefficient for having disease status 3
                       b5,              # coefficient for having disease status 4
                       b6,               # coefficient for having disease status 5
                       b7,               # coefficient for interaction (biomarker positive and treated with B)
                       PH=TRUE,          # when true, proportional hazards hold and survival times are exponential
                       # when false, treatment B subjects are simulated using a weibull with 
                       # the input shape parameter, treatment A still exponential
                       shapeB=3,           # the shape of the weibull parameter when PH is false
                       clin.trial=TRUE    # when true, will simulate a clinical trial of nevents
                       # when false, will simulate n event times
){
     # generating time between accruing each patient based on the accrual rate
     times <- rexp(n=n, rate=accru)
     # generating study entry times based on time between recruiting each patient
     entry.time <- cumsum(times)-times[1]
     
     # generating marker status for each patient, 1=marker positive, 0=marker negative
     marker.stat <- rbinom(n, size=1, prob=pos)
     
     # getting treatment assignment (1=B, 0=A)
     trtn <- rbinom(n, size=1, prob=trtB)
     
     # getting female status (1=female, 0=male)
     female <- rbinom(n, size=1, prob=female)
     
     # getting age, based on an exponential distribution (older more likely than younger)
     range <- age.max-age.min
     lambda.age <- -log(.01)/range
     age <- round(age.max - rexp(n, rate=lambda.age))
     
     ## getting each of the disease status variables
     dis1 <- rbinom(n, size=1, prob=dis1p)
     dis2 <- fifelse(dis1==0, dis2 <- rbinom(n, size=1, prob=(dis2p/(1-dis1p))), dis2 <- 0)
     dis3 <- fifelse((dis1==0 & dis2==0), dis3 <- rbinom(n, size=1, prob=(dis3p/((1-dis1p)*(1-dis2p)))), dis3 <- 0)
     dis4 <- fifelse((dis1==0 & dis2==0 & dis3==0), dis4 <- rbinom(1, size=1, prob=(dis4p/((1-dis1p)*(1-dis2p)*(1-dis3p)))), dis4 <- 0)
     dis5 <- fifelse((dis1==0 & dis2==0 & dis3==0 & dis4==0), dis5 <- 1, dis5 <- 0)
     
     dis.stat <- dis1 + (2*dis2) + (3*dis3) + (4*dis4) + (5*dis5)
     
     #  truncating the age variable and centering
     age.cent <- age - age.min
     age.cent <- pmax(age.cent,0)
     
     # generating subject ids in order of accrual/screening
     id <- seq(from=1, to=n, by=1)
     
     
     ### getting median survival times for each subject
     med.surv <- b0 + b1*log((age.cent+1)) + b2*female + b3*dis2 + b4*dis3 + b5*dis4 +
          b6*dis5 + bstar*trtn + bmark*marker.stat + b7*marker.stat*trtn
     
     ### getting survival times for each subject based on an exponential distribution and event time
     if (PH==TRUE){
          surv.time <- rexp(n=n, rate=log(2)/med.surv)
          
     }
     else {
          exp <- rexp(n=n, rate=log(2)/med.surv)
          weib <- rweibull(n=n, shape=shapeB, scale=med.surv/(log(2)^(1/shapeB)))
          surv.time <- trtn*weib + (1-trtn)*exp
          # for (i in 1:n){
          #   if(trtn[i]==1){surv.time[i] <- rweibull(n=1, shape=shapeB, scale=med.surv/(log(2)^(1/shapeB)))}
          #   else {surv.time[i] <- rexp(n=1, rate=log(2)/med.surv)}
          # }
     }
     event <- rbinom(n=n,size=1,prob=(1-LTFU))
     event.time <- entry.time + surv.time
     
     lnt <- log(surv.time)
     
     
     data <- matrix(c(id, entry.time, marker.stat, trtn, female, age,age.cent, dis.stat, dis1, dis2, dis3, 
                      dis4, dis5,event, surv.time,lnt, event.time), nrow=n)
     
     names <- c("id","entry.time","marker.stat","trtn","female","age","age.cent", "dis.stat","dis1",
                "dis2","dis3","dis4","dis5","event","surv.time","lnt", "event.time")
     colnames(data) <- names
     
     if (clin.trial==TRUE){
          
          temp <- data[data[,"event"]==1,]
          temp <- temp[order(temp[,"event.time"],decreasing=FALSE),]
          end.time <- temp[nevents,"event.time"]
          
          data <- data[data[,"entry.time"]<=end.time,]
          
          for (i in 1:nrow(data)){
               if (data[i,"event.time"]>end.time){data[i,"event"] =0
               data[i,"surv.time"]=end.time - data[i,"entry.time"]}
          }
          data <- matrix(c(data,rep(end.time,nrow(data)),rep(nrow(data),nrow(data))),nrow=nrow(data))
          colnames(data) <- c(names,"dur","enr")
     }
     
     data <- data[order(data[,"surv.time"],decreasing=FALSE),]
     
     return(data)
}



# function to estimate difference in RMST, along with standard error estimates, using the RP(3,1) flexible parametric model
# the treatment variable in the dataset must be named trtn and included last in the model statement
surv.est.rp31 <- function(data, # survival data used to model the RP model
                          tstar, # upper limit of integration for RMST calculation
                          fit,   # an RP(3,1) model fit using flexsurvspline; treatment variable must be named trtn, it must be
                          # included in the model statement last, and be 
                          # the variable allowing for change in effect over time (in the anc list)
                          alpha,  # two-sided type 1 error
                          var.typ="asymp"   ## name of the variance type to use. 2 options:
                          #"asymp" uses the correct, asymptotic variance
                          # "delta" uses the delta method variance
){
     if(length(fit$cov)>=4){
          
          # calculates the knots and lambdas from the model fit
          kmin <- as.numeric(fit$knots[1])
          k1 <- as.numeric(fit$knots[2])
          k2 <- as.numeric(fit$knots[3])
          kmax <- as.numeric(fit$knots[4])
          lam1 <- (kmax-k1)/(kmax-kmin)
          lam2 <- (kmax-k2)/(kmax-kmin)
          
          # number of non-treatment covariates included in the model
          ncovs <- (length(fit$coefficients)-6)
          
          # creating data frame to store individual results for the sandwich estimator
          res.ind <- data.frame(matrix(nrow=length(data[,1]), ncol=(20+(3*ncovs))))
          names1 <- NULL
          names0 <- NULL
          bnames <- NULL
          if (ncovs>0){
               for (t in 1:ncovs){
                    names1 <- c(names1,paste("amatc",t,"1",sep=""))
                    names0 <- c(names0,paste("amatc",t,"0",sep=""))
                    bnames <- c(bnames,paste("bmatc",t,sep=""))
               }
               
          }
          colnames(res.ind) <- c("rmst1","amatb01","amatb11","amatb21","amatb31","amatbtrt","amatbtrtint",names1,
                                 "rmst0","amatb00","amatb10","amatb20","amatb30","amatbtrt0","amatbtrtint0",names0,
                                 "bmatb0","bmatb1","bmatb2","bmatb3","bmatbtrt","bmatbtrtint",bnames)
          
          # setting the results for treatment variables to 0 when we force treatment to be 0
          res.ind$amatbtrt0 <- 0
          res.ind$amatbtrtint0 <- 0
          # getting rmst estimates for each subject with and without treatment
          
          # functions that will be integrated over t based on the model to get the estimated rmst along with elements
          # needed for the sandwich estimator
          # the first functions are for when treatment is set to 1
          
          # these are the spline basis functions
          z1 <- function(x){
               (ifelse(log(x) > k1, ((log(x)-k1)^3),0) - (lam1*(ifelse(log(x) > kmin, ((log(x)-kmin)^3),0))) - ((1-lam1)*(ifelse(log(x) > kmax, ((log(x)-kmax)^3),0))))
          }
          z1prime <- function(x){
               (ifelse(log(x) > k1, ((3/x)*((log(x)-k1)^2)),0) - (lam1*(ifelse(log(x) > kmin, ((3/x)*((log(x)-kmin)^2)),0))) - ((1-lam1)*(ifelse(log(x) > kmax, ((3/x)*((log(x)-kmax)^2)),0))))
          }
          z2 <- function(x){
               (ifelse(log(x) > k2, ((log(x)-k2)^3),0) - (lam2*(ifelse(log(x) > kmin, ((log(x)-kmin)^3),0))) - ((1-lam2)*(ifelse(log(x) > kmax, ((log(x)-kmax)^3),0))))
          }
          z2prime <- function(x){
               (ifelse(log(x) > k2, ((3/x)*((log(x)-k2)^2)),0) - (lam2*(ifelse(log(x) > kmin, ((3/x)*((log(x)-kmin)^2)),0))) - ((1-lam2)*(ifelse(log(x) > kmax, ((3/x)*((log(x)-kmax)^2)),0))))
          }
          
          # the cumulative log hazard
          cum.haz1 <- function(x,covhaz=tmp.cov.haz){
               fit$coefficients[1] + fit$coefficients[2]*log(x) +  fit$coefficients[3]*z1(x) +  fit$coefficients[4]*z2(x) +
                    fit$coefficients["trtn"] + fit$coefficients["gamma1(trtn)"]*log(x) + covhaz
          }
          # derivative of the cumulative log hazard wrt t
          lnHprime1 <- function(x, trt=trtn){
               (fit$coefficients[2]/x) +  fit$coefficients[3]*z1prime(x) +  fit$coefficients[4]*z2prime(x) + ((fit$coefficients["gamma1(trtn)"]*trt)/x)
          }
          # survival function in terms of t
          surv.ind1 <- function(x){
               exp(-exp(cum.haz1(x)))
          }
          # negative derivative of the survival function wrt Beta 0
          a.mat.b01 <- function(x){
               exp(cum.haz1(x))*surv.ind1(x)
          }
          # negative derivative of the survival function wrt Beta 1
          a.mat.b11 <- function(x){
               exp(cum.haz1(x))*surv.ind1(x)*log(x)
          }
          # negative derivative of the survival function wrt Beta 2
          a.mat.b21 <- function(x){
               exp(cum.haz1(x))*surv.ind1(x)*z1(x)
          }
          # negative derivative of the survival function wrt Beta 3
          a.mat.b31 <- function(x){
               exp(cum.haz1(x))*surv.ind1(x)*z2(x)
          }
          # negative derivative of the survival function wrt Beta 4
          a.mat.btrt <- function(x){
               exp(cum.haz1(x))*surv.ind1(x)
          }
          # negative derivative of the survival function wrt Beta 5
          a.mat.btrt.int <- function(x){
               exp(cum.haz1(x))*surv.ind1(x)*log(x)
          }
          # the function for the covariates are generated below
          
          # these functions are for when treatment is set to 0
          
          # the cumulative log hazard
          cum.haz0 <- function(x, covhaz=tmp.cov.haz){
               fit$coefficients[1] + fit$coefficients[2]*log(x) +  fit$coefficients[3]*z1(x) +  fit$coefficients[4]*z2(x) + covhaz
          }
          # derivative of the cumulative log hazard wrt t
          lnHprime0 <- function(x, trt=trtn){
               (fit$coefficients[2]/x) +  fit$coefficients[3]*z1prime(x) +  fit$coefficients[4]*z2prime(x)
          }
          # survival function in terms of t
          surv.ind0 <- function(x){
               exp(-exp(cum.haz0(x)))
          }
          # negative derivative of the survival function wrt Beta 0
          a.mat.b00 <- function(x){
               exp(cum.haz0(x))*surv.ind0(x)
          }
          # negative derivative of the survival function wrt Beta 1
          a.mat.b10 <- function(x){
               exp(cum.haz0(x))*surv.ind0(x)*log(x)
          }
          # negative derivative of the survival function wrt Beta 2
          a.mat.b20 <- function(x){
               exp(cum.haz0(x))*surv.ind0(x)*z1(x)
          }
          # negative derivative of the survival function wrt Beta 3
          a.mat.b30 <- function(x){
               exp(cum.haz0(x))*surv.ind0(x)*z2(x)
          }
          
          # the following first calculates the covariate contribution to the log cumulative hazard for each subject
          # and integrates each of the functions above to obtain the rmst and elements of the bread matrix for each subject
          # calculation of the covariate part of the bread matrix is at the bottom of this chunk
          for (i in 1:length(data[,1])){
               # getting the covariate part of the log hazard function for each subject
               tmp.cov.haz <- 0
               if (ncovs>0){
                    for (j in 5:(4+ncovs)){
                         tmp.cov.haz <- tmp.cov.haz + fit$res[j,1]*data[i,rownames(fit$res)[j]]
                    }
               }
               # output when treatment is set to 1
               res.ind[i,"rmst1"] <- integrate(surv.ind1, lower=0, upper=tstar)$value
               res.ind[i,"amatb01"] <- integrate(a.mat.b01, lower=0, upper=tstar)$value
               res.ind[i,"amatb11"] <- integrate(a.mat.b11, lower=0, upper=tstar)$value
               res.ind[i,"amatb21"] <- integrate(a.mat.b21, lower=0, upper=tstar)$value
               res.ind[i,"amatb31"] <- integrate(a.mat.b31, lower=0, upper=tstar)$value
               res.ind[i,"amatbtrt"] <- integrate(a.mat.btrt, lower=0, upper=tstar)$value
               res.ind[i,"amatbtrtint"] <- integrate(a.mat.btrt.int, lower=0, upper=tstar)$value
               # output when treatment is set to 0
               res.ind[i,"rmst0"] <- integrate(surv.ind0, lower=0, upper=tstar)$value
               res.ind[i,"amatb00"] <- integrate(a.mat.b00, lower=0, upper=tstar)$value
               res.ind[i,"amatb10"] <- integrate(a.mat.b10, lower=0, upper=tstar)$value
               res.ind[i,"amatb20"] <- integrate(a.mat.b20, lower=0, upper=tstar)$value
               res.ind[i,"amatb30"] <- integrate(a.mat.b30, lower=0, upper=tstar)$value
               
               # this is a loop to compose a function for each covariate included in the model, for each subject
               # each function is then integrated to provide the covariate part of the bread matrix
               if (ncovs>0){
                    ind <- i
                    for (k in 1:(ncovs)){
                         row <- rownames(fit$res)[(k+4)]
                         # output when treatment set to 1
                         assign(paste("a.mat.c",k,"1",sep=""),function(x, cnt=ind, rname=row){exp(cum.haz1(x))*surv.ind1(x)*data[cnt,rname]})
                         res.ind[i,(7+k)] <- integrate(paste("a.mat.c",k,"1",sep=""), lower=0, upper=tstar)$value 
                         # output when treatment set to 0
                         assign(paste("a.mat.c",k,"0",sep=""),function(x, cnt=ind, rname=row){exp(cum.haz0(x))*surv.ind0(x)*data[cnt,rname]})
                         res.ind[i,(14+ncovs+k)] <- integrate(paste("a.mat.c",k,"0",sep=""), lower=0, upper=tstar)$value
                    }
               }
          }
          
          ### point estimate for difference in rmst
          est. <- mean(res.ind[,"rmst1"])-mean(res.ind[,"rmst0"])
          
          ### calculating the asymptotic variance
          if (var.typ=="asymp"){
               
               ### constructing the A matrix, i.e., the bread in the sandwich estimator
               # starting with the covariance matrix for the betas om the upper left corner
               amat <- data.frame((1/length(data[,1]))*solve(matrix(fit$cov, nrow=nrow(fit$cov),ncol=ncol(fit$cov))))
               # getting the mean of the a matrix integrals calculated above in order of the betas from lowest subscript to highest
               # (the treatment betas are last)
               trtmeans <- NULL
               notrtmeans <- NULL
               for (i in 1:4){
                    trtmeans <- c(trtmeans,mean(res.ind[,(i+1)]))
                    notrtmeans <- c(notrtmeans,mean(res.ind[,(8+ncovs+i)]))           
               }
               if (ncovs>0){
                    for (i in 1:ncovs){
                         trtmeans <- c(trtmeans,mean(res.ind[,(i+7)]))
                         notrtmeans <- c(notrtmeans,mean(res.ind[,(14+ncovs+i)]))
                    }
               }
               for (i in 1:2){
                    trtmeans <- c(trtmeans,mean(res.ind[,(i+5)]))
                    notrtmeans <- c(notrtmeans,mean(res.ind[,(12+ncovs+i)]))
               }
               # binding the treatment=1 and treatment=0 columns to the beta covariance matrix along
               # with a column of zeros to reflect the a matrix in the paper
               amat <- cbind(amat,trtmeans,notrtmeans, rep(0,(6+ncovs)))
               amat <- rbind(amat,c(rep(0,(6+ncovs)),1,0,-1),c(rep(0,(6+ncovs)),0,1,1),c(rep(0,(6+ncovs)),-1,1,1))
               amat <- t(amat)
               
               
               # code to get score contributions from the flexsurv output
               optpars <- fit$optpars
               Y <- fit$data$Y
               X <- flexsurv:::compress.model.matrices(fit$data$mml)
               xweights <- rep(1, nrow(Y))
               bhazard <- rep(0, nrow(Y))
               rtrunc <- rep(Inf, nrow(X))
               dlist <- fit$dlist
               inits <- fit$coefficients
               
               nk <- length(fit$knots)
               dfn <- unroll.function(dsurvspline, gamma=0:(nk-1))
               pfn <- unroll.function(psurvspline, gamma=0:(nk-1))
               rfn <- unroll.function(rsurvspline, gamma=0:(nk-1))
               hfn <- unroll.function(hsurvspline, gamma=0:(nk-1))
               Hfn <- unroll.function(Hsurvspline, gamma=0:(nk-1))
               qfn <- unroll.function(qsurvspline, gamma=0:(nk-1))
               meanfn <- unroll.function(mean_survspline, gamma=0:(nk-1))
               rmstfn <- unroll.function(rmst_survspline, gamma=0:(nk-1))
               Ddfn <- unroll.function(flexsurv:::DLdsurvspline, gamma=0:(nk-1))
               DSfn <- unroll.function(flexsurv:::DLSsurvspline, gamma=0:(nk-1))
               
               dfns=list(d=dfn,p=pfn,r=rfn,h=hfn,H=Hfn,rmst=rmstfn,mean=meanfn, q=qfn,
                         DLd=Ddfn,DLS=DSfn,deriv=TRUE)
               
               aux <- list(knots=fit$knots, scale=fit$scale, timescale=fit$timescale)
               mx <- fit$mx
               fixedpars <- NULL
               
               pars <- inits
               npars <- length(pars)
               
               nbpars <- length(dlist$pars)
               pars <- as.list(pars)
               ncovs <- length(pars) - length(dlist$pars)
               if (ncovs > 0)
                    beta <- unlist(pars[(nbpars+1):npars])
               for (i in dlist$pars) {
                    if (length(mx[[i]]) > 0)
                         pars[[i]] <- pars[[i]] + X[,mx[[i]],drop=FALSE] %*% beta[mx[[i]]]
                    else
                         pars[[i]] <- rep(pars[[i]], length(Y[,"stop"]))
               }
               dead <- Y[,"status"]==1
               ddcall <- list(t=Y[dead,"stop"])
               dsccall <- list(t=Y[!dead,"stop"])
               dstcall <- list(t=Y[,"start"])
               for (i in 1:nbpars)
                    ddcall[[names(pars)[i]]] <-
                    dsccall[[names(pars)[i]]] <-
                    dstcall[[names(pars)[i]]] <-
                    dlist$inv.transforms[[i]](pars[[i]])
               for (i in seq_along(aux)){
                    ddcall[[names(aux)[i]]] <- dsccall[[names(aux)[i]]] <-
                         dstcall[[names(aux)[i]]] <- aux[[i]]
               }
               for (i in dlist$pars) {
                    ddcall[[i]] <- ddcall[[i]][dead]
                    dsccall[[i]] <- dsccall[[i]][!dead]
               }
               dd <- flexsurv:::dderiv(dfns$DLd, ddcall, X[dead,,drop=FALSE], mx, dlist)
               dscens <- flexsurv:::dderiv(dfns$DLS, dsccall, X[!dead,,drop=FALSE], mx, dlist)
               if (sum(dead) > 0) dd <- dd * xweights[dead]
               if (sum(!dead) > 0) dscens <- dscens * xweights[!dead]
               dstrunc <- flexsurv:::dderiv(dfns$DLS, dstcall, X, mx, dlist) * xweights
               res.individ <- - ( rbind(dd, dscens))
               
               res.individ <- matrix(NA, ncol = ncol(dd), nrow = nrow(X))
               res.individ[which(dead),] <- -dd
               res.individ[which(!dead), ] <- -dscens
               
               bmat.pre <- data.frame(cbind(res.individ,
                                            rmst1=(res.ind[,"rmst1"]-mean(res.ind[,"rmst1"])),rmst0=(res.ind[,"rmst0"]-mean(res.ind[,"rmst0"])),
                                            rmstdiff=(res.ind[,"rmst1"]-res.ind[,"rmst0"]-(mean(res.ind[,"rmst1"])-mean(res.ind[,"rmst0"])))))
               
               bmat <- data.frame(matrix(nrow=ncol(bmat.pre),ncol=ncol(bmat.pre)))
               for(i in 1:ncol(bmat.pre)){
                    for (j in 1:ncol(bmat.pre)){
                         bmat[i,j] <- mean(bmat.pre[,i]*bmat.pre[,j])
                    }
               }
               colnames(bmat) <- colnames(bmat.pre)
               
               cov.mat <- ((solve(data.matrix(amat))) %*% data.matrix(bmat) %*% (t(solve(data.matrix(amat)))))/length(data[,1])
               
               var <- cov.mat[ncol(bmat),ncol(bmat)]
               se <- sqrt(var)
               upper <- est.+(se*qnorm(1-(alpha/2)))
               lower <- est.-(se*qnorm(1-(alpha/2)))
               p <- 2*(1-pnorm(abs(est./se)))
          }
          else if (var.typ=="delta"){
               # log survival
               log.surv1 <- function(x){
                    -exp(cum.haz1(x))
               }
               func1 <- function(x){
                    log.surv1(x)*surv.ind1(x)
               }
               log.surv0 <- function(x){
                    -exp(cum.haz0(x))
               }
               func0 <- function(x){
                    log.surv0(x)*surv.ind0(x)
               }
               funcdiff <- function(x){
                    func1(x)-func0(x)
               }
               h1 <- function(x){
                    funcdiff(x)*log(x)
               }
               h2 <- function(x){
                    funcdiff(x)*z1(x)
               }
               h3 <- function(x){
                    funcdiff(x)*z2(x)
               }
               h5 <- function(x){
                    func1(x)*log(x)
               }
               
               ## getting the gradient values for each subject
               result <- matrix(nrow=nrow(data), ncol=(6+ncovs))
               #### integrates the gradient functions for each subject
               for (i in 1:nrow(data)){
                    # getting the covariate part of the log hazard function for each subject
                    tmp.cov.haz <- 0
                    if (ncovs>0){
                         for (j in 5:(4+ncovs)){
                              tmp.cov.haz <- tmp.cov.haz + fit$res[j,1]*data[i,rownames(fit$res)[j]]
                         }
                    }
                    result[i,1] <- integrate(funcdiff,lower=0,upper=tstar)$value
                    result[i,2] <- integrate(h1,lower=0,upper=tstar)$value
                    result[i,3] <- integrate(h2,lower=0,upper=tstar)$value
                    result[i,4] <- integrate(h3,lower=0,upper=tstar)$value
                    result[i,(ncovs+5)] <- integrate(func1,lower=0,upper=tstar)$value
                    result[i,(ncovs+6)] <- integrate(h5,lower=0,upper=tstar)$value
                    if (ncovs>0){
                         ind <- i
                         for (k in 1:ncovs){
                              row <- rownames(fit$res)[(k+4)]
                              assign(paste("cov",k,sep=""),function(x, cnt=ind, rname=row){funcdiff(x)*data[cnt,rname]})
                              result[i,(k+4)] <- integrate(f=paste("cov",k,sep=""), lower=0, upper=tstar)$value
                         }
                    }
               }
               hb <- matrix(colMeans(result),nrow=1)
               var <- (hb %*% matrix(fit$cov,nrow=(6+ncovs),ncol=6+ncovs) %*% t(hb))
               se <- sqrt(var)
               upper <- est.+(se*qnorm(1-(alpha/2)))
               lower <- est.-(se*qnorm(1-(alpha/2)))
               p <- 2*(1-pnorm(abs(est./se)))
          }
          
     }     
     if(length(fit$cov)<4) {
          est.=NA
          se=NA
          lower=NA
          upper=NA
          p=NA
     }
     
     return(c("estimate"=est.,"SE"=se, "lower"=lower,"upper"=upper,"pval"=p))
     
}



# start of code to reproduce Chen
rsum <- function(y, x, tmax) {  # sum of rectangles
     keep <- which(x < tmax)
     width <- diff(c(x[keep], tmax))
     sum(width * y[keep])
}


rowcumSums <- function(x) {
     
     addmat <- matrix(0, nrow = ncol(x), ncol = ncol(x))
     addmat[lower.tri(addmat, diag = TRUE)] <-1
     
     t(addmat %*% t(x))
     
}

### chen cox method of rmst estimation
chen.cox <- function(data,           # data frame with survival variables and covariates
                     covs=NULL,      # vector of variable names to include in the model(s); must be in data
                     tstar,          # upper limit of integration for rmst calculation
                     p=NULL          # probability of being assigned treatment; when null, uses empirical estimates
){
     data <- data.frame(data)
     form <- as.formula(paste0("Surv(surv.time, event) ~", paste(covs, collapse = "+")))
     coxfit0 <- coxph(form, data = subset(data, trtn == 0))
     coxfit1 <- coxph(form, data = subset(data, trtn == 1))
     
     etimes <- sort(data$surv.time[data$event == 1])
     etimes <- etimes[etimes <= tstar]
     Ai <- data$trtn[match(etimes, data$surv.time)]
     
     basedat <- do.call(rbind, lapply(etimes, \(tt) {
          ret <- data[, covs, drop = FALSE]
          ret$id <- 1:nrow(ret)
          ret$surv.time <- tt
          ret$event <- 1
          ret
     }))
     
     breslow0 <- matrix(predict(coxfit0, newdata = basedat, type = "expected"), 
                        nrow = nrow(data), ncol = length(etimes))
     breslow1 <- matrix(predict(coxfit1, newdata = basedat, type = "expected"), 
                        nrow = nrow(data), ncol = length(etimes))
     
     risk0 <- predict(coxfit0, newdata = data, type = "risk", reference = "zero")
     risk1 <- predict(coxfit1, newdata = data, type = "risk", reference = "zero")
     
     S0hat <- colMeans(exp(-breslow0))
     S1hat <- colMeans(exp(-breslow1))
     
     rmst0 <- rsum(c(1,S0hat), c(0,etimes), tstar)
     rmst1 <- rsum(c(1,S1hat), c(0,etimes), tstar)
     
     diffrmst <- rmst1 - rmst0
     
     
     ## individual counting process
     
     Nit <- do.call(cbind, lapply(etimes, \(tt){
          
          ifelse(data$surv.time <= tt & data$event == 1, 1, 0)
          
     }))
     
     Yit <- outer(data$surv.time, etimes, \(x,y) 1.0 * (x >= y))
     
     
     h0hattmp <- matrix(diff(c(0, etimes, tstar)), nrow = nrow(data), ncol = length(etimes) + 1, byrow = TRUE) * 
          cbind(1,exp(-breslow0)) * matrix(risk0, nrow = nrow(data), 
                                           ncol = length(etimes)+1)
     h0hat <- colMeans(rowcumSums(h0hattmp[, ncol(h0hattmp):1])[,ncol(h0hattmp):1])[-1]
     
     h1hattmp <- matrix(diff(c(0, etimes, tstar)), nrow = nrow(data), ncol = length(etimes) + 1, byrow = TRUE) * 
          cbind(1,exp(-breslow1)) * matrix(risk1, nrow = nrow(data), 
                                           ncol = length(etimes)+1)
     h1hat <- colMeans(rowcumSums(h1hattmp[, ncol(h1hattmp):1])[, ncol(h1hattmp):1])[-1]
     
     
     SS0.t.0 <- colMeans(matrix(risk0 * (data$trtn == 0), nrow = nrow(data), 
                                ncol = length(etimes), byrow = FALSE) * Yit)
     
     
     SS1.t.0 <- matrix(NA, nrow = length(etimes), ncol = length(covs))
     SS2.t.0 <- array(NA, dim = c(length(etimes), length(covs), length(covs)))
     
     SS0.t.1 <- colMeans(matrix(risk1 * (data$trtn == 1), nrow = nrow(data), 
                                ncol = length(etimes), byrow = FALSE) * Yit)
     
     
     SS1.t.1 <- matrix(NA, nrow = length(etimes), ncol = length(covs))
     SS2.t.1 <- array(NA, dim = c(length(etimes), length(covs), length(covs)))
     
     Sigma0 <- Sigma1 <- matrix(0, nrow = length(covs), ncol = length(covs))
     
     if(length(covs) > 1) { 
          dmat <- as.matrix(data[, covs]) 
          
          
     } else { 
          dmat <- as.matrix(data[, covs, drop = FALSE])
          
          
     }
     
     ZZmat <-  do.call(rbind, 
                       apply(dmat, MAR = 1, FUN = \(x) c(x %*% t(x)), 
                             simplify = FALSE))
     
     for(i in 1:length(etimes)) {
          
          SS1.t.0[i, ] <- colMeans(matrix((data$trtn == 0) * risk0 * Yit[, i], nrow = nrow(data), 
                                          ncol = length(covs)) * dmat)
          
          
          
          SS2.t.0[i,,] <- matrix(colMeans(matrix((data$trtn == 0) * risk0 * Yit[, i], nrow = nrow(data), 
                                                 ncol = length(covs)^2) * 
                                               ZZmat), 
                                 nrow = length(covs), 
                                 ncol = length(covs))
          
          SS1.t.1[i, ] <- colMeans(matrix((data$trtn == 1) * risk1 * Yit[, i], nrow = nrow(data), 
                                          ncol = length(covs)) * dmat)
          
          
          SS2.t.1[i,,] <- matrix(colMeans(matrix((data$trtn == 1) * risk1 * Yit[, i], nrow = nrow(data), 
                                                 ncol = length(covs)^2) * 
                                               ZZmat), 
                                 nrow = length(covs), 
                                 ncol = length(covs))
          
          
          Sigma0 <- Sigma0 + (1 - Ai[i]) * (SS2.t.0[i,,] / SS0.t.0[i] - 
                                                 (SS1.t.0[i,] / SS0.t.0[i]) %*% t(SS1.t.0[i,] / SS0.t.0[i]))
          Sigma1 <- Sigma1 + (Ai[i]) * (SS2.t.1[i,,] / SS0.t.1[i] - 
                                             (SS1.t.1[i,] / SS0.t.1[i]) %*% t(SS1.t.1[i,] / SS0.t.1[i]))
          
     }
     
     Sigma0 <- Sigma0 / sum(1 - data$trtn)
     Sigma1 <- Sigma1 / sum(data$trtn)
     
     
     Zbar <- SS1.t.0 / matrix(SS0.t.0, nrow = length(SS0.t.0), ncol = ncol(SS1.t.0))
     Zbar1 <- SS1.t.1 / matrix(SS0.t.1, nrow = length(SS0.t.1), ncol = ncol(SS1.t.1))
     
     denom.haz0 <- colSums(matrix(risk0 * (data$trtn == 0), nrow = length(risk0), ncol = ncol(Yit)) * Yit)
     denom.haz1 <- colSums(matrix(risk1 * (data$trtn == 1), nrow = length(risk1), ncol = ncol(Yit)) * Yit)
     g0i <- g1i <- matrix(NA, nrow = nrow(data), ncol = length(covs))
     
     
     for(i in 1:nrow(data)) {
          
          
          inmat <- (1 - Ai) * (Zbar - dmat[replicate(nrow(Zbar), i), ]) / 
               matrix(denom.haz0, nrow = length(denom.haz0), ncol = ncol(Zbar))
          
          inmat1 <- Ai * (Zbar1 - dmat[replicate(nrow(Zbar1), i), ]) / 
               matrix(denom.haz1, nrow = length(denom.haz1), ncol = ncol(Zbar1))
          
          Lcurv0 <-  rbind(0, matrix(exp(-breslow0[i, ]) * risk0[i], nrow = nrow(inmat), 
                                     ncol = ncol(inmat)) *apply(inmat, MAR = 2, cumsum))
          
          Lcurv1 <-  rbind(0, matrix(exp(-breslow1[i, ]) * risk1[i], nrow = nrow(inmat1), 
                                     ncol = ncol(inmat1)) *apply(inmat1, MAR = 2, cumsum))
          
          g0i[i,] <- apply(Lcurv0, MAR = 2, rsum, c(0, etimes), tstar)
          g1i[i,] <- apply(Lcurv1, MAR = 2, rsum, c(0, etimes), tstar)
          
          
     }
     
     g0 <- (solve(Sigma0) %*% colSums(g0i)) / sum(1 - data$trtn)
     g1 <- (solve(Sigma1) %*% colSums(g1i)) / sum(data$trtn)
     
     var.diff <- (sum(1 - data$trtn) / nrow(data)) * t(g0) %*% Sigma0 %*% g0 + 
          sum((1 - Ai) * (h0hat^2 / SS0.t.0) / denom.haz0) + 
          (sum(data$trtn) / nrow(data)) * t(g1) %*% Sigma1 %*% g1 +
          sum(Ai * (h1hat^2 / SS0.t.1) / denom.haz1) + 
          (nrow(data) - 1) / (nrow(data)) * var(apply(cbind(1,exp(-breslow1)), MAR = 1, rsum, c(0,etimes),tstar) - 
                                                     apply(cbind(1,exp(-breslow0)), MAR = 1, rsum, c(0,etimes),tstar))
     
     
     list(est = diffrmst, var = var.diff / nrow(data), 
          g0 = g0, g1 = g1, Sigma0 = Sigma0, Sigma1 = Sigma1, 
          hint0 = sum((1 - Ai) * (h0hat^2 / SS0.t.0) / denom.haz0), 
          hint1 = sum(Ai * (h1hat^2 / SS0.t.1) / denom.haz1))
}

chen.cox2 <- function(data,           # matrix with survival variables and covariates
                      covs=NULL,      # vector of variable names to include in the model(s); must be in data
                      tstar,          # upper limit of integration for rmst calculation
                      p=NULL          # probability of being assigned treatment; when null, uses empirical estimates
){
     #getting number of covariates
     ncovs <- length(covs)
     # ordering data by increasing survival time
     data <- data[order(data[,"surv.time"],decreasing=FALSE),]
     # sample sizes, total and per treatment arm
     N <- nrow(data)
     n1 <- sum(data[,"trtn"])
     n0 <- sum(1-data[,"trtn"])
     # subsetting data by treatment arm
     trt <- trt <- data[data[,"trtn"]==1,]
     notrt <- data[data[,"trtn"]==0,]
     # probability of being assigned/receiving treatment 1
     prob1 <- ifelse(is.null(p),(n1/N),p)
     # getting the model per user input covariates
     formula.chen <- as.formula(paste0("Surv(surv.time, event) ~", paste(covs, collapse = "+")))
     ## fitting the cox models separately for each treatment arm
     coxfit1 <- coxph(formula.chen, data=subset(data.frame(data), trtn==1))
     coxfit0 <- coxph(formula.chen, data=subset(data.frame(data), trtn==0))
     
     # getting event times and keeping only those before tstar and tstar itself
     times1 <- trt[trt[,"trtn"]==1 & trt[,"event"]==1 & trt[,"surv.time"] <= tstar,"surv.time"]
     times0 <- notrt[notrt[,"trtn"]==0 & notrt[,"event"]==1 & notrt[,"surv.time"] <= tstar,"surv.time"]
     
     # getting the intervals between times
     int1 <- diff(c(0,times1, tstar))
     int0 <- diff(c(0,times0, tstar))
     
     ## getting survival estimates for each arm for all subjects. rows are the estimates, column for each subject
     survest1 <- exp(-survfit(coxfit1, newdata=data.frame(data))$cumhaz)
     survest0 <- exp(-survfit(coxfit0, newdata=data.frame(data))$cumhaz)
     
     ## getting the indices of only the event times on or before tstar
     index1 <- which(trt[,"event"]==1 & (trt[,"surv.time"]<=tstar))
     index0 <- which(notrt[,"event"]==1 & (notrt[,"surv.time"]<=tstar))
     
     ## keeping only the survival estimates that correspond to events on or before tstar and including 1 for time 0
     survest1 <- rbind(rep(1,N),survest1[index1,])
     survest0 <- rbind(rep(1,N),survest0[index0,])
     
     ### calculating rmst the other way; calculate for each subject and then take mean
     rmst1 <- apply(survest1, MARGIN=2, function(x){sum(x*int1)})
     rmst0 <- apply(survest0, MARGIN=2, function(x){sum(x*int0)})
     rmstdiff <- rmst1-rmst0
     rmst.est <- mean(rmstdiff)
     
     #### getting the variance of the estimated difference in RMST; sigma function above this one
     
     # calculating S0, S1 and S2 for each event time and each treatment group
     ## just the covariates from data
     covs <- as.matrix(data[,covs])
     ## calculating the hazard multiplier due to covariates
     covhaz1 <- exp(covs %*% coxfit1$coefficients)
     covhaz0 <- exp(covs %*% coxfit0$coefficients)
     temp <- colnames(data)
     data <- cbind(data, "covhaz0"=covhaz0, "covhaz1"=covhaz1)
     colnames(data) <- c(temp,"covhaz0","covhaz1")
     # getting 1's for indicator that subject is on treatment arm, 0 else
     Ai <- data[,"trtn"]
     ## # row for each subject; column for each event time. a 1 if subject was at risk at the event time and on that arm
     Yit0 <- outer((data[,"surv.time"]*(1-Ai)), times0, \(x,y) 1.0 * (x >= y))
     Yit1 <- outer((data[,"surv.time"]*Ai), times1, \(x,y) 1.0 * (x >= y))
     
     ### covhaz for only those on the treatment arm at each event time
     covhaz0arm <- matrix(covhaz0,nrow=nrow(Yit0),ncol=ncol(Yit0))*Yit0
     covhaz1arm <- matrix(covhaz1,nrow=nrow(Yit1),ncol=ncol(Yit1))*Yit1
     
     ### getting S0, S1 and S2 for each arm
     S00 <- apply(covhaz0arm, MARGIN=2, function(x){sum(x)/N})
     
     S01 <- apply(covhaz1arm, MARGIN=2, function(x){sum(x)/N})
     
     S10 <- apply(as.matrix(covs), MARGIN=2, function(x)
     {temp <- matrix(x, nrow=nrow(covhaz0arm), ncol=ncol(covhaz0arm))*covhaz0arm
     return(as.numeric(colMeans(temp)))})
     
     S11 <- apply(as.matrix(covs), MARGIN=2, function(x)
     {temp <- matrix(x, nrow=nrow(covhaz1arm), ncol=ncol(covhaz1arm))*covhaz1arm
     return(as.numeric(colMeans(temp)))})
     
     Zbar0 <- S10/S00
     
     Zbar1 <- S11/S01
     
     S20 <- alply(covhaz0arm, 2, function(x){
          temp <- alply(covs, 1, function(y){
               (as.matrix(y) %*% t(as.matrix(y)))
          })
          temp2 <- do.call(cbind,temp)
          temp2 <- array(temp2, dim=c(dim(temp[[1]]), length(temp)))
          temp3 <- temp2*rep(x, each=ncovs^2)
          return(apply(temp3,c(1,2),sum)/N)})
     
     S21 <- alply(covhaz1arm, 2, function(x){
          temp <- alply(covs, 1, function(y){
               (as.matrix(y) %*% t(as.matrix(y)))
          })
          temp2 <- do.call(cbind,temp)
          temp2 <- array(temp2, dim=c(dim(temp[[1]]), length(temp)))
          temp3 <- temp2*rep(x, each=ncovs^2)
          return(apply(temp3,c(1,2),sum)/N)})
     
     
     ## getting sigma0 and sigma1
     temp <- do.call(cbind, S20)
     temp <- array(temp, dim=c(dim(S20[[1]]), length(S20)))
     temp2 <- temp/rep(S00, each=ncovs^2)
     
     temp3 <- alply(Zbar0, 1, function(x){temp <- as.matrix(x) %*% t(as.matrix(x))})
     temp4 <- do.call(cbind,temp3)
     temp4 <- array(temp4, dim=c(dim(temp3[[1]]), length(temp3)))
     temp5 <- temp2 - temp4
     sigma0 <- apply(temp5, c(1,2), sum)/n0
     
     temp <- do.call(cbind, S21)
     temp <- array(temp, dim=c(dim(S21[[1]]), length(S21)))
     temp2 <- temp/rep(S01, each=ncovs^2)
     
     temp3 <- alply(Zbar1, 1, function(x){temp <- as.matrix(x) %*% t(as.matrix(x))})
     temp4 <- do.call(cbind,temp3)
     temp4 <- array(temp4, dim=c(dim(temp3[[1]]), length(temp3)))
     temp5 <- temp2 - temp4
     sigma1 <- apply(temp5, c(1,2), sum)/n1
     
     ### getting the h hats at each event time and the value of the integrals with h^2/S0*dcumhaz
     # this matrix is the survival estimates at each event time in the arm * covhaz0 * time intervals
     temp0 <- matrix(int0, nrow=nrow(data), ncol=(length(times0)+1), byrow=TRUE)*
          t(survest0)*matrix(covhaz0, nrow=nrow(data), ncol=(length(times0)+1))
     # reversing order
     temp02 <- temp0[,ncol(temp0):1]
     # cumulative row sums
     temp02 <- rowcumSums(temp02)
     # drop last column
     temp02 <- temp02[,-ncol(temp02)]
     # reverse order back to original
     temp02 <- temp02[,ncol(temp02):1]
     # column means
     h0 <- colMeans(temp02)
     ### this is the value of the integral of h^2/S0 * dcumhaz
     ## note that the denominator in dcumhaz equals S0*N
     tmp <- h0/(N*S00)
     tmp2 <- h0/S00
     hint0 <- sum(tmp*tmp2)
     
     temp1 <- matrix(int1, nrow=nrow(data), ncol=(length(times1)+1), byrow=TRUE)*
          t(survest1)*matrix(covhaz1, nrow=nrow(data), ncol=(length(times1)+1))
     # reversing order
     temp12 <- temp1[,ncol(temp1):1]
     # cumulative row sums
     temp12 <- rowcumSums(temp12)
     # drop last column
     temp12 <- temp12[,-ncol(temp12)]
     # reverse order back to original
     temp12 <- temp12[,ncol(temp12):1]
     # column means
     h1 <- colMeans(temp12)
     ### this is the value of the integral of h^2/S0 * dcumhaz
     ## note that the denominator in dcumhaz equals S0*N
     tmp1 <- h1/(N*S01)
     tmp21 <- h1/S01
     hint1 <- sum(tmp1*tmp21)
     
     ### getting the values of the g0 hat and g1 hat
     
     denom0 <- S00*N
     denom1 <- S01*N
     ## g0
     inmat <- alply(covs,1,function(x){
          rbind(0,apply(((Zbar0 - matrix(as.matrix(x),nrow=nrow(Zbar0),ncol=ncol(Zbar0),byrow=TRUE))/denom0),2,cumsum))
     })
     inmat2 <- do.call(cbind,inmat)
     inmat2 <- array(inmat2, dim=c(dim(inmat[[1]]), length(inmat)))
     
     survcov <- alply(temp0,1,function(x){matrix(rep(x,ncovs),nrow=length(int0), ncol=ncovs)})
     survcov2 <- do.call(cbind,survcov)
     survcov2 <- array(survcov2, dim=c(dim(survcov[[1]]), length(survcov)))
     
     g0itmp <- inmat2*survcov2
     g0i <- apply(g0itmp,c(1,2),sum)
     
     ###g1
     inmat1 <- alply(covs,1,function(x){
          rbind(0,apply(((Zbar1 - matrix(as.matrix(x),nrow=nrow(Zbar1),ncol=ncol(Zbar1),byrow=TRUE))/denom1),2,cumsum))
     })
     inmat21 <- do.call(cbind,inmat1)
     inmat21 <- array(inmat21, dim=c(dim(inmat1[[1]]), length(inmat1)))
     
     survcov1 <- alply(temp1,1,function(x){matrix(rep(x,ncovs),nrow=length(int1), ncol=ncovs)})
     survcov21 <- do.call(cbind,survcov1)
     survcov21 <- array(survcov21, dim=c(dim(survcov1[[1]]), length(survcov1)))
     
     g1itmp <- inmat21*survcov21
     g1i <- apply(g1itmp,c(1,2),sum)
     
     g0 <- (solve(sigma0) %*% colSums(g0i))/n0
     g1 <- (solve(sigma1) %*% colSums(g1i))/n1
     
     mse <- mean((rmstdiff - rmst.est)^2)
     
     var.est <- ((n0/N)*(t(g0) %*% sigma0 %*% g0) + hint0 +
                      (n1/N)*(t(g1) %*% sigma1 %*% g1) + hint1 + mse)/N
     
     return(list("mean.diff"=rmst.est, "var.est"=var.est,"mse"=mse,"g1"=g1,
                 "g0"=g0,"sigma1"=sigma1,"sigma0"=sigma0,"hint0"=hint0,"hint1"=hint1))
}




#### wrapper function to include all RMST estimation methods into 1 funciton

RMSTdiff <- function(data,              #### data frame which includes time, event and covariate variables
                     time.var,          #### time to event variable name in data
                     event.var,         #### binary event indicator variable name in data (1=event, 0=censored)
                     trtn.var,           #### binary treatment group variable name in data (must be binary; only two groups can be compared)
                     covs=NULL,         #### vector of covariate names to adjust for in rmst estimation
                     tmax=NULL,         #### time point up to which rmst is calculated; can be NULL
                     method="pseudo",   #### method of rmst evaluation, default is the pseudo-observation technique
                     ### method options:
                     ### "pseudo" - technique using the pseudo-observations to compute rmst
                     ### "KM" - KM method without covariates from survRM2 package
                     ### "Tian" - the Tian method with covariates from the survRM2 package
                     ### "RP31" - Method using the Royston Parmar (3,1) flexible parametric model
                     ### "Chen" - estimation via the Cox model proposed by Chen and Tsiatis Biometrics 2004, model 2
                     var.type1="ajk",      ### variance type to report; only relevant for the "pseudo" method
                     ### options for pseudo are: 
                     ### "ajk" for approximate jackknife 
                     ### "sandwich" for the sandwhich estimator
                     var.type2="asymp",  ### variance type to report; only relevant for the "RP31" method
                     ### options for RP31 are:
                     ###   "asymp" for the asymptotic variance
                     ### "naive" for the delta method variance
                     alpha=.05,            ### two-sided type 1 error
                     rnd.prb1=NULL            ### probability of being assigned to trtn=1; if NULL, uses empirical estimate
                     
){
     data <- data.frame(data)
     #### renaming the variables in data for use in the below functions
     data$trtn <- data[,trtn.var]
     data$surv.time <- data[,time.var]
     data$event <- data[,event.var]
     
     ## calculates the maximum follow-up time to compute rmst in case where tmax is NULL
     if (is.null(tmax)){
          trt1 <- data[data$trtn==1,]
          trt1 <- trt1[order(trt1$surv.time),]
          trt0 <- data[data$trtn==0,]
          trt0 <- trt0[order(trt0$surv.time),]
          max1 <- max(trt1$surv.time)
          max0 <- max (trt0$surv.time)
          tmax <- min(max0,max1)
          # if ((max0 < max1) & trt0$event[nrow(trt0)]==1){
          #      tmax <- max1
          # }
          # else if ((max0 > max1) & trt1$event[nrow(trt1)]==0){
          #      tmax <- max1
          # }
          # else {tmax <- max0}
     }
     
     if (method=="pseudo"){
          formula <- as.formula(paste0("pseudom ~", paste(covs, collapse = "+"),"+ trtn"))
          ## calculating the pseudo-observations
          pseudom <- pseudomean(data$surv.time, data$event, tmax=tmax)
          ## appending pseudo-observations to the data
          data <- data.frame(cbind(data,pseudom=pseudom))
          ## fitting the gee model 
          if (var.type1=="ajk"){
               fit <- geese(formula, data=data, jack=TRUE, family=gaussian, corstr="independence", scale.fix=FALSE)
               sum.fit <- round(cbind(mean=fit$beta, SE=sqrt(diag(fit$vbeta.ajs)), Z=fit$beta/sqrt(diag(fit$vbeta.ajs)),
                                      PVal = 2-2*pnorm(abs(fit$beta/sqrt(diag(fit$vbeta.ajs))))),4)
               
               rmst <- sum.fit["trtn","mean"]
               rmst.se <- sum.fit["trtn","SE"]
               rmst.lb <- sum.fit["trtn","mean"] - (sum.fit["trtn","SE"]*qnorm(1-(alpha/2)))
               rmst.ub <- sum.fit["trtn","mean"] + (sum.fit["trtn","SE"]*qnorm(1-(alpha/2)))
               rmst.pval <- sum.fit["trtn","PVal"]
          }
          else if (var.type1=="sandwich"){
               fit <- geese(formula, data=data, family=gaussian, corstr="independence", scale.fix=FALSE)
               sum.fit <- summary(fit)$mean
               
               rmst <- sum.fit["trtn","estimate"]
               rmst.se <- sum.fit["trtn","san.se"]
               rmst.lb <- sum.fit["trtn","estimate"] - (sum.fit["trtn","san.se"]*qnorm(1-(alpha/2)))
               rmst.ub <- sum.fit["trtn","estimate"] + (sum.fit["trtn","san.se"]*qnorm(1-(alpha/2)))
               rmst.pval <- sum.fit["trtn","p"]
          }
     }
     
     else if (method=="KM"){
          if (is.null(covs)){
               fit <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, alpha=alpha)
               
               rmst <- fit$unadjusted.result["RMST (arm=1)-(arm=0)","Est."]
               rmst.se <- (fit$unadjusted.result["RMST (arm=1)-(arm=0)",3] - fit$unadjusted.result["RMST (arm=1)-(arm=0)","Est."])/qnorm(1-(alpha/2))
               rmst.lb <- fit$unadjusted.result["RMST (arm=1)-(arm=0)",2]
               rmst.ub <- fit$unadjusted.result["RMST (arm=1)-(arm=0)",3]
               rmst.pval <- fit$unadjusted.result["RMST (arm=1)-(arm=0)","p"]
          }
          else {
               rmst <- NA
               rmst.se <- NA
               rmst.lb <- NA
               rmst.ub <- NA
               rmst.pval <- NA
          }
     }
     else if (method=="Tian"){
          if (!is.null(covs)){     
               cov.mat <- data[,covs]
               fit <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, covariates=cov.mat, alpha=alpha)
               
               rmst <- fit$RMST.difference.adjusted["arm","coef"]
               rmst.se <- fit$RMST.difference.adjusted["arm","se(coef)"]
               rmst.lb <- fit$RMST.difference.adjusted["arm",5]
               rmst.ub <- fit$RMST.difference.adjusted["arm",6]
               rmst.pval <- fit$RMST.difference.adjusted["arm","p"]
               
          }
          else{
               rmst <- NA
               rmst.se <- NA
               rmst.lb <- NA
               rmst.ub <- NA
               rmst.pval <- NA  
          }
     }
     
     else if (method=="RP31"){
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          fit <- flexsurvspline(formula, anc=list(gamma1=~trtn), data=data, k=2, scale="hazard")
          result <- surv.est.rp31(data=data, tstar=tmax, fit=fit, alpha=alpha, var.typ=var.type2)
          
          rmst <- as.numeric(result["estimate"])
          rmst.se <- as.numeric(result["SE"])
          rmst.lb <- as.numeric(result["lower"])
          rmst.ub <- as.numeric(result["upper"])
          rmst.pval <- as.numeric(result["pval"])
     }
     
     else if (method=="Chen"){
          if (!is.null(covs)){
               result <- chen.cox(data=data, covs=covs, tstar=tmax, p=rnd.prb1)
               rmst <- result$est
               rmst.se <- sqrt(result$var)
               rmst.lb <- rmst - (rmst.se*qnorm(1-(alpha/2)))
               rmst.ub <- rmst + (rmst.se*qnorm(1-(alpha/2)))
               rmst.pval <- 2*(1-pnorm(abs(rmst/rmst.se)))
          }
          else{
               rmst <- NA
               rmst.se <- NA
               rmst.lb <- NA
               rmst.ub <- NA
               rmst.pval <- NA
          }
     }
     
     
     return(c("RMSTdiff"=rmst,"SE"=rmst.se, "lower.bound"=rmst.lb, "upper.bound"=rmst.ub,"pval"=rmst.pval))   
}

     
#### function to simulate 1 trial

sim.trial <- function(data,           # survival data frame or matrix
                      n,              # number of total events at final analysis
                      int.time,       # timing of interim analyses; last time always 1
                      alpha=0.05,     # two-sided total alpha
                      low.err=0.10,    # error for futility boundary calculation
                      bound.type=c(1,1), # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively  
                      RM.times=NULL,      # value or vector of time points to compare the restricted mean survival probability
                                        # if NULL, uses max follow-up time at each analysis
                      covs=NULL,
                      dataT=NULL       # dataset for calculating the true RMST at different times
){
     # getting the pvalue cutoffs for the analysis times
     x <- ldBounds(t=int.time, # the percentage of information times, last one always =1
                   iuse=bound.type,             # indicates O'brien fleming boundaries for both upper and lower boundaries
                   alpha=c(low.err, (alpha/2)))     # the alpas for the lower and upper boundaries, respectively
     ## upper bounds at each interim time
     up.bd <- x$upper.bounds
     # lower bounds
     low.bd <- x$lower.bounds
     
     # number of events needed at each timing
     n_events <- ceiling(int.time*n)
     
     # preserving original data set
     data2 <- data.frame(data)
     
     ### study times for each interim analysis
     times <- rep(NA,length(int.time))
     enr <- rep(NA,length(int.time))
     ### max follow-up time for each interim analysis
     maxt <- rep(NA,length(int.time))
     
     lr.stop <- 0
     lr.eff <- 0
     lr.fut <- 0
     
     hr.stop <- 0
     hr.eff <- 0
     hr.fut <- 0
     hr <- NA
     
     psu.sand.stop <- 0
     psu.sand.eff <- 0
     psu.sand.fut <- 0
     
     psu.ajk.stop <- 0
     psu.ajk.eff <- 0
     psu.ajk.fut <- 0
     
     km.stop <- 0
     km.eff <- 0
     km.fut <- 0
     
     tian.stop <- 0
     tian.eff <- 0
     tian.fut <- 0
     
     rp.delt.stop <- 0
     rp.delt.eff <- 0
     rp.delt.fut <- 0
     
     rp.m.stop <- 0
     rp.m.eff <- 0
     rp.m.fut <- 0
     
     if (!is.null(covs)){
          cox.stop <- 0
          cox.eff <- 0
          cox.fut <- 0
     }
     
     
     ### a loop to test for each interim analysis
     for (i in 1:length(int.time)){
          temp <- data2[data2[,"event"]==1,]
          temp <- temp[order(temp[,"event.time"],decreasing=FALSE),]
          end.time <- as.numeric(temp[n_events[i],"event.time"])
          times[i] <- end.time
          
          data <- data2[data2[,"entry.time"]<=end.time,]
          
          for (j in 1:nrow(data)){
               if (data[j,"event.time"]>end.time){data[j,"event"] =0
               data[j,"surv.time"]=(end.time - data[j,"entry.time"])}
          }
          enr[i] <- nrow(data)
          
          ### getting the times to estimate RMST at each interim analysis
          if (is.null(RM.times)){
               trt1 <- data[data$trtn==1,]
               trt0 <- data[data$trtn==0,]
               max1 <- max(trt1$surv.time)
               max0 <- max (trt0$surv.time)
               tmax <- min(max0,max1)
               maxt[i] <- tmax
               t <- tmax
          }
          else {
               t <- RM.times[i]
          }
          
          ### testing based on logrank test
          if (lr.stop==0){
               lr <- survdiff(Surv(surv.time, event) ~ trtn, data=data)
               exp_num <- lr$exp[2] - lr$obs[2]
               if (sqrt(lr$chisq) > up.bd[i] & exp_num > 0){
                    lr.stop <- i
                    lr.eff <- 1
               }
               else if (-sqrt(lr$chisq) < low.bd[i] & exp_num < 0){
                    lr.stop <- i
                    lr.fut <- 1
               }
          }
          ### testing using cox proportional hazards
          if (hr.stop==0){
               formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
               temp <- NULL
               temp <- coxph(formula,data=data)
               p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
               hr <- summary(temp)$coefficients["trtn","coef"]
               if((qnorm(1-p) > up.bd[i]) & summary(temp)$coefficients["trtn","coef"] < 0){
                    hr.stop <- i
                    hr.eff <- 1
               }
               else if ((qnorm(p) < low.bd[i]) & summary(temp)$coefficients["trtn","coef"] > 0){
                    hr.stop <- i
                    hr.fut <- 1
               }
          }  
          
          ### testing using cox pseudo method with sandwich estimate
          if (psu.sand.stop==0){
               pseudo.sand <- try(RMSTdiff(data=data,
                                           time.var="surv.time",
                                           event.var="event",
                                           trtn.var="trtn",
                                           covs=covs,
                                           tmax=t,
                                           method="pseudo",
                                           var.type1="sandwich",
                                           var.type2="asymp",
                                           alpha=.05))
               
               if(inherits(pseudo.sand, "try-error")){
                    pseudo.sand <- rep(NA,5)
               }
               if (is.na(pseudo.sand["pval"]) | is.null(pseudo.sand["pval"])){}
               else if ((qnorm(1-pseudo.sand["pval"]) > up.bd[i]) & pseudo.sand["RMSTdiff"] > 0){
                    psu.sand.stop <- i
                    psu.sand.eff <- 1
               }
               else if ((qnorm(pseudo.sand["pval"]) < low.bd[i]) & pseudo.sand["RMSTdiff"] < 0){
                    psu.sand.stop <- i
                    psu.sand.fut <- 1
               }
          }
          
          ### testing using cox pseudo method with approximate jackknife estimate
          if (psu.ajk.stop==0){
               pseudo.ajk <- try(RMSTdiff(data=data,
                                          time.var="surv.time",
                                          event.var="event",
                                          trtn.var="trtn",
                                          covs=covs,
                                          tmax=t,
                                          method="pseudo",
                                          var.type1="ajk",
                                          var.type2="asymp",
                                          alpha=.05))
               if(inherits(pseudo.ajk, "try-error")){
                    pseudo.ajk <- rep(NA,5)
               }
               if (is.na(pseudo.ajk["pval"]) | is.null(pseudo.ajk["pval"])){}
               else if ((qnorm(1-pseudo.ajk["pval"]) > up.bd[i]) & pseudo.ajk["RMSTdiff"] > 0){
                    psu.ajk.stop <- i
                    psu.ajk.eff <- 1
               }
               else if ((qnorm(pseudo.ajk["pval"]) < low.bd[i]) & pseudo.ajk["RMSTdiff"] < 0){
                    psu.ajk.stop <- i
                    psu.ajk.fut <- 1
               }
          }
          
          ### running the test using the km method
          if (km.stop==0 & is.null(covs)){
               km <- try(RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=t,
                                  method="KM",
                                  var.type1="ajk",
                                  var.type2="asymp",
                                  alpha=.05))
               if(inherits(km, "try-error")){
                    km <- rep(NA,5)
               }
               if (is.na(km["pval"]) | is.null(km["pval"])){}
               else if ((qnorm(1-km["pval"]) > up.bd[i]) & km["RMSTdiff"] > 0){
                    km.stop <- i
                    km.eff <- 1
               }
               else if ((qnorm(km["pval"]) < low.bd[i]) & km["RMSTdiff"] < 0){
                    km.stop <- i
                    km.fut <- 1
               }
          }
          
          ### running the test using the Tian method
          if (tian.stop==0 & !is.null(covs)){
               tian <- try(RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=t,
                                  method="Tian",
                                  var.type1="ajk",
                                  var.type2="asymp",
                                  alpha=.05))
               if(inherits(tian, "try-error")){
                    tian <- rep(NA,5)
               }
               if (is.na(tian["pval"]) | is.null(tian["pval"])){}
               else if ((qnorm(1-tian["pval"]) > up.bd[i]) & tian["RMSTdiff"] > 0){
                    tian.stop <- i
                    tian.eff <- 1
               }
               else if ((qnorm(tian["pval"]) < low.bd[i]) & tian["RMSTdiff"] < 0){
                    tian.stop <- i
                    tian.fut <- 1
               }
          }
          
          ### running the test using the RP31 method with delta estimate
          if (rp.delt.stop==0){
               rp31.delt <- try(RMSTdiff(data=data,
                                          time.var="surv.time",
                                          event.var="event",
                                          trtn.var="trtn",
                                          covs=covs,
                                          tmax=t,
                                          method="RP31",
                                          var.type1="ajk",
                                          var.type2="delta",
                                          alpha=.05))
               if(inherits(rp31.delt, "try-error")){
                    rp31.delt <- rep(NA,5)
               }
               if (is.na(rp31.delt["pval"]) | is.null(rp31.delt["pval"])){}
               else if ((qnorm(1-rp31.delt["pval"]) > up.bd[i]) & rp31.delt["RMSTdiff"] > 0){
                    rp.delt.stop <- i
                    rp.delt.eff <- 1
               }
               else if ((qnorm(rp31.delt["pval"]) < low.bd[i]) & rp31.delt["RMSTdiff"] < 0){
                    rp.delt.stop <- i
                    rp.delt.fut <- 1
               }
          }
          
          ### running the test using the RP31 method with M estimate
          if (rp.m.stop==0){
               rp31.asymp <- try(RMSTdiff(data=data,
                                          time.var="surv.time",
                                          event.var="event",
                                          trtn.var="trtn",
                                          covs=covs,
                                          tmax=t,
                                          method="RP31",
                                          var.type1="ajk",
                                          var.type2="asymp",
                                          alpha=.05))
               if(inherits(rp31.asymp, "try-error")){
                    rp31.asymp <- rep(NA,5)
               }
               if (is.na(rp31.asymp["pval"]) | is.null(rp31.asymp["pval"])){}
               else if ((qnorm(1-rp31.asymp["pval"]) > up.bd[i]) & rp31.asymp["RMSTdiff"] > 0){
                    rp.m.stop <- i
                    rp.m.eff <- 1
               }
               else if ((qnorm(rp31.asymp["pval"]) < low.bd[i]) & rp31.asymp["RMSTdiff"] < 0){
                    rp.m.stop <- i
                    rp.m.fut <- 1
               }
          }
          
          if (!is.null(covs)){
               if (cox.stop==0){
                    cox <- try(RMSTdiff(data=data,
                                        time.var="surv.time",
                                        event.var="event",
                                        trtn.var="trtn",
                                        covs=covs,
                                        tmax=t,
                                        method="Chen",
                                        var.type1="ajk",
                                        var.type2="asymp",
                                        alpha=.05))
                    if(inherits(cox, "try-error")){
                         cox <- rep(NA,5)
                    }
                    if (is.na(cox["pval"]) | is.null(cox["pval"])){}
                    else if ((qnorm(1-cox["pval"]) > up.bd[i]) & cox["RMSTdiff"] > 0){
                         cox.stop <- i
                         cox.eff <- 1
                    }
                    else if ((qnorm(cox["pval"]) < low.bd[i]) & cox["RMSTdiff"] < 0){
                         cox.stop <- i
                         cox.fut <- 1
                    }
               }
          }
     }
     trmst <- rep(NA,length(int.time))
     if (!is.null(0)){
          for (i in 1:length(trmst)){
               if (is.null(dataT)){
                    lam0 <- log(2)/b0
                    exp.rmst0 <- function(x){
                         exp(-lam0*x)
                    }
                    lam1 <- log(2)/(b0+bstar)
                    exp.rmst1 <- function(x){
                         exp(-lam1*x)
                    }
                    rmst0 <- integrate(exp.rmst0, lower=0, upper=maxt[i])$value
                    rmst1 <- integrate(exp.rmst1, lower=0, upper=maxt[i])$value
                    trmst[i] <- rmst1-rmst0
               }
               else {
                    fit <- rmst2(time=dataT$surv.time, status=dataT$event, arm=dataT$trtn, tau=maxt[i], alpha=.05)
                    trmst[i] <- fit$unadjusted.result["RMST (arm=1)-(arm=0)","Est."]
               }
          }
     }
     
     
     if (is.null(covs)){
          result <- c("times"=times, "enr"=enr, "maxt"=maxt, "trmst"=trmst,"ubounds"=up.bd, "lbounds"=low.bd,
                      "lr.stop"=lr.stop, "lr.eff"=lr.eff, "lr.fut"=lr.fut,
                      "hr.stop"=hr.stop, "hr.eff"=hr.eff, "hr.fut"=hr.fut, "hr"=hr,
                      "psu.sand.res"=pseudo.sand, "psu.sand.stop"=psu.sand.stop, "psu.sand.eff"=psu.sand.eff, "psu.sand.fut"=psu.sand.fut,
                      "psu.ajk.res"=pseudo.ajk, "psu.ajk.stop"=psu.ajk.stop, "psu.ajk.eff"=psu.ajk.eff, "psu.ajk.fut"=psu.ajk.fut,
                      "km.res"=km, "km.stop"=km.stop, "km.eff"=km.eff, "km.fut"=km.fut,
                      "rp.delt.res"=rp31.delt, "rp.delt.stop"=rp.delt.stop, "rp.delt.eff"=rp.delt.eff, "rp.delt.fut"=rp.delt.fut,
                      "rp.m.res"=rp31.asymp, "rp.m.stop"=rp.m.stop, "rp.m.eff"=rp.m.eff, "rp.m.fut"=rp.m.fut)
     }
     else {
          result <- c("times"=times, "enr"=enr, "maxt"=maxt, "trmst"=trmst,"ubounds"=up.bd, "lbounds"=low.bd,
                      "lr.stop"=lr.stop, "lr.eff"=lr.eff, "lr.fut"=lr.fut,
                      "hr.stop"=hr.stop, "hr.eff"=hr.eff, "hr.fut"=hr.fut, "hr"=hr,
                      "psu.sand.res"=pseudo.sand, "psu.sand.stop"=psu.sand.stop, "psu.sand.eff"=psu.sand.eff, "psu.sand.fut"=psu.sand.fut,
                      "psu.ajk.res"=pseudo.ajk, "psu.ajk.stop"=psu.ajk.stop, "psu.ajk.eff"=psu.ajk.eff, "psu.ajk.fut"=psu.ajk.fut,
                      "tian.res"=tian, "tian.stop"=tian.stop, "tian.eff"=tian.eff, "tian.fut"=tian.fut,
                      "rp.delt.res"=rp31.delt, "rp.delt.stop"=rp.delt.stop, "rp.delt.eff"=rp.delt.eff, "rp.delt.fut"=rp.delt.fut,
                      "rp.m.res"=rp31.asymp, "rp.m.stop"=rp.m.stop, "rp.m.eff"=rp.m.eff, "rp.m.fut"=rp.m.fut,
                      "cox.res"=cox, "cox.stop"=cox.stop, "cox.eff"=cox.eff, "cox.fut"=cox.fut)
     }
     return(result)
     
}
     
     
     
     
     
     
     
     
