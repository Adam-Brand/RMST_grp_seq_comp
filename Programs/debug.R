#==============================================================================
# FILENAME: debug.R
# PROJECT: RMST pseudo-observation evaluation
# PURPOSE: msci. program to aid debugging the other programs 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 4.2.0
#==============================================================================
#Notes: 





# =============================================================================


check <- readRDS("Results/sim1grpseq/PH.ss550trteff3nocov.nocoef.LTFU0.02times.rds")


#debug(simdat.cov)
Sys.time()
set.seed(1212)
data <- simdat.cov(n=1000,               #total number of observations
                   nevents=400,
                   accru=10,           # accrual rate per month
                   LTFU=0,            # loss to follow up proportion
                   pos=0.5,             # proportion of biomarker positive subjects
                   trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                   female=0.5,          # probability of being female
                   age.min=55,         # minimum age of subjects
                   age.max=85,         # maximum age of subjects
                   dis1=.2,            # probability of having disease severity 1
                   dis2=.2,            # probability of having disease severity 2
                   dis3=.3,            # probability of having disease severity 3
                   dis4=.2,            # probability of having disease severity 4
                   b0=6,              # the median survival for a male of minimum age with disease status 1
                   bstar=2,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=-1,              # the coefficient for age in the model to generate survival times
                   b2=3,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3,            # shape of weibull for trt B if PH=FALSE
                   clin.trial=FALSE
)

data <- data.frame(data)
Y <- Surv(data$surv.time, data$event)

kmfit <- survfit(Y~data$trtn)

plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,30),
     xlab="Time in Months",ylab="Survival Probability",main="Scenario 2 Type 1 Error")




debug(sim.trial)
test <- sim.trial(data,
                  n=600,
                  int.time=c(0.3,0.7,1),
                  alpha=.05,
                  low.err=0.1,
                  bound.type=c(1,1),
                  RM.times=NULL,
                  covs=NULL)

test2 <- sim.trial(data,
                  n=600,
                  int.time=c(0.3,0.7,1),
                  alpha=.05,
                  low.err=0.1,
                  bound.type=c(1,1),
                  RM.times=NULL,
                  covs="age")

data <- data.frame(data)
Y <- Surv(data$surv.time, data$event)

kmfit <- survfit(Y~data$trtn)

plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,30),
     xlab="Time in Months",ylab="Survival Probability",main="Scenario 4 Type 1 Error")
legend(22,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)


t <- seq(12,13,by=.1)
result <- rep(NA,length(t))
for (i in 1:length(t)){
     temp <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, tau=t[i])
     result[i] <- abs(temp$unadjusted.result[1,1])
}

t[which.min(result)]


test <- coxph(Surv(surv.time,event) ~trtn+female, data=data)

data <- data.frame(data)
temp <- rmst2(time=data$surv.time, status=data$event, tau=13, arm=data$trtn)

test <- data$Marker
test <- as.integer(test==data$Marker[1])
test <- NULL

data2 <- data[data[,"trtn"]==1,]

median(data[data[,"trtn"]==1,"surv.time"])
median(data[data[,"trtn"]==0,"surv.time"])

lam0 <- log(2)/9
exp.rmst0 <- function(x){
     exp(-lam0*x)
}
lam1 <- log(2)/12
exp.rmst1 <- function(x){
     exp(-lam1*x)
}
rmst0 <- integrate(exp.rmst0, lower=0, upper=Inf)$value
rmst1 <- integrate(exp.rmst1, lower=0, upper=Inf)$value
rmst.diff <- rmst1-rmst0



x<-seq(0,30,by=0.1)
y <- exp(-((log(2)/12)^3)*x^3)

weibull <- function(x){
     integrate(exp(-((log(2)/12)^3)*x^3), lower=0, upper=up)$value
}
z <- function(x){exp(-(log(2)/12)*x)}


curve(exp(-(log(2)/12)*x), from=0,to=30)
curve(exp(-((x/(12/((log(2))^(1/3))))^3)), add=TRUE)


curve(exp(-(log(2)/12)*x), from=0,to=30, col="green",lty=2, xlab="Time in Months",
      ylab="Survival Probability", main="Scenario 1")
curve(exp(-(log(2)/9)*x), from=0,to=30, col="blue",lty=1, add=TRUE)


x <- function(x){exp(-(log(2)/12)*x)}
y <- function(x){exp(-(log(2)/9)*x)}
z <- function(x){exp(-((x/(12/((log(2))^(1/3))))^3))}

integrate(x,lower=0,upper=Inf)$value - integrate(y,lower=0,upper=Inf)$value

up <- seq(from=20.64,to=20.65, by=.000001)
result <- rep(NA,length(up))

for (i in 1:length(up)){
   result[i] <- abs(integrate(y,lower=0,upper=up[i])$value - integrate(z,lower=0,upper=up[i])$value)
}
up[which.min(result)]

up <- 16.19104
abs(integrate(y,lower=0,upper=up)$value - integrate(z,lower=0,upper=up)$value)

set.seed(12)
#debug(simdat.cov)
data <- simdat.cov(n=1000,               #total number of observations
                   nevents=400,
                   accru=10,           # accrual rate per month
                   LTFU=0,            # loss to follow up proportion
                   pos=0.5,             # proportion of biomarker positive subjects
                   trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                   female=0.5,          # probability of being female
                   age.min=55,         # minimum age of subjects
                   age.max=85,         # maximum age of subjects
                   dis1=.2,            # probability of having disease severity 1
                   dis2=.2,            # probability of having disease severity 2
                   dis3=.3,            # probability of having disease severity 3
                   dis4=.2,            # probability of having disease severity 4
                   b0=12,              # the median survival for a male of minimum age with disease status 1
                   bstar=1.5,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=-2,              # the coefficient for age in the model to generate survival times
                   b2=3,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3,            # shape of weibull for trt B if PH=FALSE
                   clin.trial=FALSE
)
#saveRDS(data, file="Results/sim4/data.rds")
#data <- readRDS("Results/sim4type1/data.rds")
data <- data.frame(data)
chen.cox2(data=data,
          covs="female",
          tstar=12)
chen.cox(data=data,
          covs="female",
          tstar=12)




Y <- Surv(data$surv.time, data$event)

kmfit <- survfit(Y~data$trtn)

plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,30),
     xlab="Time in Months",ylab="Survival Probability",main="Scenario 4 Type 1 Error")
legend(22,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)
abline(v=20.2995, col="red")

t <- seq(20.298,20.3,by=.0001)
result <- rep(NA,length(t))
for (i in 1:length(t)){
    temp <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, tau=t[i])
    result[i] <- abs(temp$unadjusted.result[1,1])
}

t[which.min(result)]

data <- data.frame(data)
temp <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, tau=20.2995)



set.seed(12)
#debug(simdat.cov)
data <- simdat.cov(n=1000,               #total number of observations
                   nevents=400,
                   accru=10,           # accrual rate per month
                   LTFU=0.1,            # loss to follow up proportion
                   pos=0.5,             # proportion of biomarker positive subjects
                   trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                   female=0.5,          # probability of being female
                   age.min=55,         # minimum age of subjects
                   age.max=85,         # maximum age of subjects
                   dis1=.2,            # probability of having disease severity 1
                   dis2=.2,            # probability of having disease severity 2
                   dis3=.3,            # probability of having disease severity 3
                   dis4=.2,            # probability of having disease severity 4
                   b0=9,              # the median survival for a male of minimum age with disease status 1
                   bstar=3,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=-2,              # the coefficient for age in the model to generate survival times
                   b2=3,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3,            # shape of weibull for trt B if PH=FALSE
                   clin.trial=FALSE
)

library(foreign)
test <- data.frame(test)
test <- read.dta("ex.dta")

fit <- flexsurvspline(Surv(surv_time,event)~trtn+female+age_cent, anc=list(gamma1=~trtn), data=test, k=2, scale="hazard")
check <- RMSTdiff(data=test,
                  time.var="surv_time",
                  event.var="event",
                  trtn.var="trtn",
                  covs=c("female","age_cent"),
                  tmax=12,
                  method="RP31",
                  var.type2="delta")



write.table(test, file="ex.txt")
           
           
data <- data.frame(data)
temp <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn)
Y <- Surv(data$surv.time, data$event)

kmfit <- survfit(Y~data$trtn)

plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,30),
     xlab="Time in Months",ylab="Survival Probability",main="Scenario 3")
legend(22,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)


plot(temp$RMST.arm1$fit$time,temp$RMST.arm1$fit$surv)
lines(temp$RMST.arm0$fit$time,temp$RMST.arm0$fit$surv)


set.seed(12)
#debug(simdat.cov)
data <- simdat.cov(n=1000000,               #total number of observations
                   nevents=400,
                   accru=10,           # accrual rate per month
                   LTFU=0,            # loss to follow up proportion
                   pos=0.5,             # proportion of biomarker positive subjects
                   trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                   female=0.5,          # probability of being female
                   age.min=55,         # minimum age of subjects
                   age.max=85,         # maximum age of subjects
                   dis1=.2,            # probability of having disease severity 1
                   dis2=.2,            # probability of having disease severity 2
                   dis3=.3,            # probability of having disease severity 3
                   dis4=.2,            # probability of having disease severity 4
                   b0=12,              # the median survival for a male of minimum age with disease status 1
                   bstar=0,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=0,              # the coefficient for age in the model to generate survival times
                   b2=0,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3,            # shape of weibull for trt B if PH=FALSE
                   clin.trial=FALSE
)
saveRDS(data, file="Results/sim2/data.rds")
data <- readRDS("Results/sim2/data.rds")
data <- data.frame(data)
Y <- Surv(data$surv.time, data$event)

kmfit <- survfit(Y~data$trtn)

plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,30),
     xlab="Time in Months",ylab="Survival Probability",main="Scenario 2 Type 1 Error")
legend(22,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)
abline(v=20.61, col="red")


t <- seq(20.60,20.62,by=.001)
result <- rep(NA,length(t))
for (i in 1:length(t)){
     temp <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, tau=t[i])
     result[i] <- abs(temp$unadjusted.result[1,1])
}

t[which.min(result)]

data <- data.frame(data)
temp <- rmst2(time=data$surv.time, status=data$event, arm=data$trtn, tau=20.61)


test <- readRDS("Results/sim2/PH.ss400trteff3nocov.nocoef.LTFU0.02t")

covs <- NULL
formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
test <- coxph(formula,data=data)
p <- summary(test)$coefficients["trtn","Pr(>|z|)"]


b0 <- 18
b1 <- 2
b2 <- 3
b3 <- 0
b4 <- -1
b5 <- -3
b6 <- -6
b7 <- 0
bstar <- 0
bmark <- 0

n <- 1000000
shapeB <- 3
age.max <- 85
age.min <- 55
range <- age.max-age.min
lambda.age <- -log(.01)/range
age <- round(age.max - rexp(n, rate=lambda.age))
age.cent <- age - age.min
age.cent <- pmax(age.cent,0)
female <- rbinom(n=n, size=1, prob=0.5)
## getting each of the disease status variables
dis1p <- 0.2
dis2p <- 0.2
dis3p <- 0.3
dis4p <- 0.2
dis1 <- rbinom(n, size=1, prob=dis1p)
dis2 <- fifelse(dis1==0, dis2 <- rbinom(n, size=1, prob=(dis2p/(1-dis1p))), dis2 <- 0)
dis3 <- fifelse((dis1==0 & dis2==0), dis3 <- rbinom(n, size=1, prob=(dis3p/((1-dis1p)*(1-dis2p)))), dis3 <- 0)
dis4 <- fifelse((dis1==0 & dis2==0 & dis3==0), dis4 <- rbinom(1, size=1, prob=(dis4p/((1-dis1p)*(1-dis2p)*(1-dis3p)))), dis4 <- 0)
dis5 <- fifelse((dis1==0 & dis2==0 & dis3==0 & dis4==0), dis5 <- 1, dis5 <- 0)
trtn <- rbinom(n=n, size=1, prob=0.5)
marker.stat <- rbinom(n=n, size=1, prob=0.5)

med.surv <- b0 + b1*log((age.cent+1)) + b2*female + b3*dis2 + b4*dis3 + b5*dis4 +
    b6*dis5 + bstar*trtn + bmark*marker.stat + b7*marker.stat*trtn

event <- rep(1,1000000)

surv.time0 <- rexp(n=n, rate=log(2)/med.surv)


surv.time1 <- rweibull(n=n, shape=shapeB, scale=med.surv/(log(2)^(1/shapeB)))

surv0 <- Surv(time=surv.time0, event=event)
surv1 <- Surv(time=surv.time1, event=event)

plot(surv0)
lines(surv1)

















tmax <- NULL
data <- data.frame(data)
trt1 <- data[data$trtn==1,]
trt1 <- trt1[order(trt1$surv.time),]
trt0 <- data[data$trtn==0,]
trt0 <- trt0[order(trt0$surv.time),]
max1 <- max(trt1$surv.time)
max0 <- max (trt0$surv.time)
if ((max0 < max1) & trt0$event[nrow(trt0)]==1){
     tmax <- max1
}
if ((max0 > max1) & trt1$event[nrow(trt1)]==0){
     tmax <- max1
}
else {tmax <- max0}
covs <- "age"
undebug(chen.cox)
result <- chen.cox(data=data, covs=covs, tstar=tmax, p=rnd.prb1)



data <- data.frame(data)
covs <- NULL
formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
fit <- flexsurvspline(formula, anc=list(gamma1=~trtn), data=data, k=2, scale="hazard")
debug(surv.est.rp31)
result <- surv.est.rp31(data=data, tstar=12, fit=fit, alpha=.05)


#data <- data.frame(data)
debug(chen.cox)
### 1 run sample size 1000 took 2 hours
## 1 run ss 100 took 5 seconds
### 1 run ss 500 took 12 minutes
Sys.time()
test <- chen.cox(data, covs=c("age","female"), tstar=12)
Sys.time()

data <- readRDS("data_ex.RDS")
data2 <- data.frame(data)
debug(chen.cox2)
Sys.time()
test2 <- chen.cox2(data, covs=c("age","female"), tstar=12)
Sys.time()

test$var.est-test$hint1-test$hint0

test$hint1
test2$hint1

### checking empirical variance
set.seed(1001)
meandiff <- NULL
Sys.time()
for (i in 1:10){
     data <- simdat.cov(n=400,               #total number of observations
                        accru=10,           # accrual rate per month
                        LTFU=.02,            # loss to follow up proportion
                        pos=0.5,             # proportion of biomarker positive subjects
                        trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                        female=0.5,          # probability of being female
                        age.min=55,         # minimum age of subjects
                        age.max=85,         # maximum age of subjects
                        dis1=.1,            # probability of having disease severity 1
                        dis2=.2,            # probability of having disease severity 2
                        dis3=.3,            # probability of having disease severity 3
                        dis4=.2,            # probability of having disease severity 4
                        b0=12,              # the median survival for a male of minimum age with disease status 1
                        bstar=6,            # increase in median survival for being treated with B
                        bmark=0,             # coefficient for being biomarker positive
                        b1=0,              # the coefficient for age in the model to generate survival times
                        b2=0,                 # coefficient for being female
                        b3=0,              # coefficient for having disease status 2
                        b4=0,              # coefficient for having disease status 3
                        b5=0,              # coefficient for having disease status 4
                        b6=0,               # coefficient for having disease status 5
                        b7=0,               # coefficent for interaction between trt B and biomarker positive
                        PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                        shapeB=3            # shape of weibull for trt B if PH=FALSE
     )
     # debug(chen.cox)
     test <- chen.cox(data, covs=c("age","female"), tstar=12)
     meandiff[i] <- test$est
}
Sys.time()
## 1000 reps gives empirical variance 0.6555422; ss 100
## 1000 reps gives empirical variance 0.3206773; ss 200
var(meandiff)
## 1000 reps gives mean diff 0.9358072; ss 100
## 1000 reps gives mean diff 0.9491057; ss 200
mean(meandiff)

### using mike's code, ss 200, 10,000 reps, mean est was 0.9377548, truth is 0.9532
### empirical variance is 0.3125346
### estimated variance for 1 rep is 0.3239691


set.seed(1001)
data <- simdat.cov(n=200,               #total number of observations
                   accru=10,           # accrual rate per month
                   LTFU=.02,            # loss to follow up proportion
                   pos=0.5,             # proportion of biomarker positive subjects
                   trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                   female=0.5,          # probability of being female
                   age.min=55,         # minimum age of subjects
                   age.max=85,         # maximum age of subjects
                   dis1=.1,            # probability of having disease severity 1
                   dis2=.2,            # probability of having disease severity 2
                   dis3=.3,            # probability of having disease severity 3
                   dis4=.2,            # probability of having disease severity 4
                   b0=12,              # the median survival for a male of minimum age with disease status 1
                   bstar=6,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=0,              # the coefficient for age in the model to generate survival times
                   b2=0,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3            # shape of weibull for trt B if PH=FALSE
)
debug(chen.cox)
test1 <- chen.cox(data, covs=c("age","female"), tstar=12)
test2 <- chen.cox(data, covs=c("age","female"), tstar=12)

# start of code to reproduce Chen

## function to compute the cumulative hazard at given time t
for (i in 1:length(data[,1])){
     data$weights[i] <- length(data[data$surv.time==data$surv.time[i],1])
}
covs <- c("age","female")
ncovs <- length(covs)

data <- data[order(data$surv.time),]
trt <- data[data$trtn==1,]
notrt <- data[data$trtn==0,]

cov.form.chen <- 1
for (i in 1:length(covs)){
    
     cov.form.chen <- paste(covs[i],cov.form.chen,sep="+")
}
formula.chen <- as.formula(paste("Surv(surv.time, event)",cov.form.chen,sep="~"))

coxfit.trt <- coxph(formula.chen, data=trt, cluster=id)
coxfit.notrt <- coxph(formula.chen, data=notrt, cluster=id)

check2 <- survfit(coxfit.trt)

check <- unique(survfit(coxfit.trt, newdata=trt[1,])$surv)

summary(coxfit.trt)$coefficients
summary(coxfit.notrt)$coefficients

###############
tstar=12

# getting the combined covariate contribution for each subject for each cox model
temp.haz1 <- rep(0, length(data[,1]))
temp.haz0 <- rep(0, length(data[,1]))
for (j in 1:ncovs){
     temp.haz1 <- temp.haz1 + (summary(coxfit.trt)$coefficients[j,1]*data[,rownames(summary(coxfit.trt)$coefficients)[j]])
     temp.haz0 <- temp.haz0 + (summary(coxfit.notrt)$coefficients[j,1]*data[,rownames(summary(coxfit.notrt)$coefficients)[j]])
}
# the covariate contributions exponentiated
data$covhaz1 <- exp(temp.haz1)
data$covhaz0 <- exp(temp.haz0)



check0 <- test$surv.est0
check1 <- test$surv.est1

debug(sigmas)
test0 <- sigmas(data, covs=c("age","female"), tstar=12, trtn=0, surv.est=test$surv.est0)
test1 <- sigmas(data, covs=c("age","female"), tstar=12, trtn=1, surv.est=test$surv.est1)


debug(checn.cox)
test <- cox.cum.haz(data, covs=c("age","female"), tstar=12)


chk1 <- rnorm(1000)
tst1 <- chk1[1]
for (i in 2:length(chk1)){
     tst1 <- ((i-1)/i)*tst1 + (1/i)*chk1[i]
}
mean(chk1)
tst1

test <- data[,covs[1]]


#### start of FPM code implementation

fit <- flexsurvspline(Surv(surv.time, event)~female+age+trtn, anc=list(gamma1=~trtn),data=data.cov, k=2, scale="hazard")

undebug(surv.est.rp31)
test3 <- surv.est.rp31(data=data.cov,
                      tstar=12,
                      fit=fit)


undebug(sim.check)
test5 <- sim.check(1000)

emp.var <- var(test5$est, na.rm=T)
emp.var


x <- ldBounds(t=c(0.3,0.7,1), # the percentage of information times, last one always =1
              iuse=c(1,1),             # indicates O'brien fleming boundaries for both upper and lower boundaries
              alpha=c(0.1, 0.25))  



     