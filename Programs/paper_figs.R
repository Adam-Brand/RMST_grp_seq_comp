####### code to produce the figure for the paper Implementing RMST in R

source("Programs/source.R")


### figure for scenario 1
curve(exp(-(log(2)/12)*x), from=0,to=60, col="green",lty=2, xlab="Time in Months",
      ylab="Survival Probability", main="Scenario 1")
curve(exp(-(log(2)/9)*x), from=0,to=60, col="blue",lty=1, add=TRUE)
legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)

### figure for scenario 2

set.seed(12)
data <- simdat.cov(n=1000000,               #total number of observations
                    nevents=175,
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
data <- data.frame(data)

Y <- Surv(data$surv.time, data$event)
kmfit <- survfit(Y~data$trtn)
plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,60), xlab="Time in Months",
     ylab="Survival Probability", main="Scenario 2")
legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)


## figure for scenario 3

set.seed(12)

data3 <- simdat.cov(n=1000000,               #total number of observations
                    nevents=300,
                    accru=10,           # accrual rate per month
                    LTFU=0,            # loss to follow up proportion
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
                    shapeB=0.5,            # shape of weibull for trt B if PH=FALSE
                    clin.trial=FALSE
)
data3 <- data.frame(data3)

Y3 <- Surv(data3$surv.time, data3$event)
kmfit3 <- survfit(Y3~data3$trtn)
plot(kmfit3, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,60), xlab="Time in Months",
     ylab="Survival Probability", main="Scenario 3")
legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)


## figure for scenario 4

set.seed(12)

data4 <- simdat.cov(n=1000000,               #total number of observations
                    nevents=300,
                    accru=10,           # accrual rate per month
                    LTFU=0,            # loss to follow up proportion
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
                    bstar=6,            # increase in median survival for being treated with B
                    bmark=0,             # coefficient for being biomarker positive
                    b1=-1,              # the coefficient for age in the model to generate survival times
                    b2=3,                 # coefficient for being female
                    b3=0,              # coefficient for having disease status 2
                    b4=0,              # coefficient for having disease status 3
                    b5=0,              # coefficient for having disease status 4
                    b6=0,               # coefficient for having disease status 5
                    b7=0,               # coefficent for interaction between trt B and biomarker positive
                    PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                    shapeB=2,            # shape of weibull for trt B if PH=FALSE
                    clin.trial=FALSE
)
data4 <- data.frame(data4)

Y4 <- Surv(data4$surv.time, data4$event)
kmfit4 <- survfit(Y4~data4$trtn)
plot(kmfit4, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,60), xlab="Time in Months",
     ylab="Survival Probability", main="Scenario 4")
legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)
















