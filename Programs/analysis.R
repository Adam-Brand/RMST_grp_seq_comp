######## Analysis of results


sum.func <- function(result,events, int.time=3){
     
          bias.psu.sand <- rep(NA,length(result[,1]))
          bias.psu.ajk <- rep(NA,length(result[,1]))
          bias.km <- rep(NA,length(result[,1]))
          bias.rp.nai <- rep(NA,length(result[,1]))
          bias.rp.asy <- rep(NA,length(result[,1]))
     
          ### numbering the number of formal analyses; changing 0 to 3 when not rejecting
          for (i in 1:length(result[,1])){
               
               ### calculating bias which is est. - truth for the corresponding analysis
               if (result[i,31]==0 | result[i,31]==3){
                    bias.psu.sand[i] <- result[i,26] - result[i,12]
               }
               
               if (result[i,39]==0 | result[i,39]==3){
                    bias.psu.ajk[i] <-result[i,34] - result[i,12]
               }
               if (result[i,47]==0 | result[i,47]==3){
                    bias.km[i] <- result[i,42] - result[i,12]
               }
               if (result[i,55]==0 | result[i,55]==3){
                    bias.rp.nai[i] <- result[i,50] - result[i,12]
               }
               if (result[i,63]==0 | result[i,63]==3){
                    bias.rp.asy[i] <-result[i,58] - result[i,12]
               }
               
               ######################
               
               if (result[i,31]==1){
                    bias.psu.sand[i] <- result[i,26] - result[i,10]
               }
               
               if (result[i,39]==1){
                    bias.psu.ajk[i] <-result[i,34] - result[i,10]
               }
               if (result[i,47]==1){
                    bias.km[i] <- result[i,42] - result[i,10]
               }
               if (result[i,55]==1){
                    bias.rp.nai[i] <- result[i,50] - result[i,10]
               }
               if (result[i,63]==1){
                    bias.rp.asy[i] <-result[i,58] - result[i,10]
               }
               
               #################
               
               if (result[i,31]==2){
                    bias.psu.sand[i] <- result[i,26] - result[i,11]
               }
               
               if (result[i,39]==2){
                    bias.psu.ajk[i] <-result[i,34] - result[i,11]
               }
               if (result[i,47]==2){
                    bias.km[i] <- result[i,42] - result[i,11]
               }
               if (result[i,55]==2){
                    bias.rp.nai[i] <- result[i,50] - result[i,11]
               }
               if (result[i,63]==2){
                    bias.rp.asy[i] <-result[i,58] - result[i,11]
               }
          }
          
     
          ### taking column means
          res.mean <- colMeans(result, na.rm=TRUE)
          
          ### calculating rejection proportion for each method
          reject.psu.sand <- sum(result[,32], na.rm=TRUE)/sum(!is.na(result[,32]))
          reject.psu.ajk <- sum(result[,40], na.rm=TRUE)/sum(!is.na(result[,40]))
          reject.km <- sum(result[,48], na.rm=TRUE)/sum(!is.na(result[,48]))
          reject.rp.nai <- sum(result[,56], na.rm=TRUE)/sum(!is.na(result[,56]))
          reject.rp.asy <- sum(result[,64], na.rm=TRUE)/sum(!is.na(result[,64]))
          reject.ph <- sum(result[,23], na.rm=TRUE)/sum(!is.na(result[,23]))
          reject.lr <- sum(result[,20], na.rm=TRUE)/sum(!is.na(result[,20]))
          
          ### calculating mean and variance of the bias for each method
          mb.psu.sand <- mean(bias.psu.sand, na.rm=TRUE)
          vb.psu.sand <- var(bias.psu.sand, na.rm=TRUE)
          
          mb.psu.ajk <- mean(bias.psu.ajk, na.rm=TRUE)
          vb.psu.ajk <- var(bias.psu.ajk, na.rm=TRUE)
          
          mb.km <- mean(bias.km, na.rm=TRUE)
          vb.km <- var(bias.km, na.rm=TRUE)
          
          mb.rp.nai <- mean(bias.rp.nai, na.rm=TRUE)
          vb.rp.nai <- var(bias.rp.nai, na.rm=TRUE)
          
          mb.rp.asy <- mean(bias.rp.asy, na.rm=TRUE)
          vb.rp.asy <- var(bias.rp.asy, na.rm=TRUE)
          
          ### calculating average sample size for each method
          ss <- rep(NA,3)
          ss[1] <- ceiling(events*0.3)
          ss[2] <- ceiling(events*0.7)
          ss[3] <- events
          
          ss.psu.sand <- (ss[1]*sum(result[,31]==1) + ss[2]*sum(result[,31]==2) + ss[3]*sum(result[,31]==3 | result[,31]==0))/sum(!is.na(result[,31]))
          ss.psu.ajk <- (ss[1]*sum(result[,39]==1) + ss[2]*sum(result[,39]==2) + ss[3]*sum(result[,39]==3 | result[,39]==0))/sum(!is.na(result[,39]))
          ss.km <- (ss[1]*sum(result[,47]==1) + ss[2]*sum(result[,47]==2) + ss[3]*sum(result[,47]==3 | result[,47]==0))/sum(!is.na(result[,47]))
          ss.rp.nai <- (ss[1]*sum(result[,55]==1) + ss[2]*sum(result[,55]==2) + ss[3]*sum(result[,55]==3 | result[,55]==0))/sum(!is.na(result[,55]))
          ss.rp.asy <- (ss[1]*sum(result[,63]==1) + ss[2]*sum(result[,63]==2) + ss[3]*sum(result[,63]==3 | result[,63]==0))/sum(!is.na(result[,63]))
          ss.lr <- (ss[1]*sum(result[,19]==1) + ss[2]*sum(result[,19]==2) + ss[3]*sum(result[,19]==3 | result[,19]==0))/sum(!is.na(result[,19]))
          ss.hr <- (ss[1]*sum(result[,22]==1) + ss[2]*sum(result[,22]==2) + ss[3]*sum(result[,22]==3 | result[,22]==0))/sum(!is.na(result[,22]))
          
          
          if (ncol(result)==(55+6*int.time)){
               reject.cox <- sum(result[,72], na.rm=TRUE)/sum(!is.na(result[,64]))
               bias.cox <- NULL
               for (i in 1:length(result[,1])){
                    if (result[i,71]==0 | result[i,71]==3){
                         bias.cox[i] <- result[i,66] - result[i,12]
                    }
                    if (result[i,71]==1){
                         bias.cox[i] <- result[i,66] - result[i,10]
                    }
                    if (result[i,71]==2){
                         bias.cox[i] <- result[i,66] - result[i,11]
                    }
               }
               mb.cox <- mean(bias.cox, na.rm=TRUE)
               vb.cox <- var(bias.cox, na.rm=TRUE)
               ss.cox <- (ss[1]*sum(result[,71]==1) + ss[2]*sum(result[,71]==2) + ss[3]*sum(result[,71]==3 | result[,71]==0))/sum(!is.na(result[,71]))
          }
          
          
          if (ncol(result)==(47+6*int.time)){
               fin <- c("lr.rej"=reject.lr, "lr.ss"=ss.lr,
                        
                        "hr.rej"=reject.ph, "hr.ss"=ss.hr,
                    
                        "psu.sand.rej"=reject.psu.sand,"psu.sand.bias"=mb.psu.sand,
                        "psu.sand.vbias"=vb.psu.sand,"psu.sand.ss"=ss.psu.sand,
                        
                        "psu.ajk.rej"=reject.psu.ajk,"psu.ajk.bias"=mb.psu.ajk,
                        "psu.ajk.vbias"=vb.psu.ajk,"psu.ajk.ss"=ss.psu.ajk,
                        
                        "km.rej"=reject.km,"km.bias"=mb.km,
                        "km.vbias"=vb.km,"km.ss"=ss.km,
                        
                        "rp.nai.rej"=reject.rp.nai,"rp.nai.bias"=mb.rp.nai,
                        "rp.nai.vbias"=vb.rp.nai,"rp.nai.ss"=ss.rp.nai,
                        
                        "rp.asy.rej"=reject.rp.asy,"rp.asy.bias"=mb.rp.asy,
                        "rp.asy.vbias"=vb.rp.asy,"rp.asy.ss"=ss.rp.asy
                        )
          }
          
          else if (ncol(result)==(55+6*int.time)){
               fin <- c("lr.rej"=reject.lr, "lr.ss"=ss.lr,
                        
                        "hr.rej"=reject.ph, "hr.ss"=ss.hr,
                        
                        "psu.sand.rej"=reject.psu.sand,"psu.sand.bias"=mb.psu.sand,
                        "psu.sand.vbias"=vb.psu.sand,"psu.sand.ss"=ss.psu.sand,
                        
                        "psu.ajk.rej"=reject.psu.ajk,"psu.ajk.bias"=mb.psu.ajk,
                        "psu.ajk.vbias"=vb.psu.ajk,"psu.ajk.ss"=ss.psu.ajk,
                        
                        "km.rej"=reject.km,"km.bias"=mb.km,
                        "km.vbias"=vb.km,"km.ss"=ss.km,
                        
                        "rp.nai.rej"=reject.rp.nai,"rp.nai.bias"=mb.rp.nai,
                        "rp.nai.vbias"=vb.rp.nai,"rp.nai.ss"=ss.rp.nai,
                        
                        "rp.asy.rej"=reject.rp.asy,"rp.asy.bias"=mb.rp.asy,
                        "rp.asy.vbias"=vb.rp.asy,"rp.asy.ss"=ss.rp.asy,
                        
                        "cox.rej"=reject.cox,"cox.bias"=mb.cox,
                        "cox.vbias"=vb.cox,"cox.ss"=ss.cox
                        )
          }
          
          return(fin)
}



#### for sim1 nocov, 550; top half of table for scenario 1
sim1nocov <- readRDS("Results/sim1grpseq/PH.ss550trteff3nocov.nocoef.LTFU0.02times.rds")
#debug(sum.func)
result1top <- sum.func(sim1nocov, events=550)
result1top
#### the RP method failed to converge in 8 trials (1%)
summary(sim1nocov)
# type 1 error row
sim1type <- readRDS("Results/sim1grpseqtype1/PH.ss550trteff0nocov.nocoef.LTFU0.02times.rds")
result1t1 <- sum.func(sim1type, events=550)
result1t1
### average analysis times for scenario 1
colMeans(sim1nocov[,1:3])


#### for sim1 nocov, 550; bottom half of table for scenario 1
sim1nocov <- readRDS("Results/sim1grpseq/PH.ss550trteff3allcov.nocoef.LTFU0.02times.rds")
#debug(sum.func)
result1bot <- sum.func(sim1nocov, events=550)
result1bot
#### the RP method failed to converge in 9 trials (1%)
summary(sim1nocov)
# type 1 error row
sim1type <- readRDS("Results/sim1grpseqtype1/PH.ss550trteff0allcov.nocoef.LTFU0.02times.rds")
result1b1 <- sum.func(sim1type, events=550)
result1b1


#### for sim2 nocov, 175; top half of table for scenario 2
sim2nocov <- readRDS("Results/sim2grpseq/PH.ss175trteff3nocov.fem.age.LTFU0.02times.rds")
#debug(sum.func)
result2top <- sum.func(sim2nocov, events=175)
result2top
#### the RP method failed to converge in 2 trials (1%)
summary(sim2nocov)
# type 1 error row
sim2type <- readRDS("Results/sim2grpseqtype1/PH.ss175trteff0nocov.fem.age.LTFU0.02times.rds")
result2t1 <- sum.func(sim2type, events=175)
result2t1
### average analysis times for scenario 2
colMeans(sim2nocov[,1:3])

#### for sim2 nocov, 175; bottom half of table for scenario 2
sim2allcov <- readRDS("Results/sim2grpseq/PH.ss175trteff3allcov.age.fem.LTFU0.02times.rds")
#debug(sum.func)
result2bot <- sum.func(sim2allcov, events=175)
result2bot
#### the RP method failed to converge in 11 trials (1%)
summary(sim2allcov)
# type 1 error row
sim2type <- readRDS("Results/sim2grpseqtype1/PH.ss175trteff0allcov.age.fem.LTFU0.02times.rds")
result2b1 <- sum.func(sim2type, events=175)
result2b1



#### for sim3 nocov, 220; top half of table for scenario 3
sim3nocov <- readRDS("Results/sim3grpseq/PH.ss220trteff2nocov.fem.age.LTFU0.02times.rds")
#debug(sum.func)
result3top <- sum.func(sim3nocov, events=220)
result3top
#### the RP method failed to converge in 2 trials (1%)
summary(sim3nocov)
# type 1 error row
sim3type <- readRDS("Results/sim3grpseqtype1/PH.ss220trteff0nocov.fem.age.LTFU0.02times.rds")
result <- sum.func(sim3type, events=220)
result
### average analysis times for scenario 3
colMeans(sim3nocov[,1:3])


#### for sim3 nocov, 220; bottom half of table for scenario 3
sim3allcov <- readRDS("Results/sim3grpseq/PH.ss220trteff2allcov.age.fem.LTFU0.02times.rds")
result3bot <- sum.func(sim3allcov, events=220)
result3bot
#### the RP method failed to converge in 2 trials (1%)
summary(sim3nocov)
# type 1 error row
sim3type <- readRDS("Results/sim3grpseqtype1/PH.ss220trteff0allcov.age.fem.LTFU0.02times.rds")
result <- sum.func(sim3type, events=220)
result



#### for sim4 nocov, 150; top half of table for scenario 4
sim4nocov <- readRDS("Results/sim4grpseq/PH.ss150trteff6nocov.fem.age.LTFU0.02times.rds")
#debug(sum.func)
result4top <- sum.func(sim4nocov, events=150)
result4top
#### the RP method failed to converge in 2 trials (1%)
summary(sim4nocov)
# type 1 error row
sim4type <- readRDS("Results/sim4grpseqtype1/PH.ss150trteff0nocov.fem.age.LTFU0.02times.rds")
result <- sum.func(sim4type, events=150)
result
### average analysis times for scenario 4
colMeans(sim4nocov[,1:3])


#### for sim4 allcov, 150; bottom half of table for scenario 4
sim4allcov <- readRDS("Results/sim4grpseq/PH.ss150trteff6allcov.age.fem.LTFU0.02times.rds")
#debug(sum.func)
result4bot <- sum.func(sim4allcov, events=150)
result4bot
#### the RP method failed to converge in 2 trials (1%)
summary(sim4allcov)
# type 1 error row
sim4type <- readRDS("Results/sim4grpseqtype1/PH.ss150trteff0allcov.age.fem.LTFU0.02times.rds")
result <- sum.func(sim4type, events=150)
result




