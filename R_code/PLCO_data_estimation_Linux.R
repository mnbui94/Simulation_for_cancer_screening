#author : Mai Ngoc Bui
#parameter estimation (PLCO data - Linux version)
require(msm)
require(ucminf)
require(pbapply)
require(readxl)
require(plyr)
require(parallel)
#### FUNCTIONS ####
age_centered <- 68
setwd("~/Biometrical/C++")
dyn.load("Lik_grid_Linux.so")#load built-in C++ code


source("~/Biometrical/R_code/Estimation_functions.R")

##### MODEL A #####
#### 1.four-states models defined by R. Bhatt ####
setwd("~/Biometrical/PLCO_data")
load("ModelA_4States_PLCOdata.RData")#filename is full_data
#total no of men in each arm
length(unique(full_data[full_data$arm == 1,]$id))
length(unique(full_data[full_data$arm == 2,]$id))
#printing state table
print(statetable.msm(state, id, full_data[full_data$arm == 1,]))
print(statetable.msm(state, id, full_data[full_data$arm == 2,]))


#separate data into similar sub-group to assign to each core
subj <- list()
subj[[1]] <- names(table(full_data[full_data$arm ==1,]$id))
subj[[2]] <- names(table(full_data[full_data$arm ==2,]$id))

no_cores = 15#number of cores
index<- list()
index[[1]] <- splitIndices(length(subj[[1]]),no_cores)
index[[2]] <- splitIndices(length(subj[[2]]),no_cores)

mul_data <- list()
mul_no_obs <- list()
for(i in 1:no_cores){
  mul_data[[i]] <- full_data[full_data$id %in% c(subj[[1]][index[[1]][[i]]],
                                                 subj[[2]][index[[2]][[i]]]),]
  mul_data[[i]] <- mul_data[[i]][order(mul_data[[i]]$id),]
  mul_no_obs[[i]] <- c(0,cumsum(count(mul_data[[i]],"id")$freq))
}


#initialise values for parameters, e.g. p0 and p1.
p0 = c(-4.47,-4.34,-1.47,0.04,0.09,-0.59,-0.16)
p1 = c(-4,-4,-0.2,0.04,0.09,-2,0.1)

system.time(print(lik_grid_4States_parallel(mul_no_obs,p0,mul_data,no_cores)))#takes 0.7 sec, result = 212835
system.time(print(lik_grid_4States_parallel(mul_no_obs,p1,mul_data,no_cores)))#takes 0.7 sec, result = 217029.6

system.time(
  est1_ModelA_4states <- ucminf(p0, function(p) lik_grid_4States_parallel(mul_no_obs,p,mul_data,no_cores), control = list(maxeval = 50000),hessian = TRUE)
) #takes 343 sec = 5.7 min
print(est1_ModelA_4states$par)#estimated parameters
print(est1_ModelA_4states$value)#likelihood values = 212409.6

system.time(
  est2_ModelA_4states <- ucminf(p1, function(p) lik_grid_4States_parallel(mul_no_obs,p,mul_data,no_cores), control = list(maxeval = 50000),hessian = TRUE)
) #takes 336 sec = 5.6 min
print(est2_ModelA_4states$par)#estimated parameters
print(est2_ModelA_4states$value)#likelihood values = 212409.6


setwd("~/Biometrical/Estimation_results")
save(est1_ModelA_4states, file = "est1_ModelA_4states.RData")

#### 2.three-states (exiting from state 3) ####
setwd("~/Biometrical/PLCO_data")
load("ModelA_3States_PLCOdata.RData")

#reset the state 3,4 and 5 to 1,2 and 3.
data <- data_state3
data$state[data$state == 3] <- 1
data$state[data$state == 4] <- 2
data$state[data$state == 5] <- 3

no_obs <- c(0,cumsum(count(data,"id")$freq))
data <- data[order(data$id),]

p0 <- c( -3.648022, -1.955821, -0.05, -0.1)
p1 <- c(-4.121640, -3.951800, 0.1, 0.05)
lik_fixed_3States(no_obs, p0,data)# result = 13830.35
lik_fixed_3States(no_obs, p1,data)#result = 8719.419

est1_ModelA_3States <- ucminf(p0, function(p) lik_fixed_3States(no_obs,p,data), control = list(maxeval = 10000),hessian = TRUE)
est1_ModelA_3States$par #estimated parameters
est1_ModelA_3States$value #likelihood values = 7942.418

est2_ModelA_3States <- ucminf(p1, function(p) lik_fixed_3States(no_obs,p,data), control = list(maxeval = 10000),hessian = TRUE)
est2_ModelA_3States$par #estimated parameters
est2_ModelA_3States$value #likelihood values = 7942.418

p    <- est1_ModelA_3States$par
p.se <- sqrt(diag(solve(est1_ModelA_3States$hessian)))

print(cbind("paramaters" = c("lambda_{34}","lambda_{35}","beta_{34}","beta_{35}"), "estimated" =round(p,2),
            "SE" =round(p.se,2),"lower(95% CI)" = round(p -1.96*p.se,2),"upper(95% CI)" = round(p+1.96*p.se,2),
            "p-values"=round(1-pchisq((p/p.se)^2,df=1),2)),quote=FALSE)


setwd("~/Biometrical/Estimation_results")
save(est1_ModelA_3States, file = "est1_ModelA_3States.RData")


#### MODEL B ####
setwd("~/Biometrical/PLCO_data")
load("ModelB_3States_PLCOdata.RData")

#reset the state 2,4 and 5 to 1,2 and 3.
data <- data_state2
data$state[data$state == 2] <- 1
data$state[data$state == 4] <- 2
data$state[data$state == 5] <- 3

no_obs <- c(0,cumsum(count(data,"id")$freq))
data <- data[order(data$id),]

p0 <- c( -3.648022, -1.955821, -0.05, -0.1)
p1 <- c(-4.121640, -3.951800, 0.1, 0.05)
lik_fixed_3States(no_obs, p0,data) #result = 8403.907
lik_fixed_3States(no_obs, p1,data) #result = 4721.867

est1_ModelB_3States <- ucminf(p0, function(p) lik_fixed_3States(no_obs,p,data), control = list(maxeval = 10000),hessian = TRUE)
est1_ModelB_3States$par #estimated parameters
est1_ModelB_3States$value #likelihood values = 4265.325

est2_ModelB_3States <- ucminf(p1, function(p) lik_fixed_3States(no_obs,p,data), control = list(maxeval = 10000),hessian = TRUE)
est2_ModelB_3States$par #estimated parameters
est2_ModelB_3States$value #likelihood values = 4265.325

p    <- est1_ModelB_3States$par
p.se <- sqrt(diag(solve(est1_ModelB_3States$hessian)))

print(cbind("paramaters" = c("lambda^*_{24}","lambda^*_{25}","beta^*_{24}","beta^*_{25}"), "estimated" =round(p,2),
            "SE" =round(p.se,2),"lower(95% CI)" = round(p -1.96*p.se,2),"upper(95% CI)" = round(p+1.96*p.se,2),
            "p-values"=round(1-pchisq((p/p.se)^2,df=1),2)),quote=FALSE)


setwd("~/Biometrical/Estimation_results")
save(est1_ModelB_3States, file = "est1_ModelB_3States.RData")