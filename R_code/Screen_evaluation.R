#author : Mai Ngoc Bui
#simulation for the screened group (SCREEENING EVALUATION)

start_time <- Sys.time()
setwd("~/Biometrical/C++_code")
dyn.load("MSM_simulation_Windows.dll")

source("~/Biometrical/R_code/Simulation_functions.R")


#### MSM simulation ####
age_centered <- 68 #centering age
base_age <- 40

##parameters required for simulating control group.
par <- data.frame("lambda" = c(-4.95,-4.4,-2.5,-0.07,-2.13),
                  "beta" = c(0.02, 0.12, 0.00, 0.01, 0.17))#parameters from model A.
study_years = c(2000,2035) #year period of the study.
prop_year = 1 #proportion of each base_age.
left_trunc = 40 #left-truncation age.
#Gompertz distribution for all transitions in model A, except transition from state 2 to state 3.
hazard_distr <- c(2,2,0,2,2)


##parameters required for simulating screening group.
sub_par <- data.frame("lambda" = c(-2.4,-5.1),
                      "beta" = c(0.06, 0.09)) #parameters from model B.
screen_misc <- c(-1.61,-0.12) #misclassification parameters.
100/(1+exp(screen_misc[1] + screen_misc[2] * (69 - age_centered)))#percentages of sensitivity at age 69.
100/(1+exp(screen_misc[1] + screen_misc[2] * (55 - age_centered)))#percentages of sensitivity at age 55.
take_up = 1 #compliance rate.
#Gompertz distribution for all transitions in model B.
hazard_sub_distr <- c(2,2)



N <- 5*10^5 #population size.


data_arm2 <- control_sim_Cpp_Rseed(par,  study_years, base_age, prop_year, left_trunc, id = 1:N, hazard_distr, max_no_seeds = 2)

#### FREQ = 2 ####
screen_freq  = 2 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
Freq2_lead_time <- list()
for (j in 1:3){
  cat(paste("*** Current run, j = ",j,"*** \n"))
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)


  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  Freq2_lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  output_rep[[j]] <- list()
  #CANCER MORTALITY REDUCTION#
  output_rep[[j]][[1]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$life_time_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$life_time_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$life_time_record$state == 5) / sum(data_arm2$life_time_record$state == 5)))
  names(output_rep[[j]])[1] <- "cancer_deaths_reduction_life_time"

  output_rep[[j]][[2]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$study_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$study_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$study_record$state == 5) / sum(data_arm2$study_record$state == 5)))
  names(output_rep[[j]])[2] <- "cancer_deaths_reduction_during_study"



  #EXCESS INCIDENCE#

  output_rep[[j]][[3]] <- data.frame("total_cancer_arm1"= sum(data_arm1$life_time_record$state %in% c(2,3)),
                                     "total_cancer_arm2"= sum(data_arm2$life_time_record$state %in% c(3)))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]],
                                "excess_cancer_cases" = output_rep[[j]][[3]]$total_cancer_arm1 - output_rep[[j]][[3]]$total_cancer_arm2,"total_screen_detected" = sum(data_arm1$life_time_record$state == 2))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]], "percentages" = 100 * output_rep[[j]][[3]]$excess_cancer_cases/ output_rep[[j]][[3]]$total_screen_detected)
  names(output_rep[[j]])[3] <- "excess_incidence_life_time"

  #INTERVAL CANCER#

  output_rep[[j]][[4]] <- sum((data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] > screen_age[1]) & (data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] < screen_age[2]))
  names(output_rep[[j]])[4] <- "total_interval_cancers"


}
Freq2 <- output_rep

#### FREQ = 8 ####
screen_freq  = 8 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
Freq8_lead_time <- list()
for (j in 1:3){
  cat(paste("*** Current run, j = ",j,"*** \n"))
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  Freq8_lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  output_rep[[j]] <- list()
  #CANCER MORTALITY REDUCTION#
  output_rep[[j]][[1]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$life_time_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$life_time_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$life_time_record$state == 5) / sum(data_arm2$life_time_record$state == 5)))
  names(output_rep[[j]])[1] <- "cancer_deaths_reduction_life_time"

  output_rep[[j]][[2]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$study_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$study_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$study_record$state == 5) / sum(data_arm2$study_record$state == 5)))
  names(output_rep[[j]])[2] <- "cancer_deaths_reduction_during_study"



  #EXCESS INCIDENCE#

  output_rep[[j]][[3]] <- data.frame("total_cancer_arm1"= sum(data_arm1$life_time_record$state %in% c(2,3)),
                                     "total_cancer_arm2"= sum(data_arm2$life_time_record$state %in% c(3)))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]],
                                "excess_cancer_cases" = output_rep[[j]][[3]]$total_cancer_arm1 - output_rep[[j]][[3]]$total_cancer_arm2,"total_screen_detected" = sum(data_arm1$life_time_record$state == 2))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]], "percentages" = 100 * output_rep[[j]][[3]]$excess_cancer_cases/ output_rep[[j]][[3]]$total_screen_detected)
  names(output_rep[[j]])[3] <- "excess_incidence_life_time"

  #INTERVAL CANCER#

  output_rep[[j]][[4]] <- sum((data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] > screen_age[1]) & (data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] < screen_age[2]))
  names(output_rep[[j]])[4] <- "total_interval_cancers"

}
Freq8 <- output_rep

#### FREQ = 16 ####
screen_freq  = 16 #frequency of screening.
screen_age = c(55,69) #range of screening age.
Freq16_lead_time <- list()
output_rep <- list()
for (j in 1:3){
  cat(paste("*** Current run, j = ",j,"*** \n"))
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  Freq16_lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  output_rep[[j]] <- list()
  #CANCER MORTALITY REDUCTION#
  output_rep[[j]][[1]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$life_time_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$life_time_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$life_time_record$state == 5) / sum(data_arm2$life_time_record$state == 5)))
  names(output_rep[[j]])[1] <- "cancer_deaths_reduction_life_time"

  output_rep[[j]][[2]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$study_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$study_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$study_record$state == 5) / sum(data_arm2$study_record$state == 5)))
  names(output_rep[[j]])[2] <- "cancer_deaths_reduction_during_study"



  #EXCESS INCIDENCE#

  output_rep[[j]][[3]] <- data.frame("total_cancer_arm1"= sum(data_arm1$life_time_record$state %in% c(2,3)),
                                     "total_cancer_arm2"= sum(data_arm2$life_time_record$state %in% c(3)))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]],
                                "excess_cancer_cases" = output_rep[[j]][[3]]$total_cancer_arm1 - output_rep[[j]][[3]]$total_cancer_arm2,"total_screen_detected" = sum(data_arm1$life_time_record$state == 2))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]], "percentages" = 100 * output_rep[[j]][[3]]$excess_cancer_cases/ output_rep[[j]][[3]]$total_screen_detected)
  names(output_rep[[j]])[3] <- "excess_incidence_life_time"

  #INTERVAL CANCER#

  output_rep[[j]][[4]] <- sum((data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] > screen_age[1]) & (data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] < screen_age[2]))
  names(output_rep[[j]])[4] <- "total_interval_cancers"

}
Freq16 <- output_rep

#### Age 55 - 69 ####
screen_freq  = 4 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
Freq4_lead_time <- list()
for (j in 1:3){
  cat(paste("*** Current run, j = ",j,"*** \n"))
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  Freq4_lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  output_rep[[j]] <- list()
  #CANCER MORTALITY REDUCTION#
  output_rep[[j]][[1]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$life_time_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$life_time_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$life_time_record$state == 5) / sum(data_arm2$life_time_record$state == 5)))
  names(output_rep[[j]])[1] <- "cancer_deaths_reduction_life_time"

  output_rep[[j]][[2]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$study_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$study_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$study_record$state == 5) / sum(data_arm2$study_record$state == 5)))
  names(output_rep[[j]])[2] <- "cancer_deaths_reduction_during_study"



  #EXCESS INCIDENCE#

  output_rep[[j]][[3]] <- data.frame("total_cancer_arm1"= sum(data_arm1$life_time_record$state %in% c(2,3)),
                                     "total_cancer_arm2"= sum(data_arm2$life_time_record$state %in% c(3)))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]],
                                "excess_cancer_cases" = output_rep[[j]][[3]]$total_cancer_arm1 - output_rep[[j]][[3]]$total_cancer_arm2,"total_screen_detected" = sum(data_arm1$life_time_record$state == 2))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]], "percentages" = 100 * output_rep[[j]][[3]]$excess_cancer_cases/ output_rep[[j]][[3]]$total_screen_detected)
  names(output_rep[[j]])[3] <- "excess_incidence_life_time"

  #INTERVAL CANCER#

  output_rep[[j]][[4]] <- sum((data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] > screen_age[1]) & (data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] < screen_age[2]))
  names(output_rep[[j]])[4] <- "total_interval_cancers"

}
Freq4 <- output_rep

#### Age 50 - 64 ####
screen_age = c(50,64) #range of screening age.
screen_freq  = 4 #frequency of screening.
output_rep <- list()
Freq4_age50_lead_time <- list()
for (j in 1:3){
  cat(paste("*** Current run, j = ",j,"*** \n"))
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)


  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  Freq4_age50_lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  output_rep[[j]] <- list()
  #CANCER MORTALITY REDUCTION#
  output_rep[[j]][[1]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$life_time_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$life_time_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$life_time_record$state == 5) / sum(data_arm2$life_time_record$state == 5)))
  names(output_rep[[j]])[1] <- "cancer_deaths_reduction_life_time"

  output_rep[[j]][[2]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$study_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$study_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$study_record$state == 5) / sum(data_arm2$study_record$state == 5)))
  names(output_rep[[j]])[2] <- "cancer_deaths_reduction_during_study"



  #EXCESS INCIDENCE#

  output_rep[[j]][[3]] <- data.frame("total_cancer_arm1"= sum(data_arm1$life_time_record$state %in% c(2,3)),
                                     "total_cancer_arm2"= sum(data_arm2$life_time_record$state %in% c(3)))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]],
                                "excess_cancer_cases" = output_rep[[j]][[3]]$total_cancer_arm1 - output_rep[[j]][[3]]$total_cancer_arm2,"total_screen_detected" = sum(data_arm1$life_time_record$state == 2))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]], "percentages" = 100 * output_rep[[j]][[3]]$excess_cancer_cases/ output_rep[[j]][[3]]$total_screen_detected)
  names(output_rep[[j]])[3] <- "excess_incidence_life_time"

  #INTERVAL CANCER#

  output_rep[[j]][[4]] <- sum((data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] > screen_age[1]) & (data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] < screen_age[2]))
  names(output_rep[[j]])[4] <- "total_interval_cancers"


}
Freq4_age50 <- output_rep

#### Age 60 - 74 ####
screen_age = c(60,74) #range of screening age.
screen_freq  = 4 #frequency of screening.
output_rep <- list()
Freq4_age60_lead_time <- list()
for (j in 1:3){
  cat(paste("*** Current run, j = ",j,"*** \n"))
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)


  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  Freq4_age60_lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  output_rep[[j]] <- list()
  #CANCER MORTALITY REDUCTION#
  output_rep[[j]][[1]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$life_time_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$life_time_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$life_time_record$state == 5) / sum(data_arm2$life_time_record$state == 5)))
  names(output_rep[[j]])[1] <- "cancer_deaths_reduction_life_time"

  output_rep[[j]][[2]] <-   data.frame("cancer_deaths_screen" = sum(data_arm1$study_record$state == 5),
                                       "cancer_deaths_control" = sum(data_arm2$study_record$state == 5),
                                       "percentages" = 100 *(1 - sum(data_arm1$study_record$state == 5) / sum(data_arm2$study_record$state == 5)))
  names(output_rep[[j]])[2] <- "cancer_deaths_reduction_during_study"



  #EXCESS INCIDENCE#

  output_rep[[j]][[3]] <- data.frame("total_cancer_arm1"= sum(data_arm1$life_time_record$state %in% c(2,3)),
                                     "total_cancer_arm2"= sum(data_arm2$life_time_record$state %in% c(3)))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]],
                                "excess_cancer_cases" = output_rep[[j]][[3]]$total_cancer_arm1 - output_rep[[j]][[3]]$total_cancer_arm2,"total_screen_detected" = sum(data_arm1$life_time_record$state == 2))
  output_rep[[j]][[3]] <- cbind(output_rep[[j]][[3]], "percentages" = 100 * output_rep[[j]][[3]]$excess_cancer_cases/ output_rep[[j]][[3]]$total_screen_detected)
  names(output_rep[[j]])[3] <- "excess_incidence_life_time"

  #INTERVAL CANCER#

  output_rep[[j]][[4]] <- sum((data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] > screen_age[1]) & (data_arm1$life_time_record$age[data_arm1$life_time_record$state == 3] < screen_age[2]))
  names(output_rep[[j]])[4] <- "total_interval_cancers"

}
Freq4_age60 <- output_rep



####compare with screen_evaluation_long####
#compare your runs with the provided outputs from the author.
setwd("~/Biometrical/Figures_Tables")
sink("Compare_screen_evaluation.txt")

#freq = 2
cat(paste("#### Your run of screen_freq = 2: \n" ))
cat(paste("\n"))
print(Freq2[[1]])
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 2 from Intermediate_results folder: \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq2_screen.RData")
cat(paste("\n"))
print(output_rep[[1]])
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))


#freq = 4
cat(paste("#### Your run of screen_freq = 4: \n" ))
cat(paste("\n"))
print(Freq4[[1]])
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 4 from Intermediate_results folder: \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq4_age55_screen.RData")
cat(paste("\n"))
print(output_rep[[1]])
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))

#freq = 8
cat(paste("#### Your run of screen_freq = 8: \n" ))
cat(paste("\n"))
print(Freq8[[1]])
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 8 from Intermediate_results folder: \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq8_screen.RData")
cat(paste("\n"))
print(output_rep[[1]])
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))
#freq = 16
cat(paste("#### Your run of screen_freq = 16: \n" ))
cat(paste("\n"))
print(Freq16[[1]])
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 16 from Intermediate_results folder: \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq16_screen.RData")
cat(paste("\n"))
print(output_rep[[1]])
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))



#age 50-64
cat(paste("#### Your run of screen_freq = 4, age 50 - 64: \n" ))
cat(paste("\n"))
print(Freq4_age50[[1]])
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 4, age 50 - 64 from Intermediate_results folder: \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq4_age50_screen.RData")
cat(paste("\n"))
print(output_rep[[1]])
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))

#age 60-74
cat(paste("#### Your run of screen_freq = 4, age 50 - 64: \n" ))
cat(paste("\n"))
print(Freq4_age60[[1]])
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 4, age 50 - 64 from Intermediate_results folder: \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq4_age60_screen.RData")
cat(paste("\n"))
print(output_rep[[1]])
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))
sink()




####compare with screen_evaluation_lead_times####
#compare your runs with the provided outputs from the author.
setwd("~/Biometrical/Figures_Tables")
sink("Compare_screen_evaluation_lead_times.txt")

#freq = 2
cat(paste("#### Your run of screen_freq = 2: Statistic summary \n" ))
cat(paste("\n"))
print(summary(Freq2_lead_time[[1]]))
cat(paste("\n\n\n"))

cat(paste("#### First output of  screen_freq = 2 from Intermediate_results folder: Statistic summary \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq2_screen_lead_time.RData")
cat(paste("\n"))
print(summary(lead_time[[1]]))
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))


#freq = 4
cat(paste("#### Your run of screen_freq = 4: Statistic summary \n" ))
cat(paste("\n"))
print(summary(Freq4_lead_time[[1]]))
cat(paste("\n\n\n"))

cat(paste("#### First output of  screen_freq = 4 from Intermediate_results folder: Statistic summary \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq4_age55_screen_lead_time.RData")
cat(paste("\n"))
print(summary(lead_time[[1]]))
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))

#freq = 8
cat(paste("#### Your run of screen_freq = 8: Statistic summary \n" ))
cat(paste("\n"))
print(summary(Freq8_lead_time[[1]]))
cat(paste("\n\n\n"))

cat(paste("#### First output of  screen_freq = 8 from Intermediate_results folder: Statistic summary \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq8_screen_lead_time.RData")
cat(paste("\n"))
print(summary(lead_time[[1]]))
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))
#freq = 16
cat(paste("#### Your run of screen_freq = 16: Statistic summary \n" ))
cat(paste("\n"))
print(summary(Freq16_lead_time[[1]]))
cat(paste("\n\n\n"))

cat(paste("#### First output of  screen_freq = 16 from Intermediate_results folder: Statistic summary \n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq16_screen_lead_time.RData")
cat(paste("\n"))
print(summary(lead_time[[1]]))
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))



#age 50-64
cat(paste("#### Your run of screen_freq = 4, age 50 - 64: Statistic summary \n" ))
cat(paste("\n"))
print(summary(Freq4_age50_lead_time[[1]]))
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 4, age 50 - 64 from Intermediate_results folder: Statistic summary\n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq4_age50_screen_lead_time.RData")
cat(paste("\n"))
print(summary(lead_time[[1]]))
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))

#age 60-74
cat(paste("#### Your run of screen_freq = 4, age 60 - 74: Statistic summary \n" ))
cat(paste("\n"))
print(summary(Freq4_age60_lead_time[[1]]))
cat(paste("\n\n\n"))

cat(paste("#### First output of screen_freq = 4, age 60 - 74 from Intermediate_results folder: Statistic summary\n" ))
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq4_age60_screen_lead_time.RData")
cat(paste("\n"))
print(summary(lead_time[[1]]))
cat(paste("********************************************************\n"))
cat(paste("\n\n\n"))
sink()

setwd("~/Biometrical/C++_code")
dyn.unload("MSM_simulation_Windows.dll")
end_time <- Sys.time()
end_time - start_time



