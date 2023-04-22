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
screen_age = c(55,69) #range of screening age.
#Gompertz distribution for all transitions in model B.
hazard_sub_distr <- c(2,2)



N <- 5*10^5 #population size.

data_arm2 <- control_sim_Cpp_Rseed(par,  study_years, base_age, prop_year, left_trunc, id = 1:N, hazard_distr, max_no_seeds = 2)
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
save(data_arm2, "data_arm2_age40.RData")

#### FREQ = 2 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 2 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

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

  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(output_rep,file = "Freq2_screen.RData")
  }
}


#### FREQ = 8 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 8 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
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

  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(output_rep,file = "Freq8_screen.RData")
  }
}

#### FREQ = 16 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 16 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
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


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(output_rep,file = "Freq16_screen.RData")
  }
}

#### Age 55 - 69 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 4 #frequency of screening.
screen_age = c(55,69) #range of screening age.
output_rep <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
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


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(output_rep,file = "Freq4_age55_screen.RData")
  }
}


#### Age 50 - 64 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_age = c(50,64) #range of screening age.
screen_freq  = 4 #frequency of screening.
output_rep <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
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


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(output_rep,file = "Freq4_age50_screen.RData")
  }
}


#### Age 60 - 74 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_age = c(60,74) #range of screening age.
screen_freq  = 4 #frequency of screening.
output_rep <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
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


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(output_rep,file = "Freq4_age60_screen.RData")
  }
}


######
end_time <- Sys.time()
end_time - start_time
#take up to 9 hr for N = 10^6
#take up to 5 hr for N = 5 * 10^5