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


#### FREQ = 2 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 2 #frequency of screening.
screen_age = c(55,69) #range of screening age.
lead_time <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)



  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(lead_time,file = "Freq2_screen_lead_time.RData")
  }
}


#### FREQ = 8 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 8 #frequency of screening.
screen_age = c(55,69) #range of screening age.
lead_time <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)



  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(lead_time,file = "Freq8_screen_lead_time.RData")
  }
}

#### FREQ = 16 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 16 #frequency of screening.
screen_age = c(55,69) #range of screening age.
lead_time <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)



  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(lead_time,file = "Freq16_screen_lead_time.RData")
  }
}

#### Age 55 - 69 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_freq  = 4 #frequency of screening.
screen_age = c(55,69) #range of screening age.
lead_time <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)



  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(lead_time,file = "Freq4_age55_screen_lead_time.RData")
  }
}


#### Age 50 - 64 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_age = c(50,64) #range of screening age.
screen_freq  = 4 #frequency of screening.
lead_time <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)



  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(lead_time,file = "Freq4_age50_screen_lead_time.RData")
  }
}


#### Age 60 - 74 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
screen_age = c(60,74) #range of screening age.
screen_freq  = 4 #frequency of screening.
lead_time <- list()
for (j in 1:500){
  data_arm1<- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (j * N+1):((j + 1)*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)



  screen_detected <- data_arm1$life_time_record[data_arm1$life_time_record$state == 2,]
  clinical_diag <- data_arm1$natural_history[(data_arm1$natural_history$id %in% screen_detected$id) & (data_arm1$natural_history$state == 3),]

  lead_time[[j]] <- clinical_diag$age - screen_detected$age[screen_detected$id %in% clinical_diag$id]


  if( (j %% 50) == 0){
    cat(paste("*** j = ",j,"***"))
    save(lead_time,file = "Freq4_age60_screen_lead_time.RData")
  }
}


######
end_time <- Sys.time()
end_time - start_time


