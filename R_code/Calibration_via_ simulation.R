#author : Mai Ngoc Bui
#Calibration with Emperical data from Callender et al.
start_time <- Sys.time()
setwd("~/Biometrical/C++_code")
dyn.load("MSM_simulation_Windows.dll")#load C++ complied file


source("~/Biometrical/R_code/Simulation_functions.R")
#### MSM simulation ####
age_centered <- 68 #centering age

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
screen_freq  = 4 #frequency of screening.
#Gompertz distribution for all transitions in model B.
hazard_sub_distr <- c(2,2)



N <- 10^6 #population size.


#### BASELINE OF AGE = 55 ####
#age of all cohorts at the start of the study.
base_age = 55
#simulation for control arm
  data_arm2_age55 <- control_sim_Cpp_Rseed(par,  study_years, base_age, prop_year, left_trunc, id = 1:N, hazard_distr, max_no_seeds = 2)

#simulation for screen arm
  data_arm1_age55 <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (N+1):(2*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
save(data_arm1_age55, file = "data_arm1_simulated_age55.RData")
save(data_arm2_age55, file = "data_arm2_simulated_age55.RData")


#### BASELINE OF AGE = 40 ####
#age of all cohorts at the start of the study.
base_age = 40
#simulation for control arm

  data_arm2_age40 <- control_sim_Cpp_Rseed(par,  study_years, base_age, prop_year, left_trunc, id = 1:N, hazard_distr, max_no_seeds = 5)


#simulation for screen arm

  data_arm1_age40 <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id = (N+1):(2*N), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
save(data_arm1_age40, file = "data_arm1_simulated_age40.RData")
save(data_arm2_age40, file = "data_arm2_simulated_age40.RData")

#### BASELINE OF AGE = 55, SCREENING ONLY QUARTER OF ALL COHORTS EACH YEAR ####
#vary the screen_age to be (55,69), (56,69), (57,69) and (58,69) and simulate them separately.
N <- 10^6
base_age <- 55
data_arm1_age55_quarter_a <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age = c(55,69),  study_years, base_age, prop_year, left_trunc, id = N + ((0 * N/4 + 1): (1 * N/4)), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
data_arm1_age55_quarter_b <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age = c(56,69),  study_years, base_age, prop_year, left_trunc, id = N + ((1 * N/4 + 1): (2 * N/4)), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
data_arm1_age55_quarter_c <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age = c(57,69),  study_years, base_age, prop_year, left_trunc, id = N + ((2 * N/4 + 1): (3 * N/4)), hazard_distr, hazard_sub_distr, max_no_seeds = 2)
data_arm1_age55_quarter_d <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up,screen_freq, screen_age = c(58,69),  study_years, base_age, prop_year, left_trunc, id = N + ((3 * N/4 + 1): (4 * N/4)), hazard_distr, hazard_sub_distr, max_no_seeds = 2)

#combine these 4 datasets.
data_arm1_age55_quarter <- data_arm1_age55_quarter_a
data_arm1_age55_quarter$natural_history <- rbind(data_arm1_age55_quarter_a$natural_history,
                                                 data_arm1_age55_quarter_b$natural_history,
                                                 data_arm1_age55_quarter_c$natural_history,
                                                 data_arm1_age55_quarter_d$natural_history)
data_arm1_age55_quarter$life_time_record <- rbind(data_arm1_age55_quarter_a$life_time_record,
                                                  data_arm1_age55_quarter_b$life_time_record,
                                                  data_arm1_age55_quarter_c$life_time_record,
                                                  data_arm1_age55_quarter_d$life_time_record)
data_arm1_age55_quarter$study_record <- rbind(data_arm1_age55_quarter_a$study_record,
                                              data_arm1_age55_quarter_b$study_record,
                                              data_arm1_age55_quarter_c$study_record,
                                              data_arm1_age55_quarter_d$study_record)

setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
save(data_arm1_age55_quarter, file = "data_arm1_simulated_age55_quarter.RData")


end_time <- Sys.time()
end_time - start_time
#Time difference of 2.312708 mins

#### MIXED BASELINE OF AGES ####
#freqency of men in UK by ages (2017)
setwd("~/Biometrical/Supplement_data")
require(readxl)
male_2017 <- read_excel("ukmidyearestimates2017finalversion.xls",
                        sheet = "MYE2 - M", skip = 3)
male_2017 <- male_2017[!is.na(male_2017$Name),]
uk <- male_2017[male_2017$Name == "UNITED KINGDOM",]
age_uk <- data.frame(1:91,as.numeric(uk[5:95]))
names(age_uk) <- c("age","freq")

#distribution at baseline of age
base_age <- 55:69
age_vec <- as.factor(unlist(lapply(which(age_uk$age %in% base_age ), function(i) rep(age_uk$age[i],age_uk$freq[i]))))
prop_year <- as.vector(prop.table(table(age_vec)))

N <- 1000
mixed_data <- list()
mixed_data[[1]] <- screen_sim_Cpp_Rseed(par,sub_par,screen_misc, take_up = 0.9,screen_freq = 1, screen_age = c(55,69),  study_years = c(2000,2035), base_age, prop_year, left_trunc = 40, id = (N+1):(2*N), hazard_distr, hazard_sub_distr, max_no_seeds = 3)
mixed_data[[2]] <- control_sim_Cpp_Rseed(par,  study_years = c(2000,2035), base_age, prop_year, left_trunc = 40, id = 1:N, hazard_distr, max_no_seeds = 3)
names(mixed_data) <- c("screened_group","control_group")
setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
save(mixed_data, file = "mixed_simulated_data.RData")

#### SOJOURN TIME ####
N <- 10^7
base_age <- 40
system.time(
  data_arm2 <- control_sim_Cpp_Rseed(par,  study_years, base_age, prop_year, left_trunc, id = 1:N, hazard_distr, max_no_seeds = 2)
)
index_cancer <- which(data_arm2$natural_history$state == 3)
soj_time <- data.frame("time_state3" = data_arm2$natural_history$age[index_cancer], "time_state2"= data_arm2$natural_history$age[index_cancer - 1], "soj_time" = data_arm2$natural_history$age[index_cancer] - data_arm2$natural_history$age[index_cancer - 1])
setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
save(soj_time, file = "sojourn_time_age40.RData")



setwd("~/Biometrical/C++_code")
dyn.unload("MSM_simulation_Windows.dll")
end_time <- Sys.time()
end_time - start_time #Time difference of 1.234004 mins