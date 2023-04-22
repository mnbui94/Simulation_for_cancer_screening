#author : Mai Ngoc Bui
#function for parameter estimations.

#likelihood function for 4-states multi-state models on Windows
lik_grid_4States <- function(no_obs,par,data){
  grid_pw <- 1
  left_age <- 40 - age_centered
  min_del <- 1/365

  .Call("lik_grid_4States",as.double(par[1:3]), as.double(c(par[4:5],0)), as.double(par[6:7]), as.integer(data$state), as.double(data$centered_age), as.integer(length(no_obs) - 1), as.integer(no_obs), as.integer(data$screen),as.double(grid_pw), as.double(left_age), as.double(min_del))
}

#likelihood function for 4-states multi-state models on Linux
lik_grid_4States_parallel <- function(mul_no_obs,p,mul_data,no_cores){
  grid_pw <- 1
  left_age <- 40 - age_centered
  min_del <- 1/365

  sum(unlist(parallel::mclapply(1:no_cores, function(i).Call("lik_grid_4States",as.double(p[1:3]), as.double(c(p[4:5],0)), as.double(p[6:7]), as.integer(mul_data[[i]]$state), as.double(mul_data[[i]]$centered_age), as.integer(length(mul_no_obs[[i]]) - 1), as.integer(mul_no_obs[[i]]), as.integer(mul_data[[i]]$screen),as.double(grid_pw), as.double(left_age), as.double(min_del)),mc.cores = no_cores)))
}



#likelihood function for 3-states multi-state models
lik_fixed_3States <- function(no_obs,p,data){
  .Call("lik_fixed_3States",as.double(p[1:2]), as.double(c(p[3:4])), as.integer(data$state), as.double(data$age - age_centered), as.integer(length(no_obs) - 1), as.integer(no_obs))
}