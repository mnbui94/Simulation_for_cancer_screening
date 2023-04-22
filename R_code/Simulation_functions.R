#author : Mai Ngoc Bui
#function for simulation scheme.


#simulation for natural history function:
natural_sim_Cpp_Rseed <- function(p,  study_years, base_age, prop_year, left_trunc, id, hazard_distr){
  set.seed(id[1])
  random_seeds <- array(runif(6*length(id)),dim= c(6, length(id)))
  if(sum(hazard_distr == 0) != 0){
    set.seed(id[1])
    random_seeds[which(hazard_distr == 0),] <- rexp(N, exp(p$lambda[which(hazard_distr == 0)]))
  }

  set.seed(id[1])
  if(length(base_age) == 1){
    initial_years <- replicate(length(id),study_years[1] - base_age + left_trunc)
  } else {
    initial_years <- study_years[1] - replicate(length(id),sample(base_age, size = 1, prob = prop_year)) + left_trunc
  }
  t0 <- left_trunc - age_centered
  ans <- .Call("R_many_natural_history_R_seed", as.double(p$lambda), as.double(p$beta), as.double(t0), as.integer(initial_years), as.integer(hazard_distr), as.integer(length(id)), as.integer(id), as.double(random_seeds))
  #1st column = id; 2nd column = state; 3rd column = year; 4th column = time
  data.frame("id" = ans[,1], "state" = ans[,2], "age" = round(ans[,4] + age_centered, digits = 3),"year" = ans[,3], "time" = ans[,4])
}
#simulation for control group function:
control_sim_Cpp_Rseed <- function(p,  study_years, base_age, prop_year, left_trunc, id, hazard_distr,max_no_seeds){
  set.seed(id[1])
  random_seeds <- runif(6 * length(id) * max_no_seeds)
  if(sum(hazard_distr == 0) != 0){
    set.seed(id[1])
    random_seeds[seq(from = which(hazard_distr == 0), by = 6, length.out= max_no_seeds * length(id))] <- rexp(max_no_seeds * length(id), exp(p$lambda[which(hazard_distr == 0)]))
  }

  set.seed(id[1])
  if(length(base_age) == 1){
    initial_years <- replicate(length(id),study_years[1] - base_age + left_trunc)
  } else {
    initial_years <- study_years[1] - replicate(length(id),sample(base_age, size = 1, prob = prop_year)) + left_trunc
  }
  t0 <- left_trunc - age_centered

  ans <- .Call("R_many_simulation_control_R_seed", as.double(p$lambda), as.double(p$beta), as.double(t0), as.integer(initial_years),as.integer(study_years), as.integer(hazard_distr), as.integer(length(id)), as.integer(id), as.double(random_seeds), as.integer(max_no_seeds))


  #1st column = id; 2nd column = state; 3rd column = year; 4th column = time
  output <- list()
  output[[1]] <- data.frame("id" = ans[[1]][,1], "state" = ans[[1]][,2], "age" = round(ans[[1]][,4] + age_centered, digits = 3),"year" = ans[[1]][,3], "time" = ans[[1]][,4])
  output[[2]] <- data.frame("id" = ans[[2]][,1], "state" = ans[[2]][,2], "age" = round(ans[[2]][,4] + age_centered, digits = 3),"year" = ans[[2]][,3],"screen" = rep(0, nrow(ans[[2]])), "time" = ans[[2]][,4])
  output[[3]] <- data.frame("id" = ans[[3]][,1], "state" = ans[[3]][,2], "age" = round(ans[[3]][,4] + age_centered, digits = 3),"year" = ans[[3]][,3], "arm" = rep(2, nrow(ans[[3]])),"screen" = rep(0,nrow(ans[[3]])), "time" = ans[[3]][,4])
  names(output) <- c("natural_history", "life_time_record", "study_record")
  return(output)
}
#simulation for screened group function:
screen_sim_Cpp_Rseed <- function(p,sub_p,screen_misc, take_up,screen_freq, screen_age,  study_years, base_age, prop_year, left_trunc, id, hazard_distr, hazard_sub_distr, max_no_seeds){
  set.seed(id[1])
  random_seeds <- runif(6 * length(id) * max_no_seeds)
  random_seeds_sub <- runif(2 * length(id))
  if(sum(hazard_distr == 0) != 0){
    set.seed(id[1])
    random_seeds[seq(from = which(hazard_distr == 0)[1], by = 6, length.out= 1 *  max_no_seeds * length(id))] <- rexp(max_no_seeds * length(id), exp(p$lambda[which(hazard_distr == 0)]))
  }
  if(sum(hazard_sub_distr == 0) !=0){
    set.seed(id[1])
    random_seeds_sub[seq(from = which(hazard_sub_distr == 0)[1], by = 2, length.out=1*length(id))] <- rexp(1 * length(id), exp(sub_p$lambda[which(hazard_sub_distr == 0)]))
  }

  set.seed(id[1])
  if(length(base_age) == 1){
    initial_years <- replicate(length(id),study_years[1] - base_age + left_trunc)
  } else {
    initial_years <- study_years[1] - replicate(length(id),sample(base_age, size = 1, prob = prop_year)) + left_trunc
  }
  t0 <- left_trunc - age_centered
  ans <- .Call("R_many_simulation_screen_freq_R_seed", as.double(p$lambda), as.double(p$beta), as.double(sub_p$lambda),as.double(sub_p$beta),as.double(screen_misc), as.double(t0), as.integer(initial_years),as.integer(study_years), as.double(screen_age - age_centered), as.integer(screen_freq),as.integer(round((screen_age[2] - screen_age[1] + 1)/screen_freq)), as.double(take_up),  as.integer(hazard_distr), as.integer(hazard_sub_distr), as.integer(length(id)), as.integer(id), as.double(random_seeds),as.double(random_seeds_sub), as.integer(max_no_seeds))



  #1st column = id; 2nd column = state; 3rd column = year; 4th column = time
  output <- list()
  output[[1]] <- data.frame("id" = ans[[1]][,1], "state" = ans[[1]][,2], "age" = round(ans[[1]][,4] + age_centered, digits = 3),"year" = ans[[1]][,3], "time" = ans[[1]][,4])
  output[[2]] <- data.frame("id" =  ans[[2]][,1], "state" = ans[[2]][,2], "age" = round(ans[[2]][,4] + age_centered, digits = 3),"year" = ans[[2]][,3], "arm" = rep(1, nrow(ans[[2]])),"screen" = ans[[2]][,5], "time" = ans[[2]][,4])
  output[[3]] <- data.frame("id" = ans[[3]][,1], "state" = ans[[3]][,2], "age" = round(ans[[3]][,4] + age_centered, digits = 3),"year" = ans[[3]][,3], "arm" = rep(1, nrow(ans[[3]])), "screen" = ans[[3]][,5],"time" = ans[[3]][,4])
  names(output) <- c("natural_history", "life_time_record", "study_record")
  return(output)
}