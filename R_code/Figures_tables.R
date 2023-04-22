#author : Mai Ngoc Bui
#all figures and tables used in the manuscript
start_time <- Sys.time()
## PACKAGES ##
library(msm)
library(pbapply)
library(extraDistr)
library(readxl)
library(ggpubr)

#### LOAD DATA ####
#emperical data
setwd("~/Biometrical/Supplement_data")
emperical_arm2 <- read_excel("noMRI_cohorts.xlsx",sheet = "no-screening")
emperical_arm1 <- read_excel("noMRI_cohorts.xlsx",sheet = "age-based-screening")
total_emperical <- 341434

#routine data
suppressWarnings(all_deaths_age <- as.data.frame(read_excel("drtables17.xls",sheet = "Table 2", skip = 8))[-c(1:78,156:161),])

#simulated data
setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
load("data_arm1_simulated_age40.RData")
load("data_arm2_simulated_age40.RData")
load("data_arm1_simulated_age55.RData")
load("data_arm2_simulated_age55.RData")
load("data_arm1_simulated_age55_quarter.RData")
#### TAB. 1####
setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
load("mixed_simulated_data.RData")


setwd("~/Biometrical/Figures_Tables")
sink("table_1.txt")
id <- c(1,7,81,2)
for ( i in 1:length(id)){
  cat(paste("******* TABLE 1 - Individual", i, " *******\n"))
  print(mixed_data$control_group$natural_history[mixed_data$control_group$natural_history$id == id[i],])
}
sink()


sink("figure_3.txt")
print("******* FIGURE 3 - NATURAL HISTORY *******")
print(mixed_data$control_group$natural_history[mixed_data$control_group$natural_history$id == 81,])
print("******* FIGURE 3 - OBSERVED RECORD *******")
print(mixed_data$control_group$study_record[mixed_data$control_group$study_record$id == 81,])
cat(paste("\n\n"))
print("******* FIGURE 3 - NATURAL HISTORY *******")
print(mixed_data$control_group$natural_history[mixed_data$control_group$natural_history$id == 10,])
print("******* FIGURE 3 - OBSERVED RECORD *******")
print(mixed_data$control_group$study_record[mixed_data$control_group$study_record$id == 10,])
sink()



sink("figure_4.txt")
print("******* FIGURE 4 - NATURAL HISTORY *******")
mixed_data$screened_group$natural_history[mixed_data$screened_group$natural_history$id == 1236,]
print("******* FIGURE 4 - OBSERVED RECORD *******")
mixed_data$screened_group$study_record[mixed_data$screened_group$study_record$id == 1236,]
cat(paste("\n\n"))
print("******* FIGURE 4 - NATURAL HISTORY *******")
mixed_data$screened_group$natural_history[mixed_data$screened_group$natural_history$id == 1958,]
print("******* FIGURE 4 - OBSERVED RECORD *******")
mixed_data$screened_group$study_record[mixed_data$screened_group$study_record$id == 1958,]
sink()



#### FIG.5 & S2, TAB.4 ####
#figure 5#



total_simulated <- length(unique(data_arm1_age55$life_time_record$id))#sample size

#find the age at which cancer cases are found (round up the age).
who_has_cancer_arm1 <- data_arm1_age55$life_time_record[data_arm1_age55$life_time_record$state %in% c(2,3),]
cases_arm1 <- ceiling(who_has_cancer_arm1$age[!duplicated(who_has_cancer_arm1$id)])
who_has_cancer_arm2 <- data_arm2_age55$life_time_record[data_arm2_age55$life_time_record$state %in% c(2,3),]
cases_arm2 <- ceiling(who_has_cancer_arm2$age[!duplicated(who_has_cancer_arm2$id)])

#tranform to frequency table of age 55-89
freq_cases_arm1 <- as.data.frame(table(cases_arm1))[1:35,]
freq_cases_arm2 <- as.data.frame(table(cases_arm2))[1:35,]

#find the age at which men died from cancer (round up the age).
cancer_deaths_arm1 <- ceiling((data_arm1_age55$life_time_record[data_arm1_age55$life_time_record$state == 5,])$age)
cancer_deaths_arm2 <- ceiling((data_arm2_age55$life_time_record[data_arm2_age55$life_time_record$state == 5,])$age)


#tranform to frequency table of age 55-89
freq_cancer_deaths_arm1 <- as.data.frame(table(cancer_deaths_arm1))[1:35,]
freq_cancer_deaths_arm2 <- as.data.frame(table(cancer_deaths_arm2))[1:35,]


#find the age at which men died from cancer (round up the age).
other_deaths_arm1 <- ceiling((data_arm1_age55$life_time_record[data_arm1_age55$life_time_record$state == 4,])$age)

other_deaths_arm2 <- ceiling((data_arm2_age55$life_time_record[data_arm2_age55$life_time_record$state == 4,])$age)

#tranform to frequency table of age 55-89
freq_other_deaths_arm1 <- as.data.frame(table(other_deaths_arm1))[1:35,]
freq_other_deaths_arm2 <- as.data.frame(table(other_deaths_arm2))[1:35,]

#TAB 4 #
incidence_rate_emprical <- c(sum(emperical_arm1$pca_cases[emperical_arm1$age %in% 55:58]),
                             sum(emperical_arm1$pca_cases[emperical_arm1$age %in% 59:62]),
                             sum(emperical_arm1$pca_cases[emperical_arm1$age %in% 63:66]),
                             sum(emperical_arm1$pca_cases[emperical_arm1$age %in% 67:69])) /     total_emperical
incidence_rate_simulated <- c(freq_cases_arm1$Freq[seq(from = 1, by = 4, length.out = 4)])/total_simulated
setwd("~/Biometrical/Figures_Tables")
sink("table_4.txt")
print(data.frame("age_group" = c("55-58","59-62", "63-66","67-69"),incidence_rate_simulated,incidence_rate_emprical))
sink()


setwd("~/Biometrical/Figures_Tables")
pdf(file = "compare_simulated_emperical.pdf",
    width = 8, height = 5.5)
par(mar = c(4.5,4,2,1.5),oma = c(0,0.5,0.5,2))
layout(mat = matrix(c(1,2,3,
                      4,5,6,
                      7,7,7),
                    nrow = 3,
                    ncol = 3,byrow = TRUE),
       heights = c(3.2,3.2,0.5),
       widths = c(3,3,3))

ymax = max(cumsum(freq_cases_arm1$Freq)/total_simulated,
           c(cumsum(emperical_arm1$pca_cases/total_emperical)))
plot(55:89,cumsum(freq_cases_arm1$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative cancer incidence",ylim = c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm1$pca_cases/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_cancer_deaths_arm1$Freq)/total_simulated,
           c(cumsum(emperical_arm1$deaths_pca/total_emperical)))
plot(55:89,cumsum(freq_cancer_deaths_arm1$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (cancer)",ylim = c(0, ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm1$deaths_pca/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_other_deaths_arm1$Freq)/total_simulated,
           c(cumsum(emperical_arm1$deaths_othercauses/total_emperical)))
plot(55:89,cumsum(freq_other_deaths_arm1$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (other)",ylim = c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm1$deaths_othercauses/total_emperical)), col = "gray",lwd = 1.8)



ymax = max(cumsum(freq_cases_arm2$Freq)/total_simulated,
           c(cumsum(emperical_arm2$pca_cases/total_emperical)))
plot(55:89,cumsum(freq_cases_arm2$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative cancer incidence",ylim =c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm2$pca_cases/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_cancer_deaths_arm2$Freq)/total_simulated,
           c(cumsum(emperical_arm2$deaths_pca/total_emperical)))
plot(55:89,cumsum(freq_cancer_deaths_arm2$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (cancer)",ylim = c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm2$deaths_pca/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_other_deaths_arm2$Freq)/total_simulated,
           c(cumsum(emperical_arm2$deaths_othercauses/total_emperical)))
plot(55:89,cumsum(freq_other_deaths_arm2$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (other)",ylim =c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm2$deaths_othercauses/total_emperical)), col = "gray",lwd = 1.8)


par(mar= c(0,0,0,0))
plot(0, type="n", xlab="", ylab="",
     xlim = c(0,40),ylim= c(0,5),axes=FALSE,ann=FALSE)
legend(3.2,5,
       c(expression(paste(" Simulated data from MSM models       ")),expression(paste(" Empirical data from Callender et al. "))
       ), col = c("black", "gray"), lty = c(6,1),
       lwd = rep(2,2), merge = TRUE,cex = 1.2, horiz = TRUE,bty = "n")



mtext("Screened group",at = c(0.80,2), side = 4,
      col = "dimgray", outer = TRUE, cex = 0.9,padj = 0)
mtext("Control group",at = c(0.35,2), side = 4,
      col = "dimgray", outer = TRUE, cex = 0.9,padj = 0)
dev.off()


#FIG S2 #
#find the age at which cancer cases are found (round up the age).
who_has_cancer_arm1 <- data_arm1_age55_quarter$life_time_record[data_arm1_age55_quarter$life_time_record$state %in% c(2,3),]
cases_arm1 <- ceiling(who_has_cancer_arm1$age[!duplicated(who_has_cancer_arm1$id)])
#tranform to frequency table of age 55-89
freq_cases_arm1 <- as.data.frame(table(cases_arm1))[1:35,]

#find the age at which men died from cancer (round up the age).
cancer_deaths_arm1 <- ceiling((data_arm1_age55_quarter$life_time_record[data_arm1_age55_quarter$life_time_record$state == 5,])$age)
#tranform to frequency table of age 55-89
freq_cancer_deaths_arm1 <- as.data.frame(table(cancer_deaths_arm1))[1:35,]
#find the age at which men died from cancer (round up the age).
other_deaths_arm1 <- ceiling((data_arm1_age55_quarter$life_time_record[data_arm1_age55_quarter$life_time_record$state == 4,])$age)
#tranform to frequency table of age 55-89
freq_other_deaths_arm1 <- as.data.frame(table(other_deaths_arm1))[1:35,]
setwd("~/Biometrical/Figures_Tables")
pdf(file = "compare_simulated_emperical_updated.pdf",
    width = 8, height = 5.5)
par(mar = c(4.5,4,2,1.5),oma = c(0,0.5,0.5,2))
layout(mat = matrix(c(1,2,3,
                      4,5,6,
                      7,7,7),
                    nrow = 3,
                    ncol = 3,byrow = TRUE),
       heights = c(3.2,3.2,0.5),
       widths = c(3,3,3))

ymax = max(cumsum(freq_cases_arm1$Freq)/total_simulated,
           c(cumsum(emperical_arm1$pca_cases/total_emperical)))
plot(55:89,cumsum(freq_cases_arm1$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative cancer incidence",ylim = c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm1$pca_cases/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_cancer_deaths_arm1$Freq)/total_simulated,
           c(cumsum(emperical_arm1$deaths_pca/total_emperical)))
plot(55:89,cumsum(freq_cancer_deaths_arm1$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (cancer)",ylim = c(0, ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm1$deaths_pca/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_other_deaths_arm1$Freq)/total_simulated,
           c(cumsum(emperical_arm1$deaths_othercauses/total_emperical)))
plot(55:89,cumsum(freq_other_deaths_arm1$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (other)",ylim = c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm1$deaths_othercauses/total_emperical)), col = "gray",lwd = 1.8)



ymax = max(cumsum(freq_cases_arm2$Freq)/total_simulated,
           c(cumsum(emperical_arm2$pca_cases/total_emperical)))
plot(55:89,cumsum(freq_cases_arm2$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative cancer incidence",ylim =c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm2$pca_cases/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_cancer_deaths_arm2$Freq)/total_simulated,
           c(cumsum(emperical_arm2$deaths_pca/total_emperical)))
plot(55:89,cumsum(freq_cancer_deaths_arm2$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (cancer)",ylim = c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm2$deaths_pca/total_emperical)), col = "gray",lwd = 1.8)

ymax = max(cumsum(freq_other_deaths_arm2$Freq)/total_simulated,
           c(cumsum(emperical_arm2$deaths_othercauses/total_emperical)))
plot(55:89,cumsum(freq_other_deaths_arm2$Freq)/total_simulated, type = "l",xlab = "Age",ylab= "Cumulative mortality (other)",ylim =c(0,ymax), lty = 6,lwd = 1.8, frame.plot = FALSE)
lines(55:89, c(cumsum(emperical_arm2$deaths_othercauses/total_emperical)), col = "gray",lwd = 1.8)


par(mar= c(0,0,0,0))
plot(0, type="n", xlab="", ylab="",
     xlim = c(0,40),ylim= c(0,5),axes=FALSE,ann=FALSE)
legend(3.2,5,
       c(expression(paste(" Simulated data from MSM models       ")),expression(paste(" Empirical data from Callender et al. "))
       ), col = c("black", "gray"), lty = c(6,1),
       lwd = rep(2,2), merge = TRUE,cex = 1.2, horiz = TRUE,bty = "n")



mtext("Screened group",at = c(0.80,2), side = 4,
      col = "dimgray", outer = TRUE, cex = 0.9,padj = 0)
mtext("Control group",at = c(0.35,2), side = 4,
      col = "dimgray", outer = TRUE, cex = 0.9,padj = 0)
dev.off()

#### FIG.6####


# routine data
all_deaths_age <- all_deaths_age[-seq(from = 1, to = nrow(all_deaths_age), by = 7),1:2]
all_deaths_age <- all_deaths_age[!is.na(all_deaths_age$`All ages`),]
names(all_deaths_age) <- c("age", "routine_freq")
# compute freq per 100,000
all_deaths_age$routine_freq <- 10^5 * all_deaths_age$routine_freq/sum(all_deaths_age$routine_freq)


#frequency of age at all-causes deaths (aged between 55-109)
age_deaths_arm1 <- table(factor(ceiling((data_arm1_age40$natural_history[data_arm1_age40$natural_history$state %in% c(4,5),])$age), levels = 55:109))
age_deaths_arm2 <- table(factor(ceiling((data_arm2_age40$life_time_record[data_arm2_age40$life_time_record$state %in% c(4,5),])$age), levels = 55:109))



#merge routine data and simulated data, where freq w.r.t population size = 10^5
sum_deaths_arm1 <- data.frame(all_deaths_age, "simulation_freq" = 10^5 * as.numeric(age_deaths_arm1) / sum(as.numeric(age_deaths_arm1)))
sum_deaths_arm2 <- data.frame(all_deaths_age, "simulation_freq" = 10^5 * as.numeric(age_deaths_arm2) / sum(as.numeric(age_deaths_arm2)))


age_control <- c(sapply(sum_deaths_arm2$age, function(i) rep(i, 2)))
name_control <- rep(c("Routine data", "Simulation") ,nrow(sum_deaths_arm2))
a_deaths_control <- rbind(sum_deaths_arm2$routine_freq,sum_deaths_arm2$simulation_freq)
data_control <- data.frame(age_control = as.numeric(age_control),name_control, "a_deaths_control" = c(a_deaths_control))


age_screen <- c(sapply(sum_deaths_arm1$age, function(i) rep(i, 2)))
name_screen <- rep(c("Routine data", "Simulation") ,nrow(sum_deaths_arm1))
a_deaths_screen <- rbind(sum_deaths_arm1$routine_freq,sum_deaths_arm1$simulation_freq)
data_screen <- data.frame(age_screen = as.numeric(age_screen),name_screen, "a_deaths_screen" = c(a_deaths_screen))




setwd("~/Biometrical/Figures_Tables")
plot_control <- ggplot(data_control, aes(fill=name_control,  show.legend = FALSE, y=a_deaths_control, x=age_control)) +
  geom_bar(position="dodge", stat="identity")+
  ggtitle("All deaths by ages (Control group)") +
  xlab("Age") + ylab("Subjects per 100,000")+
  labs(fill ="")  +scale_fill_manual(values=c( "#000000", "#A9A9A9"))
plot_screen <- ggplot(data_screen, aes(fill=name_screen, y=a_deaths_screen, x=age_screen)) +
  geom_bar(position="dodge", stat="identity")+
  ggtitle("All deaths by ages (Screened group)")+
  xlab("Age") + ylab("Subjects per 100,000")+
  labs(fill ="") +scale_fill_manual(values=c( "#000000", "#A9A9A9"))

ggarrange(plot_control, plot_screen,
          labels = c("(a)", "(b)"),
          ncol = 2, nrow = 1, common.legend = T,legend="bottom")


ggsave("all_deaths_by_ages.pdf", width = 8, height = 3.7)



#### TAB.S2####


age_group <- c("55-59","60-64","65-69","70-74","75-79","80-84","85-89")
#control arm#
prostate_cancer_cases <-  round(10^5 * rowSums(matrix(emperical_arm2$pca_cases,ncol = 5, byrow = TRUE))/total_emperical)
deaths_cancer <-  round(10^5 * rowSums(matrix(emperical_arm2$deaths_pca,ncol = 5, byrow = TRUE))/total_emperical)
deaths_other_causes <-  round(10^5 * rowSums(matrix(emperical_arm2$deaths_othercauses,ncol = 5, byrow = TRUE))/total_emperical)
control_s1 <- data.frame(age_group, prostate_cancer_cases, deaths_cancer, deaths_other_causes)
#screen arm#
prostate_cancer_cases <-  round(10^5 * rowSums(matrix(emperical_arm1$pca_cases,ncol = 5, byrow = TRUE))/total_emperical)
deaths_cancer <-  round(10^5 * rowSums(matrix(emperical_arm1$deaths_pca,ncol = 5, byrow = TRUE))/total_emperical)
deaths_other_causes <-  round(10^5 * rowSums(matrix(emperical_arm1$deaths_othercauses,ncol = 5, byrow = TRUE))/total_emperical)
screen_s1 <- data.frame(age_group, prostate_cancer_cases, deaths_cancer, deaths_other_causes)

setwd("~/Biometrical/Figures_Tables")
sink("table_s2.txt")
print("******* TABLE S2 - control *******")
print(control_s1)
print("******* TABLE S2 - screen *******")
print(screen_s1)
sink()


#### FIG. S1 #####

setwd("~/Biometrical/PLCO_data")
org_data <- read.csv("210125 prostate long data for msm -May 2015.csv")
full_data <- data.frame("id" = org_data$plco_id, "state" =  org_data$state, "age" =  org_data$age,
                        "screen" = org_data$no_screen, "arm" = org_data$arm,
                        "centered_age" = org_data$age - 68)
full_data$screen[full_data$arm == 2] <- FALSE
full_data$screen[full_data$arm == 1] <- is.na(full_data[(full_data$arm == 1),]$screen)


l_data<- split(full_data, full_data$id)

entry_age_arm1 <- round(unlist(lapply(l_data[unlist(lapply(l_data, function(i) i$arm[1] == 1))], function(x) x$age[1])), digits = 0)
entry_age_arm2 <- round(unlist(lapply(l_data[unlist(lapply(l_data, function(i) i$arm[1] == 2))], function(x) x$age[1])), digits = 0)



exit_age_arm1 <- round(unlist(lapply(l_data[unlist(lapply(l_data, function(i) i$arm[1] == 1))], function(x) x$age[nrow(x)])), digits = 0)
exit_age_arm2 <- round(unlist(lapply(l_data[unlist(lapply(l_data, function(i) i$arm[1] == 2))], function(x) x$age[nrow(x)])), digits = 0)

setwd("~/Biometrical/Figures_Tables")
pdf(file = "PLCO_age_description.pdf",
    width = 10, height = 7)
par(mfrow = c(2,2),oma = c(0,0,2,0))
barplot(table(entry_age_arm1), xlab = "Age",ylab = "Frequency", main = "Age at the entry (Screened group)")
barplot(table(entry_age_arm2), xlab = "Age",ylab = "Frequency", main = "Age at the entry (Control group)")
barplot(table(exit_age_arm1), xlab = "Age",ylab = "Frequency", main = "Age at the exit (Screened group)")
barplot(table(exit_age_arm2), xlab = "Age",ylab = "Frequency", main = "Age at the exit (Control group)")
mtext("Summary of ages at the entry and the exit from the PLCO data", col = "dimgray", outer = TRUE, cex = 1.3)
dev.off()





#### TAB. S1 ####
setwd("~/Biometrical/PLCO_data")
load("ModelA_4States_PLCOdata.RData")
load("ModelA_3States_PLCOdata.RData")
load("ModelB_3States_PLCOdata.RData")

setwd("~/Biometrical/Figures_Tables")
sink("table_s1.txt")
cat(paste("***** total cohorts in screened group :\n"))
length(unique(full_data$id[full_data$arm == 1]))
cat(paste("\n\n"))

cat(paste("***** total cohorts in control group :\n"))
length(unique(full_data$id[full_data$arm == 2]))
cat(paste("\n\n\n"))

cat(paste("***** state table of screened group :\n"))
print(statetable.msm(state, id, full_data[full_data$arm == 1,]))
cat(paste("\n*** From state 2 ***\n"))
print(statetable.msm(state, id, data_state2[data_state2$arm == 1,]))
cat(paste("\n*** From state 3 ***\n"))
print(statetable.msm(state, id, data_state3[data_state3$arm == 1,]))
cat(paste("\n\n\n"))
cat(paste("***** state table of control group :\n"))
print(statetable.msm(state, id, full_data[full_data$arm == 2,]))
cat(paste("\n*** From state 3 ***\n"))
print(statetable.msm(state, id, data_state3[data_state3$arm == 2,]))
sink()


#### FIG 7 & 8 ####
setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq2_screen.RData")
Freq2 <- output_rep
load("Freq8_screen.RData")
Freq8 <- output_rep
load("Freq16_screen.RData")
Freq16 <- output_rep


load("Freq4_age55_screen.RData")
Freq4 <- output_rep
load("Freq4_age50_screen.RData")
Freq4_age50 <- output_rep
load("Freq4_age60_screen.RData")
Freq4_age60 <- output_rep



rm(output_rep)
screen_type <- c("Freq2","Freq4","Freq8","Freq16","Freq4_age50","Freq4_age60")

MR <- list()
for (i in 1:6){
  MR[[i]] <- unlist(lapply(
    eval(parse (text = screen_type[i])), function(x) x$cancer_deaths_reduction_life_time$percentages))
}

OD <- list()
for (i in 1:6){
  OD[[i]] <- unlist(lapply(
    eval(parse (text = screen_type[i])), function(x) x$excess_incidence_life_time$percentages))
}
names(MR) = screen_type
names(OD) = screen_type


setwd("~/Biometrical/Figures_Tables")
pdf(file = "screen_evaluation_compare_freq_500.pdf",
    width = 13, height = 5.5)
par(mfrow = c(1,2))
x <- data.frame(OD$Freq16[1:200], OD$Freq8[1:200], OD$Freq4[1:200], OD$Freq2[1:200])
names(x) <- c("Testing once","Every 8 years","Every 4 years", "Every 2 years")
boxplot(x, main = "Overdiagnosis", ylab = "Percentages")


x <- data.frame(MR$Freq16[1:200], MR$Freq8[1:200], MR$Freq4[1:200], MR$Freq2[1:200])
names(x) <- c("Testing once","Every 8 years","Every 4 years", "Every 2 years")
boxplot(x, main = "Cancer-specific mortality reduction", ylab = "Percentages")

dev.off()


setwd("~/Biometrical/Figures_Tables")
pdf(file = "screen_evaluation_screening_age_500.pdf",
    width = 13, height = 5.5)
par(mfrow = c(1,2))

x <- data.frame(OD$Freq4_age50[1:200], OD$Freq4[1:200], OD$Freq4_age60[1:200])
names(x) <- c("50-64","55-69","60-74")
boxplot(x, main = "Overdiagnosis", ylab = "Percentages",xlab = "Screening age group")


x <- data.frame(MR$Freq4_age50[1:200], MR$Freq4[1:200], MR$Freq4_age60[1:200])
names(x) <- c("50-64","55-69","60-74")
boxplot(x, main = "Cancer-specific mortality reduction", ylab = "Percentages",xlab = "Screening age group")

dev.off()

####TAB.S3 & S4 ####


overdiagnosis_mean <- c(signif(mean(OD$Freq16), digits = 3),
                           signif(mean(OD$Freq8), digits = 3),
                           signif(mean(OD$Freq4), digits = 3),
                           signif(mean(OD$Freq2), digits = 3))
overdiagnosis_CI <- rbind(signif(t.test(OD$Freq16)$conf.int, digits = 3),
                          signif(t.test(OD$Freq8)$conf.int, digits = 3),
                          signif(t.test(OD$Freq4)$conf.int, digits = 3),
                          signif(t.test(OD$Freq2)$conf.int, digits = 3))
OD_CI <- cbind(overdiagnosis_mean,overdiagnosis_CI)
colnames(OD_CI) <- c("mean","lower bound (95% CI)","upper bound (95% CI)")


cancer_deaths_reduction_mean <- c(signif(mean(MR$Freq16), digits = 3),
                        signif(mean(MR$Freq8), digits = 3),
                        signif(mean(MR$Freq4), digits = 3),
                        signif(mean(MR$Freq2), digits = 3))
cancer_deaths_reduction_CI <- rbind(signif(t.test(MR$Freq16)$conf.int, digits = 3),
                          signif(t.test(MR$Freq8)$conf.int, digits = 3),
                          signif(t.test(MR$Freq4)$conf.int, digits = 3),
                          signif(t.test(MR$Freq2)$conf.int, digits = 3))
MR_CI <- cbind(cancer_deaths_reduction_mean,cancer_deaths_reduction_CI)
colnames(MR_CI) <- c("mean","lower bound (95% CI)","upper bound (95% CI)")

setwd("~/Biometrical/Figures_Tables")
sink("table_s3.txt")
cat(paste("*************** Overdiagnosis *************** \n"))
print(OD_CI)
cat(paste("\n\n\n"))
cat(paste("****** Prostate cancer-specific mortality reduction ****** \n"))
print(MR_CI)
cat(paste("\n\n\n"))
sink()




overdiagnosis_mean <- c(signif(mean(OD$Freq4_age50), digits = 3),
                        signif(mean(OD$Freq4), digits = 3),
                        signif(mean(OD$Freq4_age60), digits = 3))
overdiagnosis_CI <- rbind(signif(t.test(OD$Freq4_age50)$conf.int, digits = 3),
                          signif(t.test(OD$Freq4)$conf.int, digits = 3),
                          signif(t.test(OD$Freq4_age60)$conf.int, digits = 3))
OD_CI <- cbind(overdiagnosis_mean,overdiagnosis_CI)
colnames(OD_CI) <- c("mean","lower bound (95% CI)","upper bound (95% CI)")


cancer_deaths_reduction_mean <- c(signif(mean(MR$Freq4_age50), digits = 3),
                                  signif(mean(MR$Freq4), digits = 3),
                                  signif(mean(MR$Freq4_age60), digits = 3))
cancer_deaths_reduction_CI <- rbind(signif(t.test(MR$Freq4_age50)$conf.int, digits = 3),
                                    signif(t.test(MR$Freq4)$conf.int, digits = 3),
                                    signif(t.test(MR$Freq4_age60)$conf.int, digits = 3))
MR_CI <- cbind(cancer_deaths_reduction_mean,cancer_deaths_reduction_CI)
colnames(MR_CI) <- c("mean","lower bound (95% CI)","upper bound (95% CI)")

setwd("~/Biometrical/Figures_Tables")
sink("table_s4.txt")
cat(paste("*************** Overdiagnosis *************** \n"))
print(OD_CI)
cat(paste("\n\n\n"))
cat(paste("****** Prostate cancer-specific mortality reduction ****** \n"))
print(MR_CI)
cat(paste("\n\n\n"))
sink()
#### TAB 6 ####

screen_frequency <- c("No screening", "Testing once", "Every 8 years", "Every 4 years", "Every 2 years", "Every 4 years", "Every 4 years")
testing_age <- c("N/A","55-69","55-69","55-69","55-69","50-64","60-74")

cancer_cases <- c(Freq2[[1]]$excess_incidence_life_time$total_cancer_arm2,
                  round(mean(unlist(lapply(Freq16, function(x) x$excess_incidence_life_time$total_cancer_arm1)))),
                  round(mean(unlist(lapply(Freq8, function(x) x$excess_incidence_life_time$total_cancer_arm1)))),
                  round(mean(unlist(lapply(Freq4, function(x) x$excess_incidence_life_time$total_cancer_arm1)))),
                  round(mean(unlist(lapply(Freq2, function(x) x$excess_incidence_life_time$total_cancer_arm1)))),
                  round(mean(unlist(lapply(Freq4_age50, function(x) x$excess_incidence_life_time$total_cancer_arm1)))),
                  round(mean(unlist(lapply(Freq4_age60, function(x) x$excess_incidence_life_time$total_cancer_arm1)))))

detected_cases <- c(0,
                  round(mean(unlist(lapply(Freq16, function(x) x$excess_incidence_life_time$total_screen_detected)))),
                  round(mean(unlist(lapply(Freq8, function(x) x$excess_incidence_life_time$total_screen_detected)))),
                  round(mean(unlist(lapply(Freq4, function(x) x$excess_incidence_life_time$total_screen_detected)))),
                  round(mean(unlist(lapply(Freq2, function(x) x$excess_incidence_life_time$total_screen_detected)))),
                  round(mean(unlist(lapply(Freq4_age50, function(x) x$excess_incidence_life_time$total_screen_detected)))),
                  round(mean(unlist(lapply(Freq4_age60, function(x) x$excess_incidence_life_time$total_screen_detected)))))
clinical_cases <- cancer_cases - detected_cases
interval_cancers <- c(0,
                      round(mean(unlist(lapply(Freq16, function(x) x$total_interval_cancers)))),
                      round(mean(unlist(lapply(Freq8, function(x) x$total_interval_cancers)))),
                      round(mean(unlist(lapply(Freq4, function(x) x$total_interval_cancers)))),
                      round(mean(unlist(lapply(Freq2, function(x) x$total_interval_cancers)))),
                      round(mean(unlist(lapply(Freq4_age50, function(x) x$total_interval_cancers)))),
                      round(mean(unlist(lapply(Freq4_age60, function(x) x$total_interval_cancers)))))
overdiagnosis <- c(0,
                   round(mean(unlist(lapply(Freq16, function(x) x$excess_incidence_life_time$excess_cancer_cases)))),
                   round(mean(unlist(lapply(Freq8, function(x) x$excess_incidence_life_time$excess_cancer_cases)))),
                   round(mean(unlist(lapply(Freq4, function(x) x$excess_incidence_life_time$excess_cancer_cases)))),
                   round(mean(unlist(lapply(Freq2, function(x) x$excess_incidence_life_time$excess_cancer_cases)))),
                   round(mean(unlist(lapply(Freq4_age50, function(x) x$excess_incidence_life_time$excess_cancer_cases)))),
                   round(mean(unlist(lapply(Freq4_age60, function(x) x$excess_incidence_life_time$excess_cancer_cases)))))
deaths_from_cancer <- c(Freq2[[1]]$cancer_deaths_reduction_life_time$cancer_deaths_control,
                        round(mean(unlist(lapply(Freq16, function(x) x$cancer_deaths_reduction_life_time$cancer_deaths_screen)))),
                        round(mean(unlist(lapply(Freq8, function(x) x$cancer_deaths_reduction_life_time$cancer_deaths_screen)))),
                        round(mean(unlist(lapply(Freq4, function(x) x$cancer_deaths_reduction_life_time$cancer_deaths_screen)))),
                        round(mean(unlist(lapply(Freq2, function(x) x$cancer_deaths_reduction_life_time$cancer_deaths_screen)))),
                        round(mean(unlist(lapply(Freq4_age50, function(x) x$cancer_deaths_reduction_life_time$cancer_deaths_screen)))),
                        round(mean(unlist(lapply(Freq4_age60, function(x) x$cancer_deaths_reduction_life_time$cancer_deaths_screen)))))

tab6 <- data.frame(screen_frequency, testing_age, cancer_cases, detected_cases, clinical_cases, interval_cancers, overdiagnosis, deaths_from_cancer)
tab6[tab6==0] <- "N/A"

setwd("~/Biometrical/Figures_Tables")
sink("table_6.txt")
cat(paste("************************************************* TABLE 6 *************************************************\n\n"))
print(tab6)
sink()


#### TAB 5 ####
setwd("~/Biometrical/Intermediate_results/simulated_cohorts")
load("sojourn_time_age40.RData")

age_groups <- c("50-54", "55-59", "60-64")
mean_sojourn_time <- signif(c(mean(soj_time$soj_time[!(soj_time$time_state2 < 50) & (soj_time$time_state2 < 55)]),
                       mean(soj_time$soj_time[!(soj_time$time_state2 < 55) & (soj_time$time_state2 < 60)]),
                       mean(soj_time$soj_time[!(soj_time$time_state2 < 60) & (soj_time$time_state2 < 65)])), digits= 2)

median_sojourn_time <- signif(c(median(soj_time$soj_time[!(soj_time$time_state2 < 50) & (soj_time$time_state2 < 55)]),
                              median(soj_time$soj_time[!(soj_time$time_state2 < 55) & (soj_time$time_state2 < 60)]),
                              median(soj_time$soj_time[!(soj_time$time_state2 < 60) & (soj_time$time_state2 < 65)])), digits= 2)
lower_confidence_interval <- signif(c(t.test(soj_time$soj_time[!(soj_time$time_state2 < 50) & (soj_time$time_state2 < 55)])$conf.int[1],
                                      t.test(soj_time$soj_time[!(soj_time$time_state2 < 55) & (soj_time$time_state2 < 60)])$conf.int[1],
                                      t.test(soj_time$soj_time[!(soj_time$time_state2 < 60) & (soj_time$time_state2 < 65)])$conf.int[1]), digits= 3)

upper_confidence_interval <- signif(c(t.test(soj_time$soj_time[!(soj_time$time_state2 < 50) & (soj_time$time_state2 < 55)])$conf.int[2],
                                      t.test(soj_time$soj_time[!(soj_time$time_state2 < 55) & (soj_time$time_state2 < 60)])$conf.int[2],
                                      t.test(soj_time$soj_time[!(soj_time$time_state2 < 60) & (soj_time$time_state2 < 65)])$conf.int[2]), digits= 3)

tab5 <- data.frame(age_groups, median_sojourn_time, mean_sojourn_time, lower_confidence_interval, upper_confidence_interval)

setwd("~/Biometrical/Figures_Tables")
sink("table_5.txt")
cat(paste("************************************************* TABLE 5 *************************************************\n\n"))
print(tab5)
sink()

#### TAB.2 ####

setwd("~/Biometrical/Estimation_results")
load("est1_ModelA_4States.RData")
load("est1_ModelA_3States.RData")
load("est1_ModelB_3States.RData")


ModelA <- c("lambda_{12}", "lambda_{14} = lambda_{24}","lambda_{23}","lambda_{34}", "lambda_{35}","beta_{12}","beta_{14} = beta_{24}","beta_{34}","beta_{35}","nu","gamma")
ModelB <- c("lambda^*_{24}", "lambda^*_{25}","beta^*_{24}","beta^*_{25}")

modelA_par_3states    <- round(est1_ModelA_3States$par, digits = 2)
modelA_par_3states_SE <- sqrt(diag(solve(est1_ModelA_3States$hessian)))

modelB_par_3states    <- round(est1_ModelB_3States$par, digits = 2)
modelB_par_3states_SE <- sqrt(diag(solve(est1_ModelB_3States$hessian)))

modelA_par_4states    <- round(est1_ModelA_4states$par, digits = 2)
modelA_par_4states_SE <- sqrt(diag(solve(est1_ModelA_4states$hessian)))



Estimated_values <- c(modelA_par_4states[1:3],modelA_par_3states[1:2],modelA_par_4states[4:5],modelA_par_3states[3:4],modelA_par_4states[6:7] )
SE <- round(c(modelA_par_4states_SE[1:3],modelA_par_3states_SE[1:2],modelA_par_4states_SE[4:5],modelA_par_3states_SE[3:4],modelA_par_4states_SE[6:7] ), digits = 2)
lower_95CI <- round(Estimated_values - 1.96 *SE, digits = 2)
upper_95CI <- round(Estimated_values + 1.96 *SE, digits = 2)

setwd("~/Biometrical/Figures_Tables")
sink("table_2.txt")
cat(paste("********* Estimated parameters of Model A using PLCO data *********\n\n"))
print(data.frame(ModelA, Estimated_values, SE, lower_95CI, upper_95CI))
cat(paste("\n\n\n"))
sink()




Estimated_values <- modelB_par_3states
SE <- round(modelB_par_3states_SE, digits = 2)
lower_95CI <- round(Estimated_values - 1.96 *SE, digits = 2)
upper_95CI <- round(Estimated_values + 1.96 *SE, digits = 2)

sink("table_2.txt", append = TRUE)
cat(paste("********* Estimated parameters of Model B using PLCO data *********\n\n"))
print(data.frame(ModelB, Estimated_values, SE, lower_95CI, upper_95CI))
sink()

####TAB.3####
setwd("~/Biometrical/Figures_Tables")
sink("table_3.txt")
cat(paste("********* Estimated parameters of Model A when callibrating against emperical data *********\n\n"))
print(data.frame("lambda" = c(-4.95,-4.4,-2.5,-0.07,-2.13),
                  "beta" = c(0.02, 0.12, 0.00, 0.01, 0.17)))
cat(paste("\n\n\n"))
cat(paste("********* Estimated parameters of Model B when callibrating against emperical data *********\n\n"))
print(data.frame("lambda" = c(-2.4,-5.1),
                      "beta" = c(0.06, 0.09)) )
sink()


#### FIG S3 ####

setwd("~/Biometrical/Intermediate_results/screen_evaluation_500k")
load("Freq2_screen_lead_time.RData")
freq2_mean_lead_time <- unlist(lapply(lead_time, mean))

load("Freq8_screen_lead_time.RData")
freq8_mean_lead_time <- unlist(lapply(lead_time, mean))

load("Freq16_screen_lead_time.RData")
freq16_mean_lead_time <- unlist(lapply(lead_time, mean))

load("Freq4_age55_screen_lead_time.RData")
freq4_age55_mean_lead_time <- unlist(lapply(lead_time, mean))

load("Freq4_age50_screen_lead_time.RData")
freq4_age50_mean_lead_time <- unlist(lapply(lead_time, mean))

load("Freq4_age60_screen_lead_time.RData")
freq4_age60_mean_lead_time <- unlist(lapply(lead_time, mean))

setwd("~/Biometrical/Figures_Tables")
pdf(file = "screen_evaluation_compare_lead_time_500.pdf",
    width = 13, height = 6)
par(mfrow = c(1,1))
x <- data.frame(freq16_mean_lead_time, freq8_mean_lead_time, freq4_age55_mean_lead_time, freq2_mean_lead_time, freq4_age50_mean_lead_time, freq4_age60_mean_lead_time)
names(x) <- c("Testing once (55-69)","Every 8 years (55-69)","Every 4 years (55-69)", "Every 2 years (55-69)", "Every 4 years (50-64)", "Every 4 years (60-74)")
boxplot(x, main = "Mean lead time", ylab = "Years")



dev.off()



end_time <- Sys.time()
end_time - start_time #Time difference of 22.00088 secs