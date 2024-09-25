Supplementary information / reproducible research files for the manuscript 
Title: "Non-Homogeneous multi-state Markov models: A simulation scheme for evaluating cancer screening strategies"

Authors: Mai Ngoc Bui; Ardo Van Den Hout; Rikesh Bhatt & Nora Pashayan
Code was written by Mai Ngoc Bui
In case of questions or comments please contact mai.bui@buv.edu.vn

The code in the R file "PLCO_data_estimation_Linux.R" was written/evaluated in R with the following software versions: 

	R version 4.1.2 (2021-11-01)
	Platform: x86_64-redhat-linux-gnu (64-bit)
	Running under: Rocky Linux 8.5 (Green Obsidian)

	Matrix products: default
	BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.12.so

	locale:
 	[1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 	[3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 	[5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
	[7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 	[9] LC_ADDRESS=C               LC_TELEPHONE=C
	[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

	attached base packages:
	[1] parallel  stats     graphics  grDevices utils     datasets  methods
	[8] base

	other attached packages:
	[1] data.table_1.14.2 ucminf_1.1-4      readxl_1.3.1      extraDistr_1.9.1
	[5] pbapply_1.5-0     msm_1.6.9

	loaded via a namespace (and not attached):
 	[1] compiler_4.1.2   Matrix_1.5-1     expm_0.999-6     survival_3.2-13
 	[5] Rcpp_1.0.7       cellranger_1.1.0 mvtnorm_1.1-3    splines_4.1.2
 	[9] grid_4.1.2       lattice_0.20-45




And the rest of the R code was written/evaluated in R with the following software versions:
	R version 4.2.2 (2022-10-31 ucrt)
	Platform: x86_64-w64-mingw32/x64 (64-bit)
	Running under: Windows 10 x64 (build 22621)

	Matrix products: default

	locale:
	[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8   
	[3] LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
	[5] LC_TIME=English_United Kingdom.utf8    

	attached base packages:
	[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

	other attached packages:
	[1] ucminf_1.1-4.1   ggpubr_0.5.0     ggplot2_3.4.0    readxl_1.4.1    
	[5] extraDistr_1.9.1 pbapply_1.6-0    msm_1.7         
	
	loaded via a namespace (and not attached):
 	[1] Rcpp_1.0.9       cellranger_1.1.0 pillar_1.8.1     compiler_4.2.2  
 	[5] tools_4.2.2      lifecycle_1.0.3  tibble_3.1.8     gtable_0.3.1    
 	[9] lattice_0.20-45  pkgconfig_2.0.3  rlang_1.0.6      Matrix_1.5-1    
	[13] cli_3.4.1        mvtnorm_1.1-3    expm_0.999-6     withr_2.5.0     
	[17] dplyr_1.0.10     generics_0.1.3   vctrs_0.5.1      grid_4.2.2      
	[21] tidyselect_1.2.0 glue_1.6.2       R6_2.5.1         rstatix_0.7.1   
	[25] fansi_1.0.3      survival_3.4-0   carData_3.0-5    car_3.1-1       
	[29] purrr_1.0.0      tidyr_1.2.1      magrittr_2.0.3   backports_1.4.1 
	[33] scales_1.2.1     splines_4.2.2    abind_1.4-5      colorspace_2.0-3
	[37] ggsignif_0.6.4   utf8_1.2.2       munsell_0.5.0    broom_1.0.2   


* Make sure to change the directories "~/Biometrical/..." to the suitable directory. 
* Make sure the following folder are empty for reproducibility purposes: "/Figures_Tables", "Intermediate_results/simulated_cohorts". 


Notice that:
** "Estimation_functions.R" contains all functions for parameter estimations in the case of longitudinal data.
** "Simulation_functions.R" contains all functions for data simulation.

	# Function "natural_sim_Cpp_Rseed" is used to simulate the natural history of multiple cohorts:
		*Input: p = a vector of model A parameters; study_years = a numeric vector of the year where the study starting and the year where the study ending; base_age = a numeric scalar/vector of the baseline of age; 
			prop_year = a numeric scalar/vector of the probabilities from the Categorical distribution of the base_age; left_trunc: a numeric scalar of left-truncation age; 
			id = a numeric vector of cohort identifications; hazard_distr = distribution of the hazard for model A (0 = Exponential distribution, 1 = Weibull distribution, 2 = Gompertz distribution).  
		*Output: a dataframe with columns: id, state, age, year, time (Here, time is the centering age). 

	# Function "control_sim_Cpp_Rseed" is used to simulate the control group of multiple cohorts:
		*Input: p = a vector of model A parameters; study_years = a numeric vector of the year where the study starting and the year where the study ending; base_age = a numeric scalar/vector of the baseline of age; 
			prop_year = a numeric scalar/vector of the probabilities from the Categorical distribution of the base_age; left_trunc: a numeric scalar of left-truncation age; 
			id = a numeric vector of cohort identifications; hazard_distr = distribution of the hazard for model A (0 = Exponential distribution, 1 = Weibull distribution, 2 = Gompertz distribution);
			max_no_seeds = an integer value (ususally be between 2-5).

		*Output: a list of three dataframes : natural_history, life_time_record and study_record. The format of natural_history data frame is the same with the output of the function "natural_sim_Cpp_Rseed".
			The other two dataframes have following columns: id, state, age, year, screen (dummy variable to indicate whether this individual has ever been screened or not) and time (Here, time is the centering age).

	# Function "screen_sim_Cpp_Rseed" is used to simulate the screened group of multiple cohorts:
		*Input: p = a vector of model A parameters; sub_p = a vector of model B parameters; screen_misc = a numeric vector of misclassification parameters; take_up = a scalar of proportion of compliance (between 0 and 1);
			screen_freq = a scalar of frequency of screening by years; screen_age = a numeric of vector of the age range where screening is offered; 
			study_years = a numeric vector of the year where the study starting and the year where the study ending;
			base_age = a numeric scalar/vector of the baseline of age; prop_year = a numeric scalar/vector of the probabilities from the Categorical distribution of the base_age; 
			left_trunc: a numeric scalar of left-truncation age; id = a numeric vector of cohort identifications; 
			hazard_distr = a vector of distribution of the hazard for model A (0 = Exponential distribution, 1 = Weibull distribution, 2 = Gompertz distribution); 
			hazard_sub_distr = a vector of distribution of the hazard for model B (0 = Exponential distribution, 1 = Weibull distribution, 2 = Gompertz distribution);
			max_no_seeds = an integer value (ususally be between 2-5).

		*Output: a list of three dataframes : natural_history, life_time_record and study_record. The format of natural_history data frame is the same with the output of the function "natural_sim_Cpp_Rseed".
			The other two dataframes have following columns: id, state, age, year, screen (dummy variable to indicate whether this individual has ever been screened or not) and time (Here, time is the centering age).
 

p/s: max_no_seeds is used to ensure we have enough number of seeds to ensure by the year where the study starting the simulated cohorts are alive and cancer-free. 
Higher value means larger computer resources are required, smaller value may cause R to terminate. If we set the base_age near the peak age where cancers are mostly found or deaths occuring, we need higher max_no_seeds. 





****************
REPRODUCIBILITY
****************
Proceed the following steps:

> 1. Run the R file "Install_packages.R" will install all required additional packages. 

> 2. Table 2 (Estimation for PLCO data): Run "PLCO_data_estimation_Windows.R" only give the estimation for the 3-states model because it will take extensive amount of time to get estimation for the 4-states model. 
	Instead, we recomment to run "PLCO_data_estimation_Linux.R" on a Linux platform device to enable parallel computing. 
	Here, we setted number of cores equal to 15, as our device has 16 cores. This number of core need to be adapted to the device configuration. 
	The time it takes to run "PLCO_data_estimation_Linux.R" when using 15 cores is approximately 15 mins. 

> 3. Callibration via simulation is written in  "Calibration_via_simulation.R":

> 4. Evaluation for screening strategies are written in "Screen_evaluation.R", "Screen_evaluation_long.R" and "Screen_evaluation_lead_times.R".
	We have repeat simulating a large set of cohorts for 500 times in "Screen_evaluation_long.R" and "Screen_evaluation_lead_times.R", which can take up to 8 hours. 
	All results are saved in the "Intermediate_results\screen_evaluation_500k" folder. 
	For comparison purpose, we recommend to run only ""Screen_evaluation.R" to compare the first few runs from the results we provided. 
	The result can be checked from the file "Compare_screen_evaluation.txt" and "Compare_screen_evaluation_lead_time.txt" from the "/Figures_Tables" folder.

> 5. All figures and tables using previous runs are obtained when running "Figure_tables.R", these will output to the "/Figures_Tables" folder.



The total run time running on our Windows platform device did not exceed beyond 20 minutes (excluding the time to install additional packages). 






	  
