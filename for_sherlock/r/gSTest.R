# scp ~/shexport/PCORI/generateByScenario.R kikapp@sherlock:~/PCORI/lib
args <- commandArgs(trailingOnly = T)
print(args)

source_location <- args[1]
beta_location <- args[2]
covariate_data_path <- args[3]
number_cores <- as.numeric(args[4])
range_low <- as.numeric(args[5])
range_high <- as.numeric(args[6])
results_write_directory <- args[7]
sim_results_name <- args[8]
log_write_directory <- args[9]

# 
# source_location <- "~/shexport/PCORI/genSurv/r/generateTimeAndFit.R"
# beta_location <- "~/shexport/PCORI/genSurv/r/betas.R"
# covariate_data_path <- "~/shexport/PCORI/genSurv/data/data_1000.csv"
# number_cores <- as.numeric("2")
# range_low <- as.numeric("1")
# range_high <- as.numeric("2")
# results_write_directory <- "~/shexport/PCORI/genSurv/output/test"
# sim_results_name <- "testResults.csv"
# log_write_directory <- "~/shexport/PCORI/PCORI/genSurv/log/test"



# source_location <- "/share/PI/manishad/PCORI/genSurv/r/generateTimeAndFit.R"
# beta_location <- "/share/PI/manishad/PCORI/genSurv/r/betas.R"
# covariate_data_path <- "/scratch/PI/manishad/PCORI/genSurv/covariates/data_1000.csv"
# number_cores <- as.numeric("2")
# range_low <- as.numeric("1")
# range_high <- as.numeric("1000")
# results_write_directory <- "/scratch/PI/manishad/PCORI/genSurv/output/test"
# sim_results_name <- "testResults.csv"
# log_write_directory <- "/scratch/PI/manishad/PCORI/genSurv/log/test"

# 
# if(is.na(output_location)) {output_location <- "/scratch/users/kikapp/dump"}
# if(is.na(n_cores)) {n_cores <- 16}
# if(is.na(range_low)) {range_low <- 0}
# if(is.na(range_high)) {range_high <- 5}
# if(is.na(source_location)) {source_location <- "~/PCORI/lib"}
# if(is.na(scenario_location)) {scenario_location <- "~/PCORI/lib"}
# if(is.na(scenario_row)) {scenario_row <- 1}

#load betas
source(beta_location)
print(beta)

#initialize seed
the_seed <- 5

source(source_location)
