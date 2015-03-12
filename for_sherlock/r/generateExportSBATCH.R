source("/Volumes/QSU/Datasets/sherlock/r/generateSbatchPSwitch.r")
library(stringr)


n_jobs <- 50
jobs <- expand.grid(segnum = 1:n_jobs,
                                  jobtime = "1:00:00",
                                  quality = "normal",
                                  node_number = 1,
                                  mem_per_node = 16000,
                                  mailtype =  "ALL",
                                  user_email = "kikapp@stanford.edu",
                                  tasks_per_node = 1,
                                  cpus_per_task = 8,
                                  path_to_r_script = "/share/PI/manishad/PCORI/genSurv/r/gSurvMain.R",
                                  stringsAsFactors = F)

source_location <- paste0("/share/PI/manishad/PCORI/genSurv/r/generateTimeAndFit.R")
beta_location <- paste0("/share/PI/manishad/PCORI/genSurv/r/betas.R")
covariate_data_path <- paste0("/share/PI/manishad/sim_23_sher/")
n_cores <- 8
n_segments <- n_jobs
segment_index <- jobs$segnum
results_write_directory <- paste0("/scratch/PI/manishad/PCORI/genSurv/output/modelFit")
sim_results_name <- "survtime.csv"
log_write_directory <- paste0("/scratch/PI/manishad/PCORI/genSurv/output/log/")
data_write_directory <- paste0("/scratch/PI/manishad/PCORI/genSurv/output/covariates/")

jobs$jobname <- paste0("gs", str_pad(jobs$segnum, 4, side = "left", pad = "0"))
jobs$outfile <- paste0("ogs", str_pad(jobs$segnum, 4, side = "left", pad = "0"), ".out")
jobs$errorfile <- paste0("egs", str_pad(jobs$segnum, 4, side = "left", pad = "0"), ".err")

jobs$args_to_r_script <- paste("--args ", source_location, beta_location, covariate_data_path,
                               n_cores, n_segments, segment_index, results_write_directory,
                               sim_results_name, log_write_directory, data_write_directory, sep = " ")

jobs$write_path <- paste0("~/shexport/PCORI/genSurv/sbatch/", 
                          str_pad(jobs$segnum, 4, side = "left", pad = "0"), 
                          "_genSurv.sbatch")
jobs$server_sbatch_path =  paste0("/share/PI/manishad/PCORI/genSurv/sbatch/", 
                                  str_pad(jobs$segnum, 4, side = "left", pad = "0"), 
                                  "_genSurv.sbatch")

runfile_path <-  paste0("~/shexport/PCORI/genSurv/r/sbatchRunGenSurv.R")
files <- generateSbatch(jobs, runfile_path = runfile_path, run_now = F)


#transfer files to sherlock (must have kinit credentials)
system("scp ~/shexport/PCORI/genSurv/sbatch/*.sbatch kikapp@sherlock:/share/PI/manishad/PCORI/genSurv/sbatch/")
system("scp ~/shexport/PCORI/genSurv/r/*.R kikapp@sherlock:/share/PI/manishad/PCORI/genSurv/r/")
