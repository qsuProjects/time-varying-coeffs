#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=job_t
#################  
#a file for job output, you can check job progress
#SBATCH --output=rm_t.out
#################
# a file for errors from the job
#SBATCH --error=rm_t.err
#################
#time you think you need; default is one hour
#SBATCH --time=0:10:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#use the QSU job partition
#SBATCH --partition=manishad
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#max is 64000
#SBATCH --mem=4000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=mmathur@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=2
#################
#SBATCH --cpus-per-task=1
#now run normal batch commands

ml load R

# arguments:
# 1. number of clusters
# 2. model type (right or wrong)
# 3. right model formula path
# 4. wrong model formula path
# 5. results write path
# 6. R functions path
# 7. path to data folder

srun R -f /share/PI/manishad/tvc/fitModels/parallelized_wrapper.R --args PUT ARGS IN!
