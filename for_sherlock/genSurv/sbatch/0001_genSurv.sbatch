#!/bin/bash
#################
#set partition 
#SBATCH -p manishad
#set a job name  
#SBATCH --job-name=gs0001
#################  
#a file for job output, you can check job progress
#SBATCH --output=ogs0001.out
#################
# a file for errors from the job
#SBATCH --error=egs0001.err
#################
#time you think you need; default is one hour
#SBATCH --time=1:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#SBATCH --mem=16000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=mmathur@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=1
#################
#SBATCH --cpus-per-task=8
#now run normal batch commands

ml load R

-- arg 1: source location
-- arg 2: beta location
-- arg 3: covariate data path
-- arg 4: number cores
-- arg 5: number segments
-- arg 6: segment index
-- arg 7: results write directory
-- arg 8: sim results name
-- arg 9: log write directory
-- arg 10: data write directory

srun R -f /share/PI/manishad/tvc/genSurv/r/gSurvMain.R --args  /share/PI/manishad/tvc/genSurv/r/generateTimeAndFit.R /share/PI/manishad/tvc/genSurv/r/betas.R /share/PI/manishad/tvc/genSurv/input_data 8 50 1 /scratch/PI/manishad/tvc/genSurv/output/modelFit survtime.csv /share/PI/manishad/tvc/genSurv/output/log/ /share/PI/manishad/tvc/genSurv/output/generated_data/
