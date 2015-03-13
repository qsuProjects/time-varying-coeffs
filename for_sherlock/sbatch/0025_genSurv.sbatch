#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=gs0025
#################  
#a file for job output, you can check job progress
#SBATCH --output=ogs0025.out
#################
# a file for errors from the job
#SBATCH --error=egs0025.err
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
#SBATCH  --mail-user=kikapp@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=1
#################
#SBATCH --cpus-per-task=8
#now run normal batch commands

ml load R
srun R -f /share/PI/manishad/PCORI/genSurv/r/gSurvMain.R --args  /share/PI/manishad/PCORI/genSurv/r/generateTimeAndFit.R /share/PI/manishad/PCORI/genSurv/r/betas.R /share/PI/manishad/sim_23_sher/ 8 50 25 /scratch/PI/manishad/PCORI/genSurv/output/modelFit survtime.csv /scratch/PI/manishad/PCORI/genSurv/output/log/ /scratch/PI/manishad/PCORI/genSurv/output/covariates/
