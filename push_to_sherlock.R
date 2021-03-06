# this pushes the files from local GitHub repo to Sherlock
# note: will give lots of permissions errors, but it does work
system( "scp -r ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/* mmathur@sherlock:/share/PI/manishad/tvc" )



######### GENCOV ######### 
# sbatch /share/PI/manishad/tvc/genCov/sbatch/1.sbatch -p manishad

# look at genCov results
# cd /share/PI/manishad/tvc/genCov/output/generated_data

# look at genCov sbatches
# cd /share/PI/manishad/tvc/genCov/sbatch

# clear folder results
# rm output/generated_data/*
# rm sbatch/rm*
  


######### GENSURV ######### 

# push just betas
system( "scp ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/genSurv/r/betas.R mmathur@sherlock:/share/PI/manishad/tvc/genSurv/r" )

# sbatch -p manishad /share/PI/manishad/tvc/genSurv/sbatch/0001_genSurv.sbatch
# sbatch /share/PI/manishad/tvc/genSurv/sbatch/0001_genSurv.sbatch

# look at genSurv sbatches
# cd /share/PI/manishad/tvc/genSurv/sbatch

# look at genSurv results
# cd /share/PI/manishad/tvc/genSurv/output/generated_data

# look at model fit
# nano /share/PI/manishad/tvc/genSurv/output/modelFit/_job_0001core_01_survtime.csv

# look at genSurv log
# cd /share/PI/manishad/tvc/genSurv/output/log

# look at generated survival times
# scp mmathur@sherlock:/share/PI/manishad/tvc/genSurv/surv_times_temp.csv ~/Desktop

# copy completed dataset back to local machine
system( "scp mmathur@sherlock:/share/PI/manishad/tvc/genSurv/output/generated_data/* ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data " )

# copy all output back to local machine
system( "scp -r mmathur@sherlock:/share/PI/manishad/tvc/genSurv/output/* ~/Desktop " )

# clear folder results
rm egs*
rm ogs*
rm output/generated_data/*
rm output/log/*
rm output/modelFit/*
rm output/surv*
  
  
######### SPLIT DATA #########   
  
# folder
# cd /share/PI/manishad/tvc/splitData

# push just sbatch
system( "scp ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/sbatch/example.sbatch mmathur@sherlock:/share/PI/manishad/tvc/splitData/sbatch" )

# copy completed dataset back to local machine
system( "scp mmathur@sherlock:/share/PI/manishad/tvc/genSurv/output/generated_data/* ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data " )

# sbatch -p manishad /share/PI/manishad/tvc/splitData/sbatch/example.sbatch







  