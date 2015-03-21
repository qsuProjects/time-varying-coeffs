# this pushes the files from local GitHub repo to Sherlock
# note: will give lots of permissions errors, but it does work
system( "scp -r ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/* mmathur@sherlock:/share/PI/manishad/tvc" )

# run them
# sbatch /share/PI/manishad/tvc/genCov/sbatch/1.sbatch -p manishad
# sbatch /share/PI/manishad/tvc/genSurv/sbatch/0001_genSurv.sbatch -p manisha

# look at genCov results
# cd /share/PI/manishad/tvc/genCov/output/generated_data

# look at genCov sbatches
# cd /share/PI/manishad/tvc/genCov/sbatch

# copy completed dataset back to local machine
system( "scp mmathur@sherlock:/share/PI/manishad/tvc/genSurv/datasets/* ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data " )