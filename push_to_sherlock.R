# this pushes the files from local GitHub repo to Sherlock
# note: will give lots of permissions errors, but it does work
system( "scp -r ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/* mmathur@sherlock:/share/PI/manishad/tvc/genSurv" )

# run them
# sbatch /share/PI/manishad/tvc/genSurv/sbatch/0001_genSurv.sbatch

# copy completed dataset back to local machine
system( "scp mmathur@sherlock:/share/PI/manishad/tvc/genSurv/datasets/* ~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data " )