
########################### LOAD COMMAND-LINE ARGUMENTS ###########################

# load command line arguments
#args = commandArgs(trailingOnly = TRUE)
#print(args)

# EDIT THIS :)
# see below (function definition) for meaning of each
n.clusters = as.numeric(args[1])
model.type = args[2]
data.path = args[3]  # path to folder with all datasets
results.write.path = args[4]
functions.path = args[5]


############# LOCAL TEST
n.clusters=2
functions.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/fitModels/r/fit_models.R"

type="right"
right.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/right_model_formula.txt"
wrong.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/wrong_model_formula.txt"
results.write.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/results"
data.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/input_data"
#####################



# from the dataset path, get vector of names of all the datasets
dataset.paths = list.files(data.path, full.names=TRUE)

# load libraries
library(plyr)
library(doSNOW)


########################### FIT MODELS IN PARALLEL ###########################


# set up parallelization backend
cl = makeCluster(n.clusters)
registerDoSNOW(cl)

# these should be equal
n.clusters
getDoParWorkers()

#give each worker an ID 
#this is optional, but occasionally is useful 
#in situations where you need workers to write data
#every worker writing to the same file will result in 
#uncaught errors and partially overwritten data
#worker ids enable each worker to write to its own set of 
#files
#clusterApply(cl, seq(along=cl), function(id) WORKER.ID <<- paste0("worker_", id))

########## PARALLELIZATION NOT YET TESTED!

l_ply( dataset.paths, .parallel=T, function(.item, type, right.model.form.path, wrong.model.form.path, results.write.path, functions.path) {  
  # each element of dataset.paths is passed in parallel to function() as .item
  
  source(functions.path)
  
  # TEST
  #num = sample(0:100, 1)
  #write.csv( num, paste("~/Desktop/", .item, "_test.csv", sep="") )
  
  library(survival)
  
  # TEST
  fit_one_model( .dataset.path=.item, .type=type, .right.model.form.path=right.model.form.path,
                 .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path )
  
  
  # write the appropriate line into its own file
 # fit_one_model( .dataset.path=.item, .type=type, .right.model.form.path=right.model.form.path,
  #              .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path )
  
}, type, right.model.form.path, wrong.model.form.path, results.write.path, functions.path)




