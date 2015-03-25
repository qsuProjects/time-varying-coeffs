
########################### LOAD COMMAND-LINE ARGUMENTS ###########################

# load command line arguments
#args = commandArgs(trailingOnly = TRUE)
#print(args)

# EDIT THIS :)
# see below (function definition) for meaning of each
# model.type = args[1]
# data.path = args[2]

# needed arguments:
#  logical for are we doing right or wrong model?
#  path to specific dataset (split data for right model; true data for wrong model)
#  results.write.path: for model results
# path to true betas so that we can compare the models

# read in true data
#d.true = read.csv(data.path)

# from the dataset path, get vector of names of all the datasets
#dataset.paths = 


library(survival)


# NOT TESTED - GULP!!
# get model fit results for each dataset
l_ply( dataset.paths, .parallel=T, function(.item, type, right.model.form.path, wrong.model.form.path, results.write.path, name.prefix) {  
  # each element of dataset.names is passed in parallel to function() as .item

  # write the appropriate line into its own file
  fit_one_model( .dataset.path=.item, .type=type, .right.model.form.path=right.model.form.path,
                 .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path, .name.prefix=name.prefix )
  
}, type, right.model.form.path, wrong.model.form.path, results.write.path, name.prefix)



########################### FUNCTION: COLLAPSE "TRUE" DATA TO "OBSERVED" DATA ########################### 

# LOCAL TEST - WORKS :)
#data = read.csv("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/SURV_2015-03-12_dataset.csv")

# keep only rows with event
collapse_data = function(data, id.var.name="id", event.var.name="d") {
  d2 = data[ data[[event.var.name]] == 1, ]
  return(d2)
}



################### OBSOLETE??
########################### FUNCTION: INITIALIZE RESULTS FILE ###########################

# create a file to store results
# rows = 2 * number of datasets generated (1 row for right model, 1 for wrong)
# prepopulate a column for model type (right vs. wrong)

# columns:
#  dataset name (prepopulated to avoid overwriting)
#  model type (right vs. wrong)
#  bias
#  MSE
#  coverage (0/1)

# arguments:
#  results.write.path
#  n.datasets: number of datasets that were generated

init_results = function(.dataset.names, .model.type) {
  
  r = as.data.frame( matrix(nrow=length(.dataset.names), ncol=6) )
  names(r) = c("dataset.name", "date", "model.type", "bias", "mse", "coverage")
  
  r$dataset.name = .dataset.names
  r$date = Sys.Date()
  r$model.type = .model.type
  
  return(r)
}


########################### FUNCTION: FIT MODEL ###########################

# given a dataset, fit the specified model (right or wrong) and write results to a single csv file

# .dataset.path: path to the specific dataset we are going to use
# .type: "right" or "wrong" model
# .right.model.form.path: path to a text file holding formula for right model
# .wrong.model.form.path: path to a text file holding formula for wrong model
# .name.prefix

fit_one_model = function( .dataset.path, .type, .right.model.form.path, .wrong.model.form.path, .results.write.path, .name.prefix ) {
  
  # read in dataset
  data = read.csv(.dataset.path)
  
  # based on the arguments, read in the Cox formula
  formula = readLines( switch( .type, "right"=.right.model.form.path, "wrong"=.wrong.model.form.path ) )

  # fit model
  rs = coxph( eval( parse(text=formula) ), data=data )
  coef = rs$coefficients
  CI.low = summary(rs)$conf.int[,3]
  CI.high = summary(rs)$conf.int[,4]
    
  # return the row for the model results
  vec = c( as.character(Sys.Date()), coef, CI.low, CI.high, .dataset.path )
  row = data.frame( matrix(nrow=1, ncol=length(vec)) )
  row[1,] = vec
  
  # improve names of row
  names(row) = c( "date.completed", paste( names( coef ), "_coef", sep="" ), 
                  paste( names( CI.low ), "_lowCI", sep="" ),  paste( names( CI.high ), "_highCI", sep="" ),
                  "dataset" )
  
  # write the row
  file.name = .name.prefix  # EDIT THIS LATER 
  write.csv(row, paste(.results.write.path, file.name, sep="/") )
}


# LOCAL TEST - WORKS :)
#.dataset.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/collapsed_data.csv"
#.type = "wrong"
#.right.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/right_model_formula.txt"
#.wrong.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/wrong_model_formula.txt"
#.results.write.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/results"
















