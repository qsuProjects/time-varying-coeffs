
########################### FUNCTION: UPDATE NAME OF DATASET ###########################

# given a file path, extract the name of the dataset and update it

update_name = function(.file.path, .name.prefix) {
  # extract the part of the dataset name that starts with "dataset"
  start = regexpr("dataset", .file.path, fixed=TRUE)
  dataset.name = substr(.file.path, start=start, stop=nchar(.file.path) )
  
  # add the name prefix
  return( paste(.name.prefix, dataset.name, sep="_") )
}

# TEST - WORKS :)
#file.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/split_dataset_1.csv"
#name.prefix = "right_results"
#update_name(file.path, name.prefix)


########################### FUNCTION: COLLAPSE "TRUE" DATA TO "OBSERVED" DATA ########################### 

# LOCAL TEST - WORKS :)
#data = read.csv("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/SURV_2015-03-12_dataset.csv")

# keep only rows with event
collapse_data = function(data, id.var.name="id", event.var.name="d") {
  d2 = data[ data[[event.var.name]] == 1, ]
  return(d2)
}


########################### FUNCTION: FIT MODEL ###########################

# given a dataset, fit the specified model (right or wrong) and write results to a single csv file

# .dataset.path: path to the specific dataset we are going to use
# .type: "right" or "wrong" model
# .right.model.form.path: path to a text file holding formula for right model
# .wrong.model.form.path: path to a text file holding formula for wrong model
# .name.prefix

fit_one_model = function( .dataset.path, .type, .right.model.form.path, .wrong.model.form.path, .results.write.path ) {
  # read in dataset
  # if fitting right model, assume the dataset path is already the split data that we're going to use
  # if fitting wrong model, assume the dataset path is the true data that must first be collapsed
  data = switch( .type, "right"=read.csv(.dataset.path),
                 "wrong"=collapse_data( read.csv(.dataset.path), id.var.name="id", event.var.name="d" ) )
  
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
  
  # write the row and append descriptive name prefix
  name.prefix = switch( .type, "right"="right_results", "wrong"="wrong_results")
  file.name = update_name(.dataset.path, name.prefix)
  write.csv(row, paste(.results.write.path, file.name, sep="/") )
}


##### LOCAL TEST - WORKS :)
#.dataset.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/SURV_2015-03-12_dataset_1.csv"
#.type = "wrong"
#.right.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/right_model_formula.txt"
#.wrong.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/wrong_model_formula.txt"
#.results.write.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/results"

#.dataset.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/split_dataset_1.csv"
#.type = "right"
#.right.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/right_model_formula.txt"
#.wrong.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/wrong_model_formula.txt"
#.results.write.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/results"

#fit_one_model( .dataset.path, .type, .right.model.form.path, .wrong.model.form.path,
#                          .results.write.path )











