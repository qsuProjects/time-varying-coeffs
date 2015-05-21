
########################### LOAD COMMAND-LINE ARGUMENTS ###########################

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
print(args)

# see below (function defintion) for meaning of each
true.data.path = args[1]
write.path.collapsed = args[2]
write.path.split = args[3]
name.prefix = args[4]



#### FOR LOCAL TEST ONLY :) ####
true.data.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test/SURV_2015-04-29_job_1_covariates.csv"
write.path.collapsed = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test"
write.path.split = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test"
name.prefix = "name.prefix"

true.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test/true_model_formula.txt"
right.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test/right_model_formula.txt"
wrong.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test/wrong_model_formula.txt"
noninteractive.model.form.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test/noninteractive_model_formula.txt"

results.write.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test/output"
################################

# read in true data
td = read.csv(true.data.path)

library(survival)
  

########################### FUNCTION: SPLIT DATASET ON FOLLOW-UP TIMES ########################### 

##### Split a dataframe of survival times by event or follow-up times. #####
# assumes there is only one row per subject
# returns a list with two elements: split.data (dataset with all the original variables)
#  and split.times (a vector of the times at which data were split)

# d = data frame
# id.var = name of id variable
# followup.var = name of followup variable
# event.var = name of binary event variable (must be coded 0/1)
# split.at = "followup" to split at followup times or "event" to split at event times 
# na.rm = remove subjects with missing followup or event variables?

event_split = function(d, id.var, followup.var, event.var, split.at, na.rm=TRUE,
                       name.prefix="", write.path) {
  
  # remove NAs in followup and event variables if specified
  # droplevels to get rid of unused levels in subsetted data
  if (na.rm) d = droplevels( d[ !is.na(d[[followup.var]]) & !is.na(d[[event.var]]), ] )

  # make vector of splitting times based on whether user specifies splitting at followup... 
  # ...or event times.
  if (split.at=="followup") time.vec = c(0, sort( unique(d[[followup.var]] ) ) )  # find and sort unique followup times
  if (split.at=="event") {
    temp = d[ d[[event.var]]==1, ]  # subset dataframe to only subjects with event
    time.vec = c(0, sort( unique( temp[[followup.var]] ) ) ) 
  }
  
  # split dataframe into a list where each list element is a different subject (row)
  split.data = split(d, d[[id.var]])
  
  d1 = lapply(split.data, function(x) expand_subject(x, time.vec, id.var, followup.var, event.var, split.at) )
  d2 = do.call( "rbind", lapply( d1, function(x) data.frame( as.list(x) ) ) )
  
  
  # write data file and splitting times
  if (!is.na(write.path)) {
    string1 = paste(Sys.Date(), name.prefix, "split_data.csv", sep="_")
    write.csv(d2, paste(write.path, string1, sep="/"))
    
    string2 = paste(Sys.Date(), name.prefix, "split_times.csv", sep="_")
    write.csv(time.vec, paste(write.path, string2, sep="/"))
  }
  
  return( list(split.data=d2, split.times=time.vec) ) 
}



########################### FUNCTION: EXPAND SINGLE ROW ON SPLIT TIMES ########################### 

###### Function (internal): Given a row from dataframe, expand it on provided splitting times. #####

expand_subject = function(subject, time.vec, id.var, followup.var, event.var, split.at) {
  followup = subject[[followup.var]]
  time.vec.2 = time.vec[time.vec <= followup]  # remove any times that are after the subject's followup
  if (split.at == "event") time.vec.2 = unique( c(time.vec.2, followup) )  # if subject had event, add a final record ending at the end of their followup time
  
  n.records = length(time.vec.2) - 1 # how many records will this subject have?
  
  # initialize dataframe for this subject's expanded rows
  n.Vars = ncol(subject)
  expanded = as.data.frame( matrix(nrow=n.records, ncol=(n.Vars + 2) ) )
  names(expanded) = c( names(subject), "start.time", "stop.time")
  
  for( i in 1:n.records ) {
    start.time = time.vec.2[i]
    stop.time = time.vec.2[i+1]
    
    expanded[i,] = cbind(subject, start.time, stop.time)
  }
  
  # for subjects who had the event, entries should be 0 event until the subject's last entry
  had.event = (subject[[event.var]] == 1)  # did the subject have the event?
  if (had.event) {
    # change the event variable to 0 wherever the stop time is not the maximum stop time
    expanded[[event.var]][ expanded$stop.time != max(expanded$stop.time) ] = 0
  }
  return(expanded)
}


########################### FUNCTION: COLLAPSE "TRUE" DATA TO "OBSERVED" DATA ########################### 

# LOCAL TEST - WORKS :)
#data = read.csv("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/SURV_2015-03-12_dataset.csv")

# keep only rows with event
collapse_data = function(data, id.var.name="id", event.var.name="d", name.prefix="", write.path=NA) {
  d2 = data[ data[[event.var.name]] == 1, ]
  
  # write data file
  if (!is.na(write.path)) {
    string = paste(Sys.Date(), name.prefix, "collapsed_data.csv", sep="_")
    write.csv(d2, paste(write.path, string, sep="/"))
  }
  return(d2)
}

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
file.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data/split_dataset_1.csv"
name.prefix = "right_results"
update_name(file.path, name.prefix)



########################### FUNCTION: FIT MODEL ###########################

# given a dataset, fit the specified model (right or wrong) and write results to a single csv file

# .data: dataset name (e.g., "dataset_1")
# .type: "right", "wrong", or "true" model
# .right.model.form.path: path to a text file holding formula for right model
# .wrong.model.form.path: path to a text file holding formula for wrong model
# .name.prefix

fit_one_model = function( .data, .data.name, .type, .noninteractive.model.form.path, .true.model.form.path,
                          .right.model.form.path, .wrong.model.form.path, .results.write.path ) {

  # based on the arguments, read in the Cox formula
  form = readLines( switch( .type, "noninteractive"=.noninteractive.model.form.path,
                            "right"=.right.model.form.path,
                            "true"=.true.model.form.path,
                            "wrong"=.wrong.model.form.path ) )
  
  # fit model
  rs = coxph( eval( parse(text=form) ), data=.data )
  coef = rs$coefficients
  CI.low = confint(rs)[,1]
  CI.high = confint(rs)[,2]
  #CI.low = summary(rs)$conf.int[,3]
  #CI.high = summary(rs)$conf.int[,4]

  # PH violation
  PH = cox.zph(rs, transform="log")$table
  PH.dose.p = PH["start.dose","p"]
  
  if (nrow(PH) > 1) {  # if there's only 1 coefficient, it doesn't report a global p-value, so just report dose p again
    PH.global.p = PH["GLOBAL","p"]
  } else {
    PH.global.p = PH.dose.p
  }
  
  
  #browser()
  
  # return the row for the model results
  vec = c( as.character(Sys.Date()), as.character(.type), coef, CI.low, CI.high, PH.global.p, PH.dose.p, .data.name, form )
  row = data.frame( matrix(nrow=1, ncol=length(vec)) )
  row[1,] = vec
 
  # improve names of row
  names(row) = c( "date.completed", "model", paste( names( coef ), "_coef", sep="" ),
                  paste( names( coef ), "_lowCI", sep="" ),  paste( names( coef ), "_highCI", sep="" ),
                  "PH.global.p", "PH.dose.p",
                  "dataset", "form" )
  
  # write the row and append descriptive name prefix
  #type.name = switch( .type, "true"="true_results", "right"="right_results", "wrong"="wrong_results")
  type.name = paste( .type, "_results", sep="" )
  file.name = paste(.data.name, "_", type.name, ".csv", sep="")
  write.csv(row, paste(.results.write.path, file.name, sep="/") )
}


########################### RUN ########################### 

# make collapsed data
cd = collapse_data(data=td, id.var.name="id", event.var.name="d",
              name.prefix=name.prefix, write.path=write.path.collapsed)

# make split data
sd = event_split(d=cd, id.var="id", followup.var="t", event.var="d", split.at="followup", 
            name.prefix=name.prefix, write.path=write.path.split)$split.data

# fit true model
fit_one_model( .data=td, .data.name="data.name", .type="true", .noninteractive.model.form.path=noninteractive.model.form.path,
               .true.model.form.path=true.model.form.path,
                 .right.model.form.path=right.model.form.path,
           .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path )

# fit right model
fit_one_model( .data=sd, .data.name="data.name", .type="right", .noninteractive.model.form.path=noninteractive.model.form.path,
               .true.model.form.path=true.model.form.path,
               .right.model.form.path=right.model.form.path,
               .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path )

# fit wrong model
fit_one_model( .data=cd, .data.name="data.name", .type="wrong", .noninteractive.model.form.path=noninteractive.model.form.path,
               .true.model.form.path=true.model.form.path,
               .right.model.form.path=right.model.form.path,
               .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path )


# fit noninteractive model
fit_one_model( .data=cd, .data.name="data.name", .type="noninteractive", .noninteractive.model.form.path=noninteractive.model.form.path,
               .true.model.form.path=true.model.form.path,
               .right.model.form.path=right.model.form.path,
               .wrong.model.form.path=wrong.model.form.path, .results.write.path=results.write.path )



