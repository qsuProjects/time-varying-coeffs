
########################### LOAD COMMAND-LINE ARGUMENTS ########################### 

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
print(args)

# see below (function defintion) for meaning of each
n = as.numeric(args[1])
dose.param.path = args[2]
max.t = as.numeric(args[3])
seed = as.numeric(args[4])
name.prefix = args[5]
write.path = args[6]


########################### READ IN DOSE PARAMETERS ########################### 

dose.params = read.csv(dose.param.path, head=TRUE); print(dose.params)
doses = dose.params$dose
doses.p = dose.params$doses.p

########################### FUNCTION: SIMULATE COVARIATES ########################### 

# n: number of subjects
# doses: vector of possible treatment doses administered at beginning of study
# doses.p: vector of probabilities of each dose
# max.t: maximum observation time (will truncate later)
# seed: seed value for simulation


simulate_Xs = function(n, doses, doses.p, max.t, seed=NA, name.prefix="", write.path=NA) {
  
  # initialize dataframe
  d = data.frame(matrix(ncol = 5, nrow = n*max.t))
  names(d) = c("id", "start.dose", "t0", "t", "dose.X.t")
  d$id = rep(1:n, each=max.t)
  
  # make start and stop times
  d$t0 = rep( 0:(max.t-1), n )
  d$t = rep( 1:max.t, n )
  
  # draw start doses for each subject
  if (!is.na(seed)) set.seed(seed)
  d$start.dose = rep( sample(doses, size=n, replace=TRUE, prob=doses.p), each=max.t )
  
  # compute current dose (interaction variable)
  d$dose.X.t = d$start.dose * d$t

  # write data file
  if (!is.na(write.path)) {
    string = paste(Sys.Date(), name.prefix, "covariates.csv", sep="_")
    write.csv(d, paste(write.path, string, sep="/"))
  }
  return(d)
}



########################### RUN ########################### 

d = simulate_Xs( n=n, doses=doses, doses.p=doses.p, max.t=max.t, seed=seed, name.prefix=name.prefix, write.path=write.path )



########################### LOCAL TEST ########################### 
# WORKS :)
#n = 1000
#doses = c(1, 2, 3)
#doses.p = rep( 1/length(doses), length(doses) )
#max.t = 20
#seed = 1
#write.path = "/Users/mmathur/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for-sherlock"

#d = simulate_Xs( n, doses, doses.p, max.t, seed, write.path )



