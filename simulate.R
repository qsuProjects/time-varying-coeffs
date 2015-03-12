
########################### FUNCTION: SIMULATE COVARIATES ########################### 

# n: number of subjects
# doses: vector of possible treatment doses administered at beginning of study
# doses.p: vector of probabilities of each dose
# max.t: maximum observation time (will truncate later)
# seed: seed value for simulation

simulate_Xs = function(n, doses, doses.p, max.t, seed=NA) {
  
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
  
  return(d)
}

# TEST - works :)
n = 1000
doses = c(1, 2, 3)
doses.p = rep( 1/length(doses), length(doses) )
max.t = 20
seed = 1
d = simulate_Xs( n, doses, doses.p, max.t, seed )






