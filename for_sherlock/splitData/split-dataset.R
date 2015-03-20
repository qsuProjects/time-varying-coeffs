#TESTTTTTT

########################### SPLIT DATASET ON FOLLOW-UP TIMES ########################### 


##### Function: Split a dataframe of survival times by event or follow-up times. #####
# assumes there is only one row per subject
# returns a dataset with all the original variables plus start.time and stop.time.

# d = data frame
# id.var = name of id variable
# followup.var = name of followup variable
# event.var = name of binary event variable (must be coded 0/1)
# split.at = "followup" to split at followup times or "event" to split at event times 
# na.rm = remove subjects with missing followup or event variables?

event_split = function(d, id.var, followup.var, event.var, split.at, na.rm=TRUE) {
  
  # remove NAs in followup and event variables if specified
  # droplevels to get rid of unused levels in subsetted data
  if (na.rm) d = droplevels( d[ !is.na(d[[followup.var]]) & !is.na(d[[event.var]]), ] )
  
  #browser()
  
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
  
  return( list(split.data=d2, split.times=time.vec) ) 
}



###### Function (internal): Given a row from dataframe, expand it on provided splitting times. #####

expand_subject = function(subject, time.vec, id.var, followup.var, event.var, split.at) {
  
  #browser()
  
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


# toy example
( fake = data.frame( id=c(1,2,3), dosage=c(50,25,100), followup=c(2,4,1), event=c(1,0,1)) )
( split = event_split(d=fake, id.var="id", followup.var="followup", event.var="event", split.at="followup", na.rm=TRUE) )
split$cumulative.dosage = split$dosage * split$stop.time



