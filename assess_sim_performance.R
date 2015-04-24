
###################### READ IN DATA ###################### 

# read in file with survival times copied from Sherlock
setwd("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/test_data")
d = read.csv("SURV_2015-03-12_dataset.csv")

# number of subjects
length(unique(d$id))

library(survival)


###################### SEE IF WE'RE HITTING INTENDED BETAS ###################### 

rs1 = coxph( Surv(t0, t, d) ~ start.dose + dose.X.t + cluster(id), data=d ); summary(rs1)

rs2 = coxph( Surv(t0, t, d) ~ start.dose + cluster(id), data=d ); summary(rs2)
