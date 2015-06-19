
########################### SIMULATE ENTIRE DATASET ###########################

library(ggplot2)
library(survival)

# simulate dataset with just one row per S  
setwd("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/genCov/r")
source("simulate.R")
  
n = 500
doses = c(0, 1)
doses.p = rep( 1/length(doses), length(doses) )
max.t = 1
seed = 1

setwd("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test")
source("split-dataset.R")

x = simulate_Xs( n, doses, doses.p, max.t, seed, write.path )

# remove unneeded columns
x = x[,-c(3:5)]


########################### GENERATE SURVIVAL TIMES ###########################

# try generating survival times
# true event time
td = x
td$surv.time=NA

# Weibull
# use Schemper's scenario C
lambda0 = 2
lambda1 = 1
nu0 = 1
nu1 = 2
td$surv.time[td$start.dose==0] = rweibull(n=length(td$start.dose[x$start.dose==0]), shape=nu0, scale=lambda0) 
td$surv.time[td$start.dose==1] = rweibull(n=length(td$start.dose[x$start.dose==1]), shape=nu1, scale=lambda1) 

# exponential
#lambda0=.5
#lambda1=1
#td$surv.time[td$start.dose==0] = rexp(n=length(td$start.dose[x$start.dose==0]), rate=lambda0)
#td$surv.time[td$start.dose==1] = rexp(n=length(td$start.dose[x$start.dose==1]), rate=lambda1)

# everyone has event
td$event=1

#td$surv.time = td$surv.time * 10  # in order not to have too many that are < 1



########################### LOOK DESCRIPTIVELY AT SIMULATED DATA ###########################

# look at distribution of survival times
ggplot(data=td, aes(x=surv.time)) + geom_density(fill="grey") + theme_bw() +
  ggtitle("") + 
  scale_x_continuous(limits=c(0,6))

aggregate(surv.time ~ start.dose, td, summary )


# KM plot
m1 = survfit( Surv(surv.time) ~ as.factor(start.dose), data=td)
colors = c("red", "#2ca25f")
plot(m1, xlim=c(0,3), ylim=c(0,1), conf.int=FALSE, mark.time=FALSE, col=colors, xlab="Days of follow-up",
     ylab="Cumulative survival probability", lwd=3, main="")
legend("topright", c("Dose=0", "Dose=1"), lty=1, lwd=3, col=colors)


# plot hazards
library(muhaz)
haz0 = with( td[td$start.dose==0,], muhaz(surv.time))
haz1 = with( td[td$start.dose==1,], muhaz(surv.time))

plot(haz0, xlim=c(0,3), ylim=c(0,6), lwd=2, col=colors[2], xlab="Time")
lines(haz1, lwd=2, col=colors[1])

legend("topleft", legend=c("Dose=0", "Dose=1"), 
       lty=1, col=colors, lwd=2)

# yes, looks like in Schemper!

# plot HR
hr = haz1$haz.est / haz0$haz.est
plot(hr, ylim=c(0,10))



########################### FIT THE MODELS ###########################

# true average HR (per Schemper Table 1): 1.93

##### Noninteractive Cox ######
rs.n = coxph( Surv(surv.time) ~ start.dose, data=td ); summary(rs.n)
(a = cox.zph(rs.n) )
# violates PH :)
# estimated HR: 2.95; pretty bad


##### Wrong Cox ######
rs.w = coxph( Surv(surv.time) ~ 
              start.dose + start.dose:surv.time, data=td ); summary(rs.w)
(a = cox.zph(rs.w) )
# violates PH
# estimated HR: completely crazy


##### Right Cox ######

# make split data using Therneau function
#split.tt = tmerge(td, td, id=td$id, tstop=surv.time)
# did not work - need to try again

# make split data (split on event times) using my function
split.mm.2 = event_split(d=td, id.var="id", followup.var="surv.time", event.var="event", split.at="event", 
                       name.prefix="name.prefix", write.path=NA)$split.data
dim(split.mm.2)

# right model with frailties
rs.r = coxph(Surv(start.time, stop.time, event) ~ start.dose + start.dose:start.time + frailty(id), data=split.mm.2 ); summary(rs.r)
(a = cox.zph(rs.r) )
# does NOT violate PH :D





# with frailties; no interaction
rs.temp = coxph(Surv(start.time, stop.time, event) ~ start.dose + frailty(id), data=split.mm.2 ); summary(rs.temp)
# still off by same amount

# pooled Cox; interaction
rs.temp = coxph(Surv(start.time, stop.time, event) ~ start.dose + start.dose:start.time + cluster(id), data=split.mm.2 ); summary(rs.temp)




