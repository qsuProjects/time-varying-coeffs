
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                              SIMULATE DATA                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

########################### GENERATE X ###########################

library(ggplot2)
library(survival)

n = 500
doses = c(0, 1)
doses.p = rep( 1/length(doses), length(doses) )
max.t = 1
seed = 1
scenario = "B"  # which Schemper scenario to use (C = linearly diverging hazards; B=nonlinearly diverging)


#setwd("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/for_sherlock/splitData/local-test")
#source("split-dataset.R")

x = simulate_Xs( n, doses, doses.p, max.t, seed, write.path )

# remove unneeded columns
x = x[,-c(3:5)]


########################### GENERATE SURVIVAL TIMES ###########################

# try generating survival times
# true event time
td = x
td$surv.time=NA

# Weibull
# use Schemper's scenario
if (scenario=="C") {
  lambda0 = 2
  nu0 = 1
  nu1 = 2
  lambda1 = 1
}

if (scenario=="B") {
  lambda0 = 1.5
  nu0 = 3
  nu1 = 1
  lambda1 = 1
}

td$surv.time[td$start.dose==0] = rweibull(n=length(td$start.dose[x$start.dose==0]), shape=nu0, scale=lambda0) 
td$surv.time[td$start.dose==1] = rweibull(n=length(td$start.dose[x$start.dose==1]), shape=nu1, scale=lambda1) 

# everyone has event
td$event=1



########################### LOOK DESCRIPTIVELY AT SIMULATED DATA ###########################

##### Look at Distribution of Survival Times #####
ggplot(data=td, aes(x=surv.time)) + geom_density(fill="grey") + theme_bw() +
  ggtitle("") + 
  scale_x_continuous(limits=c(0,6))

# look at median survival times by dose
aggregate(surv.time ~ start.dose, td, summary )

# does it match theoretical medians?
lambda0*log(2)^(1/nu0)  # group 0
lambda1*log(2)^(1/nu1)  # group 1
# yes

###### KM Plot #####
m1 = survfit( Surv(surv.time) ~ as.factor(start.dose), data=td)
colors = c("red", "#2ca25f")
plot(m1, xlim=c(0,3), ylim=c(0,1), conf.int=FALSE, mark.time=FALSE, col=colors, xlab="Time",
     ylab="Cumulative survival probability", lwd=3, main="")
legend("topright", c("Dose=0", "Dose=1"), lty=1, lwd=3, col=colors)


##### Plot Hazards #####
library(muhaz)
haz0 = with( td[td$start.dose==0,], muhaz(surv.time))
haz1 = with( td[td$start.dose==1,], muhaz(surv.time))

plot(haz0, xlim=c(0,3), ylim=c(0,6), lwd=2, col=colors[2], xlab="Time")
lines(haz1, lwd=2, col=colors[1])

legend("topleft", legend=c("Dose=0", "Dose=1"), 
       lty=1, col=colors, lwd=2)

# yes, looks like in Schemper!



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                           FIT MODELS                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

########################### NONINTERACTIVE COX ##########################

rs.n = coxph( Surv(surv.time) ~ start.dose, data=td ); summary(rs.n)
(a = cox.zph(rs.n) )


########################### WRONG INTERACTION ##########################

rs.w = coxph( Surv(surv.time) ~ 
              start.dose + start.dose:log(surv.time), data=td ); summary(rs.w)
(a = cox.zph(rs.w) )
# violates PH horribly


########################### RIGHT INTERACTION: MINE ##########################

# THIS IS THE SAME AS TERRY'S VERSION BELOW :)

# make split data (split on event times) using my function
#split.mm.2 = event_split(d=td, id.var="id", followup.var="surv.time", event.var="event", split.at="event", 
#                       name.prefix="name.prefix", write.path=NA)$split.data
#dim(split.mm.2)

# right model with frailties
#rs.r = coxph(Surv(start.time, stop.time, event) ~ start.dose + start.dose : log(stop.time) +
#               frailty(id), data=split.mm.2 ); summary(rs.r)
#(a = cox.zph(rs.r) )
# does NOT violate PH :D


########################### RIGHT INTERACTION: TERRY THERNEAU ##########################
rs.r = coxph(Surv(surv.time) ~ start.dose + tt(start.dose), data=td,
             tt = function(x, t, ...) x * log(t) ); summary(rs.r)
(a = cox.zph(rs.r) )


########################### SPLINE ##########################
rs.spl = coxph( Surv(surv.time) ~ start.dose + tt(start.dose), data=td,
                  tt = function(x, t, ...) x * pspline(t) ); summary(rs.spl)
(a = cox.zph(rs.spl) )




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                        MODEL RESULTS                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

########################### EQUATIONS FOR EACH MODEL #########################

# QUOTE THE FUNCTION PART IF PUTTING IN DATAFRAME BELOW

# HR based on Weibull parameters
eq.true = "( (nu1 / lambda1)*(x / lambda1)^(nu1 - 1) ) / ( (nu0 / lambda0)*(x / lambda0)^(nu0 - 1) )"

eq.n = "exp( coef(rs.n)[1] )"

eq.w = "exp( coef(rs.w)[1] + coef(rs.w)[2]*log(x) )"

eq.r = "exp( coef(rs.r)[1] + coef(rs.r)[2]* log(x) )"




functionize = function(string) {
  paste("function(x) {", string, "}")
}


########################### PLOT #1 #########################


ggplot( data.frame( x=seq(0,3,.1) ), aes(x=x) ) +
  stat_function( fun=eval(parse(text=functionize(eq.true) ) ), color="black", lwd=lwd ) +
  stat_function( fun=eval(parse(text=functionize(eq.w) ) ), color="red", lwd=lwd ) +
  stat_function( fun=eval(parse(text=functionize(eq.n) ) ), color="orange", lwd=lwd ) +
  stat_function( fun=eval(parse(text=functionize(eq.r) ) ), color="green", lwd=lwd ) +
  theme_bw() +
  xlab("Time") + ylab("Estimated HR") + 
  scale_y_continuous(limits=c(0,20)) +
  scale_x_continuous( breaks=seq(0, 3, .2) ) +
  theme(axis.title = element_text(size=16) ) +
  ggtitle( paste("nu0=", nu0, ", nu1=", nu1, ", lambda0=", lambda0, ", lambda1=", lambda1, sep="") )




########################### DIFFERENCES FROM TRUE MODEL #########################

# difference from true model
funct_diff = function(string) {
  paste("function(x) { (", string, "-", eq.true, ")^2 }") 
}


# try integrating them as a means of comparison
for (i in c(eq.n, eq.r, eq.w)) {
  expr = funct_diff(i)
  print(i)
  print( sqrt( as.numeric( integrate( eval(parse(text=expr)), lower=.01, upper=3)$value ) ) )
}






###### not using!!!

########################### PUT IN DATAFRAME ###########################


# initialize dataframe for results
scenario = c("B", "C")
model = c("true", "nonint", "wrong.int", "right.int")

results = as.data.frame( expand.grid(scenario, model) ); names(results) = c("scenario", "model")
results$equation = NA
results$schoen.p = NA
results$color = NA


for (s in "B") {  # scenarios
  for (m in model) {  # models
    
    # equation to plot
    results$equation[results$scenario==s & results$model==m] = 
     switch(m, "true"=eq.true, "nonint"=eq.n, "wrong.int"=eq.w, "right.int"=eq.r)
    
    # test of PH violation
    if (m != "true") {
      cox = switch(m, "nonint"=rs.n, "wrong.int"=rs.w, "right.int"=rs.r)
      a = cox.zph(cox)[[1]]
      schoen.p = a[ dim(a)[1], dim(a)[2] ]
    }
    results$schoen.p[results$scenario==s & results$model==m] = schoen.p
    
    # color for plot
    results$color[results$scenario==s & results$model==m] = 
      switch(m, "true"="black", "nonint"="blue", "wrong.int"="red", "right.int"="green")
  }  
}






########################### PLOT FITTED HRs OVER TIME ###########################

lwd=1.5


ggplot( data.frame( x=seq(0,3,.1) ), aes(x=x) ) +
    xlab("Time") + ylab("Estimated HR") + 
    scale_y_continuous(limits=c(0,20)) +
    scale_x_continuous( breaks=seq(0, 3, .2) ) +
    theme(axis.title = element_text(size=16) )

for (i in nrow(results)) {
  eq.string = eval( parse(text=results$equation[i]) )
  string = stat_function( fun=, color=results$color[i], lwd=lwd )
  eval( parse(text = ) )

}



