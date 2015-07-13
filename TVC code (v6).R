
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                              SIMULATE DATA                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(ggplot2)
library(survival)
#library(pspline)


########################### FUNCTION: SIMULATE DATA ########################### 

# n: number of subjects
# doses: vector of possible treatment doses administered at beginning of study
# doses.p: vector of probabilities of each dose
# scenario: "B" or "C"

sim_one_dataset = function(scenario, n, doses, doses.p) {
  # simulate X
  d = data.frame(matrix(ncol = 3, nrow = n))
  names(d) = c("id", "start.dose", "surv.time")
  d$id = 1:n
  d$start.dose = rep( sample(doses, size=n, replace=TRUE, prob=doses.p), each=max.t )
  
  # simulate Y from Weibull
  # use Schemper's scenarios
  if (scenario=="C") {
    lambda0 <<- 2
    nu0 <<- 1
    nu1 <<- 2
    lambda1 <<- 1
  }
  
  if (scenario=="B") {
    lambda0 <<- 1.5
    nu0 <<- 3
    nu1 <<- 1
    lambda1 <<- 1
  }
  
  # simulate from Weibull
  d$surv.time[d$start.dose==0] = rweibull(n=length(d$start.dose[d$start.dose==0]), shape=nu0, scale=lambda0) 
  d$surv.time[d$start.dose==1] = rweibull(n=length(d$start.dose[d$start.dose==1]), shape=nu1, scale=lambda1) 
  
  # everyone has event
  d$event=1
  
  return(d)
}


########################### FUNCTION: FIT MODELS ###########################

fit_models = function(d, lower.bound, upper.bound, scenario) {
  
  #browser()
  # noninteractive
  # use super-assignment so that results will be available for plotting afterward
  rs.n <<- coxph( Surv(surv.time) ~ start.dose, data=d )

  # wrong interaction
  rs.w <<- coxph( Surv(surv.time) ~ 
                  start.dose + start.dose:log(surv.time), data=d )
  
  # right interaction
  rs.r <<- coxph(Surv(surv.time) ~ start.dose + tt(start.dose), data=d,
               tt = function(x, t, ...) x * log(t) )
  
  # spline interaction
  rs.s <<- coxph(Surv(surv.time) ~ start.dose + tt(start.dose), data=d,
                 tt = function(x, t, ...) x * pspline(t) )

  # integrate difference from true model
  new.row = scenario  # for starters, the row just says what scenario it is
  for (i in c(str.n, str.r, str.w)) {
    expr = funct_diff(i)
    new.row = c( new.row, 
                 ( sqrt( as.numeric( integrate( eval(parse(text=expr)),
                                                lower=lower.bound, upper=upper.bound)$value ) ) ) )
  }
  return(new.row)
}


# difference from true model (the thing to integrate)
funct_diff = function(string) {
  paste("function(x) { (", string, "-", str.true, ")^2 }") 
}




########################### DECIDE INTEGRATION BOUNDS ###########################

#theoretical 10 and 90 percentiles
low.perc = 0.10
hi.perc = 0.90

# scenario B
lambda0 = 1.5
nu0 = 3
nu1 = 1
lambda1 = 1
grp0 = qweibull( c(low.perc, hi.perc), shape=nu0, scale=lambda0)
grp1 = qweibull( c(low.perc, hi.perc), shape=nu1, scale=lambda1)
bounds.B = colMeans( rbind(grp0, grp1) )  # average of 5 and 95 percentiles across the 2 subject groups


# scenario C
lambda0 = 2
nu0 = 1
nu1 = 2
lambda1 = 1
grp0 = qweibull( c(low.perc, hi.perc), shape=nu0, scale=lambda0)
grp1 = qweibull( c(low.perc, hi.perc), shape=nu1, scale=lambda1)
bounds.C = colMeans( rbind(grp0, grp1) )  # average of 5 and 95 percentiles across the 2 subject groups



########################### RUN SIMULATIONS ###########################

n = 500
reps = 100
doses = c(0, 1)
doses.p = rep( 1/length(doses), length(doses) )
scenarios = c("B", "C")  # which Schemper scenario to use (C = linearly diverging hazards; B=nonlinearly diverging)

# equations for each model
str.true = "( (nu1 / lambda1)*(x / lambda1)^(nu1 - 1) ) / ( (nu0 / lambda0)*(x / lambda0)^(nu0 - 1) )"
str.n = "exp( coef(rs.n)[1] )"
str.w = "exp( coef(rs.w)[1] + coef(rs.w)[2]*log(x) )"
str.r = "exp( coef(rs.r)[1] + coef(rs.r)[2]* log(x) )"


##### Scenario B #####
B = as.data.frame( matrix(nrow=1, ncol=4) )
names(B) = c("scenario", "int.n", "int.r", "int.w")

for (i in 1:reps) {
  d = sim_one_dataset(scenario="B", n, doses, doses.p)
  B = rbind(B, fit_models(d, lower.bound=bounds.B[1], upper.bound=bounds.B[2], scenario="B") )
}

B = B[-1,]  # remove dumb NA row
B[,2:4] = apply(B[,2:4], 2, as.numeric)

mean = apply(B[,2:4], 2, mean)
lo = apply(B[,2:4], 2, function(x) quantile(x, probs=.025))  # 2.5 percentile
hi = apply(B[,2:4], 2, function(x) quantile(x, probs=.975))  # 97.5 percentile 

rs = as.data.frame( cbind( mean, lo, hi ) )


##### Scenario C #####
C = as.data.frame( matrix(nrow=1, ncol=4) )
names(C) = c("scenario", "int.n", "int.r", "int.w")

for (i in 1:reps) {
  d = sim_one_dataset(scenario="C", n, doses, doses.p)
  C = rbind(C, fit_models(d, lower.bound=bounds.C[1], upper.bound=bounds.C[2], scenario="C") )
}

C = C[-1,]  # remove dumb NA row
C[,2:4] = apply(C[,2:4], 2, as.numeric)

mean = apply(C[,2:4], 2, mean)
lo = apply(C[,2:4], 2, function(x) quantile(x, probs=.025))  # 2.5 percentile
hi = apply(C[,2:4], 2, function(x) quantile(x, probs=.975))  # 97.5 percentile 

rs = rbind( rs, cbind(mean, lo, hi) )
rs = as.data.frame(rs)
rs$model = rep( c("NIM", "RM", "WM"), 2)
rs$scenario = rep( c("Decreasing HR", "Increasing HR"), each=3)


# write dataset
setwd("~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/Sim results")
write.csv(rs, paste(Sys.Date(), "tvc_sim_results.csv", sep="_") )



########################### PLOT MODEL RESULTS ###########################

pd <- position_dodge(0.4)
lwd=1.2

ggplot(rs, aes(x=model, y=mean, group=scenario, color=scenario, shape=scenario, fill=scenario ) ) +
  geom_point(size=7, position=pd) +
  geom_errorbar( aes(ymin=lo, ymax=hi, width=0.1), position=pd, lwd=1.5 ) +
  theme_bw() +
  scale_shape_manual(values=c(25,24), name="Scenario") +
  scale_color_manual(values=c("black", "orange"), name="Scenario") +
  scale_fill_manual(values=c("black", "orange"), name="Scenario") +
  xlab("Model") + ylab("Distance from truth") +
  scale_y_continuous(breaks=seq(0,60,10)) +
  theme(axis.title = element_text(size=22), axis.text = element_text(size=20), 
        legend.text = element_text(size=18), legend.title = element_text(size=18) )





########################### LOOK DESCRIPTIVELY AT SIMULATED DATA ###########################


# do just a single dataset
d = sim_one_dataset("B", n=500, doses, doses.p)
fit_models(d, lower.bound=0.01, upper.bound=3, "C")


##### Look at Distribution of Survival Times #####
ggplot(data=d, aes(x=surv.time)) + geom_density(fill="grey") + theme_bw() +
  ggtitle("") + 
  scale_x_continuous(limits=c(0,6))

# look at median survival times by dose
aggregate(surv.time ~ start.dose, d, summary )

# does it match theoretical medians?
lambda0*log(2)^(1/nu0)  # group 0
lambda1*log(2)^(1/nu1)  # group 1
# yes


###### KM Plot #####
m1 = survfit( Surv(surv.time) ~ as.factor(start.dose), data=d)
colors = c("red", "#2ca25f")
plot(m1, xlim=c(0,3), ylim=c(0,1), conf.int=FALSE, mark.time=FALSE, col=colors, xlab="Time",
     ylab="Cumulative survival probability", lwd=3, main="")
legend("topright", c("Dose=0", "Dose=1"), lty=1, lwd=3, col=colors)


##### Plot Hazards #####
library(muhaz)
haz0 = with( d[d$start.dose==0,], muhaz(surv.time))
haz1 = with( d[d$start.dose==1,], muhaz(surv.time))

plot(haz0, xlim=c(0,3), ylim=c(0,6), lwd=2, col=colors[2], xlab="Time")
lines(haz1, lwd=2, col=colors[1])

legend("topleft", legend=c("Dose=0", "Dose=1"), 
       lty=1, col=colors, lwd=2)

# yes, looks like in Schemper!






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                        PRETTY PLOT                                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# experiment with plotting
C = sim_one_dataset(scenario="C", n, doses, doses.p)
fit_models(C, lower.bound=bounds.C[1], upper.bound=bounds.C[2], scenario="C")
#mm.split.C = event_split(C, id.var="id", followup.var="surv.time", event.var="event", split.at="event", na.rm=TRUE,
#                         name.prefix="", write.path=NA)$split.data
# C is for scenario C



########################### PLOT #1 #########################


# ADD THE MARGINAL DISTRIBUTION OF SURVIVAL TIMES UNDER X-AXIS

# equations for HR as function of time for each model
eq.true = function(x) { ( (nu1 / lambda1)*(x / lambda1)^(nu1 - 1) ) / ( (nu0 / lambda0)*(x / lambda0)^(nu0 - 1) ) }
eq.n = function(x) { exp( coef(rs.n)[1] ) }
eq.w = function(x) { exp( coef(rs.w)[1] + coef(rs.w)[2]*log(x) ) }
eq.r = function(x) { exp( coef(rs.r)[1] + coef(rs.r)[2]* log(x) ) }


lwd=1.5

# temporary variable of cumulative survival for use in rug plot
num.alive = function(x) length( d$surv.time[d$surv.time >= x] )
v = Vectorize( FUN=num.alive )
S = data.frame(time=d$surv.time, num.alive=v(d$surv.time))
S = S[order(S$num.alive, decreasing=TRUE),]

d$cum.surv = function(x) {length(d$surv.time[d$surv.time<])

                          
                          
# AWESOME PLOT INCLUDING SPLINE!!

output <- data.frame(DxT = seq(min(d$start.dose) * min(d$surv.time),
                               max(d$start.dose * d$surv.time),
                               0.01))

output = data.frame(DxT = seq(1 * min(d$surv.time),
                              1 * max(d$surv.time),
                              0.01))
output$HR <- exp(rs.s$coefficients[1] + pspline(output$DxT) %*% rs.s$coefficients[-1])

to.plot = output[output$DxT <= 3,]

main = ggplot( to.plot, aes(x=DxT, y=HR) ) +
  geom_line(color="blue", lwd=lwd) +
  stat_function( fun=eq.true, color="black", lwd=lwd ) +
  stat_function( fun=eq.w, color="red", lwd=lwd ) +
  stat_function( fun=eq.n, color="orange", lwd=lwd ) +
  stat_function( fun=eq.r, color="green", lwd=lwd ) +
  theme_bw() +
  xlab("Time") + ylab("Estimated HR (Dose 1 vs. 0)") + 
  scale_y_continuous(limits=c(0,20)) +
  scale_x_continuous( breaks=seq(0, 3, .2) ) +
  theme(axis.title = element_text(size=16) ) +
  
  # "rug" from generated dataset to show surv time density
  # works but isn't quite what I wanted
  #geom_rug( data=d[d$surv.time<3,], mapping=aes(x=num.alive, y=start.dose),
   #         alpha=0.4, sides="b") +

  # try for cumulative survival instead
  # kind works! 
  # BOOKMARK: need to make it like a 1D heatmap :)
  # look into geom_bin2d
  #geom_step( data=d[d$surv.time<3,], aes(x=surv.time, y=..y..), stat="ecdf" ) +
  
  # WORK? I THINK?
  geom_rug( data=d[d$surv.time<3,], aes(x=surv.time, y=..y..), stat="ecdf",
                     alpha=0.4, sides="b") +
  
  ggtitle( paste("nu0=", nu0, ", nu1=", nu1, ", lambda0=", lambda0, ", lambda1=", lambda1, sep="") )
plot(main)






# this one works
rug = ggplot( d, aes(x=surv.time, y=start.dose)) + geom_point() +
  scale_x_continuous(limits=c(0,3)) + geom_rug(col=rgb(.5,0,0,alpha=.2))
plot(rug)

# try to use empty ggplot
# http://stackoverflow.com/questions/31069324/adding-points-from-other-dataset-to-ggplot2
rug = ggplot() +
  scale_x_continuous(limits=c(0,3)) + geom_rug( data=d, mapping=aes(x=surv.time, y=start.dose),
                                                col=rgb(.5,0,0,alpha=.3))
plot(rug)


# marginal density of survival times for below X-axis
# try using this answer: http://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
library(gridExtra)
hist = ggplot( d, aes(x=surv.time)) + geom_density(fill="black") + theme_bw() +
  scale_x_continuous(limits=c(0,3))
plot(hist)

grid.arrange(main, hist, ncol=1, nrow=2, widths=c(4, 1), heights=c(1, 4))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                        OTHER THINGS                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

########################### TIMECOX PLOT ###########################

library(timereg)
tc = timecox( Surv(surv.time) ~ start.dose, data=d )
# looks VERY close to the truth.

plot( cox.zph(rs.n) )




