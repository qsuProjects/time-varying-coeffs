if(FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}

require(msm)
require(survival)
require(plyr)

setwd("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/")

scenario <- "ALLCHECK"
#load data
peeps <- read.table(file = "2014-10-16_covariates_for_kk.csv", sep = ",", header = T, stringsAsFactors = F)
peeps[1:4,]

peeps$race_black <- 0
peeps$race_other <- 0

peeps$race_black[peeps$racecatexp == 1] <- 1
peeps$race_other[peeps$racecatexp == 2] <- 1

peeps$log_cd4 <- log(peeps$cd4)
# peeps$log_vln <- log(peeps$vln)

n_subjects <- length(unique(peeps$id))

# summary(peeps)
##define beta
betas <- data.frame(male = 0.15,
                    d_abac = 0.7, 
                    d_ataz = 0.5,
                    d_dida = 0.3,
                    d_efav = 0,
                    d_emtr = 0,
                    #                     d_indi = 0,
                    d_lami = 0.5,
                    d_lopi = 0,
                    d_nelf = 0,
                    #                     d_nevi = 0,
                    d_rito = 0,
                    #                     d_saqu = 0,
                    d_stav = 0,
                    d_teno = 0,
                    d_zido = 0,
                    age = 0.02,
                    bmi = 0.05,
                    log_cd4 = 0.09,
                    log_vln = 0.09,
                    bps = 0.005,
                    bpd = 0.005,
                    ldl = 0.005,
                    hdl = 0.002,
                    trig = 0.0002,
                    race_black = -0.2,
                    race_other = 0
)


# covars <- peeps[ , c(names(.beta), "id", "time.point")]
# covars$linpred <- (as.matrix(covars[ ,names(.beta)]) %*% t(.beta))
# #covars$frailty <- rep(rgamma(n = .n_subjects, shape = 20, rate = 20), each = 350)
# covars$frailty <- rep(rnorm(n = .n_subjects, mean = 0, sd = 3 ), each = 350)
# covars$xB <-  exp(covars$linpred + covars$frailty)
# 
# covars_list <- dlply(covars, .(id))

n_clusters <- 4
library(plyr)
require(doSNOW)
cl<-makeCluster(n_clusters)
registerDoSNOW(cl)

n_clusters
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(id, scenario) {
  WORKER.ID <<- paste0("core_", id)
  set.seed(10*id)
  cat("coef,se,p,var,type,rep\n", file = paste0("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/generatedData/", 
                                                WORKER.ID, "_", scenario, "_results.csv"), append = F)
  
  
}, scenario)
# .rep <- 6
# .all_data <- peeps
# .n_subjects <- n_subjects
# .beta <- betas
l_ply(1:1000, .parallel = T, function(.rep, .all_data, .n_subjects, .beta, .scenario) {
  require(msm)
  require(survival)
  require(plyr)
  require(coxme)
  require(geepack)
  require(gee)
  source("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/coxmeCoefs.R")
  #   WORKER.ID <- paste0("core_", .rep)
  print(.rep)
  
  
  #   bootstrap_sample <- sample(1:length(unique(.all_data$id)), size = .n_subjects, replace = T)
  #   bootstrap_sample[1:10]
  #   bootstrap_indices <- ldply(bootstrap_sample, function(..id) {
  #   #      print(..id)
  #     return( data.frame(index = ((..id-1)*(350) + 1):((..id)*(350))) )
  #   })
  #   bootstrap_df <- .all_data[bootstrap_indices$index,]
  bootstrap_df <- .all_data
  bootstrap_df$bootstrap_id <- bootstrap_df$id
  #   bootstrap_df$bootstrap_id <- rep(1:.n_subjects, each = 350)
  
  #   covars[1:10,]
  #generate time-dependent lambda
  covars <- bootstrap_df[ , c(names(.beta), "id", "bootstrap_id", "time.point")]
  covars$linpred <- (as.matrix(covars[ ,names(.beta)]) %*% t(.beta))
  #covars$frailty <- rep(rgamma(n = .n_subjects, shape = 20, rate = 20), each = 350)
  covars$frailty <- rep(rnorm(n = .n_subjects, mean = 0, sd = 0.5 ), each = 350)
  covars$xB <-  exp(covars$linpred + covars$frailty)
  
  #   ddply(covars, .(bootstrap_id), function(.df) {return(length(unique(.df$frailty)))})
  summary(covars$frailty)
  summary(covars$linpred)
  #   covars[1,]
  #put data frame in list form by id
  covars_list <- dlply(covars, .(bootstrap_id))
  
  #define simulation parameters
  ZBbar <- mean(covars$linpred + covars$frailty)
  nu <- 20
  mediangoal <- 75
  
  median <- mediangoal
  lambda <- (log(2)/exp(ZBbar))*median^(-nu)
  
  # the g function is defined as the inverse of the baseline cummulative hazard from
  ## a Weibull with shape nu and scale lambda defined above
  g <- function(x){
    ((1/lambda)*x)^(1/nu)
  }
  g.inv <- function(x){
    lambda*(x^nu)
  }  
  
  t <- 0:349
  t.diff <- (t[-1] - t[1:(length(t) - 1)])[-(length(t) - 1)]
  g.inv.t <- g.inv(t)
  g.inv.t.diff <- (g.inv(t[-1]) - g.inv(t[1:(length(t) - 1)]))[-(length(t) - 1)]
  
  #CREATING THE BOUNDS OF TRUNCATION
  t.max <- 350
  t.min <- 1
  
  g.inv.t.max <- g.inv(t.max)
  g.inv.t.min <- g.inv(t.min)
  
  
  #K function applies ACCEPT-REJECT algorithm
  k <- function(x, m, M, rates, t){
    ifelse(x <= m | x >= M, 0, dpexp(x, rates, t))
  }
  
  #define survival time generation function
  gen.y <- function(.x, .g_inv_t,  .g_inv_t_min, .g_inv_t_max) {
    .x1 <- .x$xB
    #     print(.x$bootstrap_id[1])
    #       print(length(.x$xB))
    #   print(length(.g_inv_t))
    .d <- ppexp(.g_inv_t_max, .x1, .g_inv_t) - ppexp(.g_inv_t_min, .x1, .g_inv_t)
    .M <- 1 / .d
    .r <- 60
    .count<-0
    #counter of times repeat is run
    while (.count<1000) {
      .count <- .count+1
      .y <- rpexp(.r, .x1, .g_inv_t)
      .u <- runif(.r)
      .t <- .M * (k(.y, .g_inv_t_min, .g_inv_t_max, .x1, .g_inv_t) / .d / dpexp(.y, .x1, .g_inv_t))
      .y <- .y[.u <= .t][1]
      if (!is.na(.y)) {break}
    }
    .y
  }
  
  y <- ldply(covars_list, gen.y, g.inv.t,  g.inv.t.min, g.inv.t.max)
  sum(is.na(y$V1))
  y$g.y <- g(y[ ,2])
  #   summary( y$g.y)
  #   hist( y$g.y)
  
  if (sum(is.na(y$V1)) == 0) {
    
    #CREATING DATASET
    data_none <- list()
    # i <- 1
    for (i in 1:.n_subjects) {
      #       print(i)
      id.temp <- rep(i, ceiling(y$g.y[i]))
      time.temp <- c(1:ceiling(y$g.y[i]))
      time0.temp <- 0:ceiling(y$g.y[i] - 1)
      d.temp <- c(rep(0, length(time.temp) - 1), 1)
      z.temp <- covars_list[[i]][1:ceiling(y$g.y[i]), c(names(.beta), "frailty", "time.point")]
      data.temp <- z.temp
      data.temp$id <- id.temp
      data.temp$t <- time.temp
      data.temp$t0 <- time0.temp
      data.temp$d_none <- d.temp 
      data_none[[i]] <- data.temp
    }
    #     peeps[peeps$id == 95,][1:10,]
    # data_none[[1]]
    data_none <- ldply(data_none)
    
    if (.rep <= 4) {
      data_write <- data_none
      data_write$rep <- .rep
      write.table(data_write, file = paste0("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/generatedData/", 
                                            WORKER.ID, "_", .scenario, "_generatedData.csv"), 
                  append = F, row.names = F, col.names = F, sep = ",")
    }
    #fit models
    
    cox_formula <- as.formula(paste0("Surv(t0, t, d_none) ~ ", paste0(names(.beta), collapse = " + ")))
    results <- data.frame((summary(coxph(cox_formula, data = data_none))$coef[ ,c(1,3,5)]))
    names(results) <- c("coef", "se", "p")
    results$var <- names(.beta)
    results$type <- "Unpooled Cox"
    results$rep <- .rep
    
    pooled_cox_formula <- as.formula(paste0("Surv(t0, t, d_none) ~ ", paste0(names(.beta), collapse = " + "), " + cluster(id) "))
    pooled_results <- data.frame(summary(coxph(pooled_cox_formula, data = data_none), )$coef[ ,c(1,4,6)])
    names(pooled_results) <- c("coef", "se", "p")
    pooled_results$var <- names(.beta)
    pooled_results$type <- "Pooled Cox"
    pooled_results$rep <- .rep
#     summary(data_none$d_saqu)
#     summary(data_none$d_indi)
#     [c(1,7:26)]
    coxme_formula <- as.formula(paste0("Surv(t0, t, d_none) ~ ", paste0(names(.beta), collapse = " + "), " + (1|id) "))
    coxme_results <- data.frame(coxmeCoefs(coxme(coxme_formula, data = data_none)))
    #     names(data_none)
#     names(beta)[11]
    names(coxme_results) <- c("coef", "se", "p")
    coxme_results$var <- row.names(coxme_results)
    coxme_results$type <- "Frailty Cox"
    coxme_results$rep <- .rep
#     summary(data_none)
    
    pooled_poisson_formula <- as.formula(paste0("d_none ~ ", paste0(names(.beta), collapse = " + ")))
    pooled_poisson_results <- data.frame( summary(geeglm(data = data_none, formula = pooled_poisson_formula, id = id, 
                                                         corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)])
    names(pooled_poisson_results) <- c("coef", "se", "p")
    pooled_poisson_results$var <- row.names(pooled_poisson_results)
    pooled_poisson_results$type <- "Pooled Poisson"
    pooled_poisson_results$rep <- .rep
    
    pooled_logistic_formula <- as.formula(paste0("d_none ~ ", paste0(names(.beta), collapse = " + ")))
    pooled_logistic_results <- data.frame( summary(geeglm(data = data_none, formula = pooled_logistic_formula, id = id, 
                                                          corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)])
    names(pooled_logistic_results) <- c("coef", "se", "p")
    pooled_logistic_results$var <- row.names(pooled_logistic_results)
    pooled_logistic_results$type <- "Pooled Logistic"
    pooled_logistic_results$rep <- .rep
    
    #     data_none[1:3,]
    to_write <- rbind(results, pooled_results, coxme_results, pooled_poisson_results, pooled_logistic_results)
    write.table(to_write, file = paste0("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/generatedData/", 
                                        WORKER.ID, "_", .scenario, "_results.csv"), 
                append = T, row.names = F, col.names = F, sep = ",")
  } else {
    print("Failed to Converge")
  }
  gc()
}, peeps, n_subjects, betas, scenario)

#setwd("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/")
all_results[1,]
all_results <- openFilesInDirectory("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/generatedData/", 
                                    match_string = paste0(scenario, "_results.csv"), header = F, delim_str = ",", merge = T)
names(all_results) <- c("coef","se","p","var","type","rep", "loaded_file_name")
all_results$coef <- as.numeric(all_results$coef)
all_results$se <- as.numeric(all_results$se)
all_results$p <- as.numeric(all_results$p)
# all_results$var <- "bps"
table(all_results$type)

all_results[1,]
length(unique(all_results$rep))

betacol <- data.frame(beta = t(betas))
betacol$var <- row.names(betacol)
all_results <- merge(x = all_results, y = betacol, by.x = "var", by.y = "var", all.x = T)
all_results[1,]


# .df <- all_results[1:10,]
summarized_results <- ddply(all_results, .(var, type), function(.df) {
  
  return(data.frame(true = .df$beta[1],
                    mean_estimate = mean( .df$coef, na.rm = T ),
                    mean_se = mean( .df$se, na.rm = T ),
                    n_converge = sum( !is.na( .df$coef ) ),
                    bias = mean( .df$coef, na.rm = T ) - .df$beta[1],
                    std_bias = 100 *( mean( .df$coef, na.rm = T ) - .df$beta[1] ) / mean( .df$se, na.rm = T ),
                    mse = mean( ( .df$coef - .df$beta ) ^ 2, na.rm = T),
                    rmse = sqrt(mean( ( .df$coef - .df$beta ) ^ 2, na.rm = T)),
                    reject_prob = mean( .df$p < 0.05, na.rm = T ),
                    stringsAsFactors = F )
  )})

summarized_results$var_order <- paste0(summarized_results$var, " (beta = ", summarized_results$true, ")")
summarized_results$var_order <- factor(summarized_results$var_order, levels = unique(summarized_results$var_order[order(summarized_results$std_bias)]))
summarized_results[1,]
summarized_results$type_plot <- gsub(" ", "\n", summarized_results$type)
summarized_results$type_plot <- factor(summarized_results$type_plot, levels = unique(summarized_results$type_plot)[c(5,2,1,3,4)])
summarized_results[ ,c("var", "true", "mean_estimate", "mean_se", "bias", "reject_prob", "std_bias")]

pdf(paste0("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/plot/model_results_", scenario, ".pdf"), 
    height = 18, width = 6)

ggplot(data = summarized_results, aes(x = std_bias, y = var_order)) +
  geom_point() +
  facet_grid( type_plot ~ .) +
  geom_vline(xintercept = c(-40, 40), color = "red") +
  scale_x_continuous("Standardized Bias", limits = c(-100,100)) +
  scale_y_discrete("Variable")

ggplot(data = summarized_results, aes(x = rmse, y = var_order)) +
  geom_point() +
  facet_grid( type_plot ~ .) +
  scale_x_continuous("RMSE") +
  scale_y_discrete("Variable")

ggplot(data = summarized_results[summarized_results$mean_se < 30, ], aes(x = mean_se, y = var_order)) +
  geom_point() +
  facet_grid( type_plot ~ .) +
  scale_x_continuous("Mean Model Standard Error") +
  scale_y_discrete("Variable")

ggplot(data = summarized_results, aes(x = bias, y = var_order)) +
  geom_point() +
  facet_grid( type_plot ~ .) +
  scale_x_continuous("Bias") +
  scale_y_discrete("Variable")


ggplot(data = summarized_results, aes(x = reject_prob, y = var_order)) +
  geom_point() +
  geom_vline(xintercept = c(0.05, 0.8), color = "red") +
  facet_grid( type_plot ~ .) +
  scale_x_continuous("Rejection Probability") +
  scale_y_discrete("Variable")

dev.off()


survival_times <- openFilesInDirectory("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/generatedData/", 
                                       match_string = paste0(scenario, "_generatedData.csv"), 
                                       header = F, delim_str = ",", merge = T)
data_none[1,]

survival_times[1:1,]
survival_times <- survival_times[ ,27:31]
names(survival_times) <- c("id", "t", "t0", "d", "rep")
names(survival_times) <- c(names(betas), "frailty", "row", "id", "t", "t0", "d", "rep", "source_file")

library(dplyr)
survival_times_tow <- select(survival_times, -row)
survival_times_tow <- select(survival_times_tow, -source_file)
survival_times_tow <- filter(survival_times_tow, rep == 1)
survival_times_tow <- select(survival_times_tow, -rep)
survival_times_tow[1,]
write.table(survival_times_tow, file = "generatedData/covariatesSurvivalTimes.csv",
            sep = ",", row.names = F, col.names = T)

survival_times[survival_times$rep == 1, ]
survival_times <- survival_times[survival_times$d == 1, ]


survival_times <- ddply(survival_times, .(rep), function(.df) {
  .df$median = median(.df$t)
  return(.df)
})

pdf(paste0("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/plot/survival_times", scenario, ".pdf"))

ggplot(data = survival_times, aes(x = t, fill = factor(rep))) + 
  geom_density(alpha = 0.25) +
  geom_vline(xintercept = unique(survival_times$median)) +
  scale_x_continuous("Survival Time", breaks = c(75, seq(0,350, 50)), limits = c(0, 350)) +
  scale_y_continuous("Density") +
  scale_fill_brewer("Replication\nNumber", palette = "Set1")


ggplot(data = survival_times, aes(y = t, x = factor(rep))) + 
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete("Replication\nNumber") +
  scale_y_continuous("Survival Time", breaks = c(75, seq(0,350, 50)), limits = c(0, 350)) 

dev.off()
