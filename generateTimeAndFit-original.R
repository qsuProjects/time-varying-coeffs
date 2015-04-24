if(FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}


###########################################
#this script loads simulated covariate data,
#uses those data to generate survival times,
#merges those survival times with covariate data,
#then fits various models on results

#It is written to be sourced from a script 
##executed on a multi-core system.

#Sourcing script should initialize:

#number_cores: Integer, number of cores to use
#the_seed: Integer, the seed
#covariate_data_path: String, path and filename for covariate data
#results_write_directory: String, where to write results
#beta: Data frame, one row, each column has a corresponding column
##in the covariate file. Names column names should be identical.
##values are the values to use to generate linear predictor
##for use in survival time generation algorithm
#range_low: Integer, first rep to do
#range_high: Integer, last rep to do
#sim_results_name: String, file name for simulation results
#log_write_directory: where to write log files
#data_write_directory: where to write data files


require(msm)
require(survival)
require(plyr)


#set up parallelization
require(doSNOW)
cl<-makeCluster(number_cores)
registerDoSNOW(cl)

number_cores
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(.id, .results_write_directory, .log_write_directory, .sim_results_name, .segment_index) {
  set.seed(.segment_index*(2^.id))
  
  #load required packages
  require(msm)
  require(survival)
  require(plyr)
  require(coxme)
  require(geepack)
  require(gee)
  require(stringr)
  
  WORKER.ID <<- paste0("job_", str_pad(.segment_index, 4, side = "left", pad = "0"),
                       "_core_", str_pad(.id, 2, side = "left", pad = "0"))
  LOG.FILE.PATH <<- paste0(.log_write_directory, "/", WORKER.ID, ".txt")
  
  MODEL.OUTPUT.FILE.PATH <<- paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name)
  
  cat("coef,se,p,var,type,proportion_censored,rep\n", file = MODEL.OUTPUT.FILE.PATH, append = F)
  
}, results_write_directory, log_write_directory, sim_results_name, segment_index)


##main simulation loop
l_ply(range_low:range_high, .parallel = T, function(.rep, .covariate_data_path, .beta, .results_write_directory, 
                                                    .sim_results_name, .log_write_directory, .data_write_directory) {
  
  # .coxme_object <- coxme(.coxme_formula, data = .data_none)
  # class(bdsmatrix::diag(.coxme_object$var))
  # as.vector(.coxme_object$var)
  
  # see coxme:::print.coxme
  coxmeCoefs <- function (.coxme_object) 
  {
    .tmp <- ""
    .beta <- .coxme_object$coefficients
    .nvar <- length(.beta)
    .nfrail <- nrow(.coxme_object$var) - .nvar
    .omit <- .coxme_object$na.action
    if (.nvar > 0) {
      .se <- sqrt(bdsmatrix::diag(.coxme_object$var)[.nfrail + 1:.nvar])
      .tmp <- cbind(.beta, 
                    .se, 
                    signif(1 - pchisq((.beta/.se)^2, 1), 2))
      dimnames(.tmp) <- list(names(.beta), c("coef", 
                                             "se(coef)", 
                                             "p"))
    }
    return(.tmp)
  }
  
  .rep_str <- str_pad(.rep, 4, pad = "0")
  # .covariate_data_path <- "~/shexport/PCORI/genSurv/data/"
  # .covariate_data_path <- covariate_data_path
  # .current_data_file <- "2015-02-01_job_10_dataset_1"
  # .beta <- betas
  # .rep <- (range_low:range_high)[14]
  
  #load data
  .current_data_file <- list.files(.covariate_data_path)[.rep]
  cat(paste0("Reading covariates from ", .current_data_file, " JOB: ", WORKER.ID,
             .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  .covariates <- read.table(file = paste0(.covariate_data_path, "/", .current_data_file), 
                            sep = ",", header = T, stringsAsFactors = F)
  
  cat(paste0("Covariates read from ", .current_data_file, " JOB: ", WORKER.ID,
             .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  
  #find number of observations per person
  .n_observations_per <- unique(table(.covariates$id))
  
  ### process covariates ####
  
  #viral load
  #split into quintiles
  .covariates$log_vln_quint <- cut(.covariates$log_vln, 
                                   breaks = quantile(.covariates$log_vln, seq(from = 0, to = 1, by = 0.2)),
                                   labels = c("first", "second", "third", "fourth", "fifth"))
  .covariates$log_vln_2 <- .covariates$log_vln_3 <- .covariates$log_vln_4 <- .covariates$log_vln_5 <- 0
  .covariates$log_vln_2[.covariates$log_vln_quint == "second"] <- 1
  .covariates$log_vln_3[.covariates$log_vln_quint == "third"] <- 1
  .covariates$log_vln_4[.covariates$log_vln_quint == "fourth"] <- 1
  .covariates$log_vln_5[.covariates$log_vln_quint == "fifth"] <- 1
  
  #cd4
  #use cutpoints from paper
  .covariates$cd4_cuts <- cut(.covariates$cd4, 
                              breaks = c(-1e6, 50, 100, 200, 350, 1e6),
                              labels = c("lt_50", "50_100", "100_200", "200_350", "350+"))
  .covariates$ind_cd4_50_100 <- 
    .covariates$ind_cd4_100_200 <-
    .covariates$ind_cd4_200_350 <- 
    .covariates$ind_cd4_350_500 <- 0
  
  .covariates$ind_cd4_50_100[.covariates$cd4_cuts == "50_100"] <- 1
  .covariates$ind_cd4_100_200[.covariates$cd4_cuts == "100_200"] <- 1
  .covariates$ind_cd4_200_350[.covariates$cd4_cuts == "200_350"] <- 1
  .covariates$ind_cd4_350_500[.covariates$cd4_cuts == "350+"] <- 1
  
  #bmi
  #use cutpoints from paper
  
  #use cutpoints from paper
  .covariates$bmi_cuts <- cut(.covariates$bmi, 
                              breaks = c(-1e6, 20, 25, 30, 1e6),
                              labels = c("lt_20", "20_25", "25_30", "gt30"))
  .covariates$ind_bmi_lt_20 <- 
    .covariates$ind_bmi_25_30 <-
    .covariates$ind_bmi_gt_30 <- 0
  
  .covariates$ind_bmi_lt_20[.covariates$bmi_cuts == "lt_20"] <- 1
  .covariates$ind_bmi_25_30[.covariates$bmi_cuts == "25_30"] <- 1
  .covariates$ind_bmi_gt_30[.covariates$bmi_cuts == "gt30"] <- 1
  
  #generate dummy variables for race
  .covariates$race_black <- 0
  .covariates$race_other <- 0
  
  .covariates$race_black[.covariates$racecatexp == 1] <- 1
  .covariates$race_other[.covariates$racecatexp == 2] <- 1
  
  #make age constant for each subject
  .covariates <- ddply(.covariates, .(id), function(.df) {
    .df$age <- .df$age[1]
    return(.df)
  })
  
  #determine how many subjects are in the data
  .n_subjects <- length(unique(.covariates$id))
  
  #generate time-dependent lambda
  # .covariates <- .covariates[ , c(names(.beta), "id")]
  .covariates$linpred <- as.matrix(.covariates[ , names(.beta)]) %*% t(.beta)
  .covariates$frailty <- rep(rnorm(n = .n_subjects, mean = 0, sd = 0.5 ), each = .n_observations_per)
  
  .covariates$xB <-  exp(.covariates$linpred + .covariates$frailty)
  
  #put data frame in list form by id
  .covars_list <- dlply(.covariates, .(id))
  
  #define simulation parameters
  .ZBbar <- mean(.covariates$linpred)
  .nu <- 2
  .mediangoal <- 50
  #nu = 10 and mediangoal = 50 work
  .lambda <- (log(2)/exp(.ZBbar))*.mediangoal^(-.nu)
  
  # the g function is defined as the inverse of the baseline cummulative hazard from
  ## a Weibull with shape nu and scale lambda defined above
  .g <- function(x){
    ((1/.lambda)*x)^(1/.nu)
  }
  .g_inv <- function(x){
    .lambda*(x^.nu)
  }  
  
  .t <- 0:(.n_observations_per-1)
  .t_diff <- (.t[-1] - .t[1:(length(.t) - 1)])[-(length(.t) - 1)]
  .g_inv_t <- .g_inv(.t)
  .g_inv_t_diff <- (.g_inv(.t[-1]) - .g_inv(.t[1:(length(.t) - 1)]))[-(length(.t) - 1)]
  
  #CREATING THE BOUNDS OF TRUNCATION
  .t_max <- .n_observations_per
  .t_min <- 1
  
  .g_inv_t_max <- .g_inv(.t_max)
  .g_inv_t_min <- .g_inv(.t_min)
  
  
  #K function applies ACCEPT-REJECT algorithm
  .k <- function(..x, ..m, ..M, ..rates, ..t){
    ifelse(..x <= ..m | ..x >= ..M, 0, dpexp(..x, ..rates, ..t))
  }
  
  #define survival time generation function
  .gen_y <- function(.x, .g_inv_t,  .g_inv_t_min, .g_inv_t_max) {
    .x1 <- .x$xB
    .d <- ppexp(.g_inv_t_max, .x1, .g_inv_t) - ppexp(.g_inv_t_min, .x1, .g_inv_t)
    .M <- 1 / .d
    .r <- 60
    .count<-0
    #counter of times repeat is run
    while (.count<1000) {
      .count <- .count+1
      .y <- rpexp(.r, .x1, .g_inv_t)
      .u <- runif(.r)
      .t <- .M * (.k(.y, .g_inv_t_min, .g_inv_t_max, .x1, .g_inv_t) / .d / dpexp(.y, .x1, .g_inv_t))
      .y <- .y[.u <= .t][1]
      if (!is.na(.y)) {break}
    }
    print(.count)
    .y
  }
  
  cat(paste0("Covariates processed for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  .survival_times <- ldply(.covars_list, .gen_y, .g_inv_t,  .g_inv_t_min, .g_inv_t_max)
  cat(paste0("Survival times generated for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  .survival_times$g_y <- .g(.survival_times[ ,2])
  summary(.survival_times$g_y)
  hist(.survival_times$g_y)
  
  if (sum(is.na(.survival_times$V1)) == 0) {
    
    #create uncensored model-ready dataset
    .data_none <-  ldply(1:.n_subjects, function(..subject, ..survival_times, ..covars_list) {
      ..survival_time <- ceiling(..survival_times$g_y[..subject])
      ..to_return <- ..covars_list[[..subject]][1:..survival_time, ]
      ..to_return$id <- ..subject
      ..to_return$t <- c(1:..survival_time)
      ..to_return$t0 <- 0:(..survival_time - 1)
      ..to_return$d <- c( rep(0, ..survival_time - 1), 1)
      ..to_return$proportion_censored = 0
      
      return(..to_return)
      
    }, .survival_times, .covars_list)
    
    .data_none$source_file <- .current_data_file
    
    cat(paste0("Uncensored data set generated for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    write.table(.data_none, file = paste0(.data_write_directory, "/SURV_", .current_data_file), sep = ",", row.names = F, col.names = T)
    cat(paste0("Uncensored data written for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    #     #create randomly censored model-ready dataset
    #     .data_censored <-  llply(1:3, function(..censor, ..survival_times, ..covars_list, ..n_subjects, 
    #                                            ..data_write_directory, ..log_write_directory, ..rep) {
    #       ..proportion_censored <- c(0.2, 0.5, 0.8)[..censor]
    #       ..events <- rbinom(.n_subjects, size = 1, p = (1 - ..proportion_censored) )
    #       
    #       ..censored_set <- ldply(1:..n_subjects, function(...subject, ...survival_times, ...covars_list, ...events, ...proportion_censored) {
    #         ...survival_time <- ceiling(...survival_times$g_y[...subject])
    #         ...to_return <- ...covars_list[[...subject]][1:...survival_time, ]
    #         ...to_return$id <- ...subject
    #         ...to_return$t <- c(1:...survival_time)
    #         ...to_return$t0 <- 0:(...survival_time - 1)
    #         ...to_return$d <- ...events[...subject]
    #         ...to_return$proportion_censored = ...proportion_censored
    #         return(...to_return)
    #         
    #       }, ..survival_times, ..covars_list, ..events, ..proportion_censored)
    #       
    #       write.table(..censored_set, file = paste0(..data_write_directory, "/", WORKER.ID, "_", ..rep, "_cens", ..proportion_censored*100,".csv"),
    #                   sep = ",", row.names = F, col.names = T)
    #       cat(paste0(..proportion_censored, " censored data written for rep ", ..rep, " at ", Sys.time(), "\n"), file = paste0(..log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    #       
    #     }, .survival_times, .covars_list, .n_subjects, .data_write_directory, .log_write_directory, .rep)
    #     
    #     cat(paste0("Censored data sets generated for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    #     
    #initialize model container list
    .uncensored_results <- list()
    
    #fit uncensored models
    #     .unpooled_cox_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(names(.beta), collapse = " + ")))
    #     .unpooled_cox_results <- data.frame(summary(coxph(.unpooled_cox_formula, data = .data_none))$coef[ ,c(1,3,5)])
    #     names(.unpooled_cox_results) <- c("coef", "se", "p")
    #     .unpooled_cox_results$var <- row.names(.unpooled_cox_results)
    #     .unpooled_cox_results$type <- "Unpooled Cox"
    #     .unpooled_cox_results$proportion_censored <- 0
    #     
    #     .uncensored_results[[1]] <- .unpooled_cox_results
    
    .pooled_cox_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(names(.beta), collapse = " + "), " + cluster(id) "))
    .pooled_cox_results <- data.frame(summary(coxph(.pooled_cox_formula, data = .data_none), )$coef[ ,c(1,4,6)])
    names(.pooled_cox_results) <- c("coef", "se", "p")
    .pooled_cox_results$var <- row.names(.pooled_cox_results)
    .pooled_cox_results$type <- "Pooled Cox"
    .pooled_cox_results$proportion_censored <- 0
    
    .uncensored_results[[1]] <- .pooled_cox_results
    
    cat(paste0("Pooled Cox fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    ################################      
    .coxme_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(names(.beta), collapse = " + "), " + (1|id) "))
    .coxme_results <- data.frame(coxmeCoefs(coxme(.coxme_formula, data = .data_none)))
    names(.coxme_results) <- c("coef", "se", "p")
    .coxme_results$var <- row.names(.coxme_results)
    .coxme_results$type <- "Frailty Cox"
    .coxme_results$proportion_censored <- 0
    
    .uncensored_results[[2]] <- .coxme_results
    cat(paste0("Frailty Cox fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    ################################  
    
    .pooled_poisson_formula <- as.formula(paste0("d ~ ", paste0(names(.beta), collapse = " + ")))
    .pooled_poisson_results <- data.frame( summary(geeglm(data = .data_none, formula = .pooled_poisson_formula, id = id, 
                                                          corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)])
    names(.pooled_poisson_results) <- c("coef", "se", "p")
    .pooled_poisson_results$var <- row.names(.pooled_poisson_results)
    .pooled_poisson_results$type <- "Pooled Poisson"
    
    .uncensored_results[[3]] <- .pooled_poisson_results
    
    cat(paste0("Pooled Poisson fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    .pooled_logistic_formula <- as.formula(paste0("d ~ ", paste0(names(.beta), collapse = " + ")))
    .pooled_logistic_results <- data.frame( summary(geeglm(data = .data_none, formula = .pooled_logistic_formula, id = id, 
                                                           corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)])
    names(.pooled_logistic_results) <- c("coef", "se", "p")
    .pooled_logistic_results$var <- row.names(.pooled_logistic_results)
    .pooled_logistic_results$type <- "Pooled Logistic"
    
    .uncensored_results[[4]] <- .pooled_logistic_results
    
    cat(paste0("Pooled Logistic fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    .uncensored_results <- ldply(.uncensored_results)
    
    cat(paste0("Uncensored models fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    #     .censored_results <- ldply(.data_censored, function(..data_censored, ..unpooled_cox_formula, ..pooled_cox_formula,
    #                                                         ..coxme_formula, ..pooled_poisson_formula, ..pooled_logistic_formula) {
    #       
    #       ..proportion_censored <- ..data_censored$proportion_censored[1]
    #       
    #       ..censored_results <- list()
    #       
    #       ..unpooled_cox_results <- data.frame(summary(coxph(..unpooled_cox_formula, data = ..data_censored))$coef[ ,c(1,3,5)])
    #       names(..unpooled_cox_results) <- c("coef", "se", "p")
    #       ..unpooled_cox_results$var <- row.names(..unpooled_cox_results)
    #       ..unpooled_cox_results$type <- "Unpooled Cox"
    #       ..unpooled_cox_results$proportion_censored <- ..proportion_censored
    #       
    #       ..censored_results[[1]] <- ..unpooled_cox_results
    #       
    #       ..pooled_cox_results <- data.frame(summary(coxph(..pooled_cox_formula, data = ..data_censored), )$coef[ ,c(1,4,6)])
    #       names(..pooled_cox_results) <- c("coef", "se", "p")
    #       ..pooled_cox_results$var <- row.names(.pooled_cox_results)
    #       ..pooled_cox_results$type <- "Pooled Cox"
    #       ..pooled_cox_results$proportion_censored <- ..proportion_censored
    #       
    #       ..censored_results[[2]] <- ..pooled_cox_results
    #       
    #       ..coxme_results <- data.frame(coxmeCoefs(coxme(..coxme_formula, data = ..data_censored)))
    #       names(..coxme_results) <- c("coef", "se", "p")
    #       ..coxme_results$var <- row.names(..coxme_results)
    #       ..coxme_results$type <- "Frailty Cox"
    #       ..coxme_results$proportion_censored <- ..proportion_censored
    #       
    #       ..censored_results[[3]] <- ..coxme_results
    #       
    #       ..pooled_poisson_results <- data.frame( summary(geeglm(data = ..data_censored, formula = ..pooled_poisson_formula, id = id, 
    #                                                             corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)])
    #       names(..pooled_poisson_results) <- c("coef", "se", "p")
    #       ..pooled_poisson_results$var <- row.names(.pooled_poisson_results)
    #       ..pooled_poisson_results$type <- "Pooled Poisson"
    #       ..pooled_poisson_results$proportion_censored <- ..proportion_censored
    #       
    #       ..censored_results[[4]] <- ..pooled_poisson_results
    #       
    #       
    #       ..pooled_logistic_results <- data.frame( summary(geeglm(data = ..data_censored, formula = ..pooled_logistic_formula, id = id, 
    #                                                              corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)])
    #       names(..pooled_logistic_results) <- c("coef", "se", "p")
    #       ..pooled_logistic_results$var <- row.names(..pooled_logistic_results)
    #       ..pooled_logistic_results$type <- "Pooled Logistic"
    #       ..pooled_logistic_results$proportion_censored <- ..proportion_censored
    #       
    #       ..censored_results[[5]] <- ..pooled_logistic_results
    #       
    #       
    #       ..censored_results <- ldply(..censored_results)
    #       
    #       ldply(..censored_results)
    #     }, .unpooled_cox_formula, .pooled_cox_formula, .coxme_formula, .pooled_poisson_formula, .pooled_logistic_formula)
    #     
    #     cat(paste0("Censored models fit for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    #     
    .to_write <-  .uncensored_results
    # rbind(.uncensored_results, .censored_results)
    .to_write$rep <- .rep
    
    write.table(.to_write, file = MODEL.OUTPUT.FILE.PATH, 
                append = T, row.names = F, col.names = F, sep = ",")
    # paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name), append = F)
    cat(paste0("Results written for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    cat(paste0("-------------------END-REP------------------------\n"), file = LOG.FILE.PATH, append = T)
    
  } else {
    cat(paste0("Convergence fail for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  }
  gc()
}, covariate_data_path, betas, results_write_directory, sim_results_name, log_write_directory, data_write_directory)
