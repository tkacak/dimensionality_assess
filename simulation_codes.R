#### Data Generation Code
#Author: Tugay Kacak

##---Load required library and functions##
library(mvtnorm) 
library(dplyr)
library(tidyr)
load_mat <- function(k, p, vpf, pf, sf,vpf_unique = TRUE){
  
  switch(pf,
         "low" = {
           pf_low = 0.35 
           pf_upper = 0.5
         },
         "medium" = {
           pf_low = 0.5
           pf_upper = 0.65
         },
         "high" = {
           pf_low = 0.65
           pf_upper = 0.8
         })
  switch(sf,
         "none" = {
           sf_low = 0 
           sf_upper = 0
         },
         "low" = {
           sf_low = 0
           sf_upper = 0.1
         },
         "medium" = {
           sf_low = 0.1
           sf_upper = 0.2
         })
  
  x <- runif(p, pf_low, pf_upper) 
  y <- runif(p*(k-1) , sf_low, sf_upper)
  
  i <- 1:(p)
  if(vpf_unique){
    j <- rep(1:k, each=vpf) 
  }
  else{
    j <- rep(1:k, times=vpf) 
  }
  
  L <- matrix(NA, p, k) 
  L[cbind(i, j)] <- x 
  L[is.na(L)] <- y
  L
}

#Generate ordinal datasets with normal and different two levels skewed distributions
ordinal_simulation <- function( N,vpf, k, rho, pf, sf, distr, cat){
  # simulate data
  # simulate loading pattern
  p = vpf*k
  L <- load_mat(k, p, vpf, pf, sf)
  fcor <- matrix(rho, k,k) + diag(1-rho,k,k)
  Sigma <- L%*%fcor%*%t(L) + diag(diag(diag(p)-L%*%fcor%*%t(L)))
  
  # simulate data with specific distribution depending on simulation condition
  
  dat <- switch(distr,
                "ordinal" = {
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      qnorm((1:(cat-1))/cat),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat
                },
                "moderately" = {
                  
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      sort(sample(c(-1,1),1)*qnorm((1:(cat-1))/(1.5*cat))),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat
                },
                "severely" = {
                  
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      sort(sample(c(-1,1),1)*qnorm((1:(cat-1))/(2.5*cat))),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat}
  )#go switch function
  list(dat = dat, L = L, p = p, k = k, N = N, distr = distr)
}

#Fully-crossed design
cond <- expand.grid(cat = c(2,5),N = c(250, 1000), vpf = c(4,7), k = c(1,2,4), 
                    rho = c(0,0.45), pf = c("low", "medium","high"), sf = c("none","medium"), distr = c("ordinal","moderately","severely"))
#Removing nonsense conditions
simulation_design <- list(data_simulation = 
                            cond[!(cond$pf == "high" & cond$sf == "medium") & !(cond$k == 1 & cond$sf == "medium") & !(cond$k == 1 & cond$rho > 0),])
#Conditions list
current_condition <- simulation_design$data_simulation

#Directory for saving datasets
output_folder <- "D:/Ranalysis/makale/dimensionality_assessment_of_ordinal_items_with_common_factor_retention_criteria/datasets"
analyzes_folder <- "D:/Ranalysis/makale/dimensionality_assessment_of_ordinal_items_with_common_factor_retention_criteria/analyzes"

n_conditions=nrow(current_condition)
n_rep=500

#Generate Datasets and Save Them to "output_folder"
for (i in 1:n_conditions) {
  set.seed(1234+i)
  # Generate the replications for the current condition
  sim_data_i <- lapply(1:n_rep, function(j) {
    set.seed(1234+i*1000+j)
    ordinal_simulation(
      cat = current_condition[i, 1],
      N = current_condition[i, 2],
      vpf = current_condition[i, 3],
      k = current_condition[i, 4],
      rho = current_condition[i, 5],
      pf = current_condition[i, 6],
      sf = current_condition[i, 7],
      distr = current_condition[i, 8]
    )$dat
  })
  
  # Construct the filename dynamically based on the condition number
  filename <- paste0(output_folder, "/sim_data_condition_", i, ".RDS")
  
  # Save the dataset for the current condition
  saveRDS(sim_data_i, file = filename)
  
  # Remove the dataset from memory to free up space
  rm(sim_data_i)
}

# Source Empirical Kaiser Criterion
source("D:/Ranalysis/makale/dimensionality_assessment_of_ordinal_items_with_common_factor_retention_criteria/empkc.R")

#Multicore Analyzes
library(doParallel)
library(foreach)
#select cores
ncores = detectCores()
n.cluster = makeCluster(ncores - 4)
registerDoParallel(cl=n.cluster)

#Estimations with Parallel Programming
est <- foreach(i = 1:n_conditions) %dopar% {
  # Load the dataset for the current condition
  sim_data <- readRDS(paste0(output_folder, "/sim_data_condition_", i, ".RDS"))
  ekc <- lapply(1:n_rep, function(j) EKC(data = sim_data[[j]],restricted = T))
  pa_pca <- lapply(1:n_rep, function(j) EFAtools::PARALLEL(x = sim_data[[j]],n_datasets = 500,percent = 95,eigen_type = "EFA",decision_rule = "percentile")$n_fac_EFA)
  kgc <- lapply(1:n_rep,function(j) (EFAtools::KGC(x = sim_data[[j]],eigen_type = "EFA",method = "PAF")$n_fac_EFA))
  cd  <- lapply(1:n_rep,function(j) (EFAtools::CD(x = sim_data[[j]],N_samples = 100,N_pop = 10000,alpha = 0.3,cor_method = "pearson")$n_factors))
  condition_results <- 
       list(
       pa_pca = pa_pca,
       ekc = ekc,
       cd = cd,
       kgc = kgc)
  filename <- paste0(analyzes_folder, "/analysis_condition_", i, ".RDS")
  saveRDS(condition_results, file = filename)
}
stopCluster(n.cluster)

est <- lapply(1:n_conditions, function(i)
  readRDS(paste0(analyzes_folder,"/analysis_condition_",i,".RDS")))

est_combined <- lapply(1:n_conditions, function(i) {
  # Extract all replications for the given condition and convert them to a dataframe
  condition_results <- do.call(rbind, lapply(1:n_rep, function(j) {
    data.frame(
      ekc = est[[i]][["ekc"]][[j]],
      pa_pca = est[[i]][["pa_pca"]][[j]],
      cd = est[[i]][["cd"]][[j]],
      kgc = est_kgc[[i]][["kgc"]][[j]]
    )
  }))
  
  # Return the dataframe for this condition
  return(condition_results)
})

#saveRDS(est_combined,file = "estimations.RDS")
est_combined <- readRDS(file = "estimations.RDS")

# Define the list of method names
methods <- c("ekc", "cd", "pa_pca","kgc")

####---------------- ACCURACY ANALYZES ----------------####
# Step 1: Calculate the accuracy for each method
acc <- lapply(1:nrow(current_condition), function(i) {
  true_k <- data.frame(current_condition$k) # True number of factors
  # Compare the results dataframe with the true values
  ifelse(est_combined[[i]] == true_k[i, ], 1, 0)
})
# Step 2: Calculate the accuracy percentage for each method
acc_percent <- lapply(1:nrow(current_condition), function(i) {
  (colSums(acc[[i]],na.rm = T) / n_rep) * 100
})
# Step 3: Let's prepare results for plotting.
acc <- do.call(rbind,acc_percent)
acc <- cbind(acc,current_condition)

####---------------- BIAS ANALYZES ----------------####
# Step 1: Calculate the <absolute bias> for each method
abs_bias <- lapply(1:nrow(current_condition), function(i) {
  true_k <- data.frame(current_condition$k) # True number of factors
  # Compare the results df with the true number of factors
  colSums(abs(est_combined[[i]]-true_k[i,]),na.rm = T)/n_rep
})
# Step 2: Let's prepare results for plotting.
abs_bias <- do.call(rbind,abs_bias)
abs_bias <- cbind(abs_bias,current_condition)

# Reshape the data from wide to long format
bias_long_df <- abs_bias %>%
  pivot_longer(cols = c("ekc", "cd","pa_pca","kgc"), 
               names_to = "method",
               values_to = "bias_value")

acc_long_df <- acc %>%
  pivot_longer(cols = c("ekc", "cd","pa_pca","kgc"), 
               names_to = "method",
               values_to = "acc_value")

