source("data & fitting.R")

AIC_import = c("n", "beta", "thresholds", "cor_X", "active_p", "simulation_data_folder",
               "AIC_model_folder", "ordinal_data", "AIC_ordinal_fitting")
AIC_method = function(seed_num_data){
  if (!file.exists(simulation_data_folder)) {
    dir.create(simulation_data_folder, recursive = TRUE)
  }
  if (!file.exists(AIC_model_folder)) {
    dir.create(AIC_model_folder, recursive = TRUE)
  }
  if (file.exists(paste0(simulation_data_folder, seed_num_data, ".rds"))){
    simulation_data = readRDS(paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  else{
    simulation_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_data)
    saveRDS(simulation_data, file = paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  fitting_data = subset(simulation_data, select = -Y_star)
  if (file.exists(paste0(AIC_model_folder, seed_num_data, ".rds"))){
    AIC_model = readRDS(paste0(AIC_model_folder, seed_num_data, ".rds"))
  }
  else{
    AIC_model = AIC_ordinal_fitting(fitting_data)
    saveRDS(AIC_model, file = paste0(AIC_model_folder, seed_num_data, ".rds"))
  }
  selected_AIC = labels(terms(AIC_model))
  fitting_data_filtered = fitting_data[, selected_AIC, drop = FALSE]
  Y_star_pre = as.matrix(fitting_data_filtered, nrow=n, byrow=TRUE) %*% AIC_model$beta
  Y_star_pre = as.vector(Y_star_pre[,1])
  tau_hat = AIC_model$alpha
  Y_hat = rep(0, n)
  for (l in 1:n){
    if (Y_star_pre[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_pre[l] < tau_hat))
    }
  }
  Acc_train = mean(Y_hat == fitting_data[, (p+1)])
  
  result_AIC = rep(0, 5)
  if (identical(selected_AIC, paste0("X", active_p))){
    result_AIC[1] = 1
  }
  result_AIC[2] = sum(selected_AIC %in% paste0("X", active_p))
  result_AIC[3] = (length(selected_AIC) - sum(selected_AIC %in% paste0("X", active_p)))
  if (all(paste0("X", active_p) %in% selected_AIC)){
    result_AIC[4] = 1
  }
  result_AIC[5] = Acc_train
  return(list(
    result_AIC = result_AIC,
    selected_AIC = selected_AIC,
    beta_hat = AIC_model$beta,
    tau_hat = tau_hat
  ))
}

BIC_import = c("n", "beta", "thresholds", "cor_X", "active_p", "simulation_data_folder",
                 "BIC_model_folder", "ordinal_data", "BIC_ordinal_fitting")
BIC_method = function(seed_num_data){
  if (!file.exists(simulation_data_folder)) {
    dir.create(simulation_data_folder, recursive = TRUE)
  }
  if (!file.exists(BIC_model_folder)) {
    dir.create(BIC_model_folder, recursive = TRUE)
  }
  if (file.exists(paste0(simulation_data_folder, seed_num_data, ".rds"))){
    simulation_data = readRDS(paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  else{
    simulation_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_data)
    saveRDS(simulation_data, file = paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  fitting_data = subset(simulation_data, select = -Y_star)
  if (file.exists(paste0(BIC_model_folder, seed_num_data, ".rds"))){
    BIC_model = readRDS(paste0(BIC_model_folder, seed_num_data, ".rds"))
  }
  else{
    BIC_model = BIC_ordinal_fitting(fitting_data)
    saveRDS(BIC_model, file = paste0(BIC_model_folder, seed_num_data, ".rds"))
  }
  selected_BIC = labels(terms(BIC_model))
  fitting_data_filtered = fitting_data[, selected_BIC, drop = FALSE]
  Y_star_pre = as.matrix(fitting_data_filtered, nrow=n, byrow=TRUE) %*% BIC_model$beta
  Y_star_pre = as.vector(Y_star_pre[,1])
  tau_hat = BIC_model$alpha
  Y_hat = rep(0, n)
  for (l in 1:n){
    if (Y_star_pre[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_pre[l] < tau_hat))
    }
  }
  Acc_train = mean(Y_hat == fitting_data[, (p+1)])
  
  result_BIC = rep(0, 5)
  if (identical(selected_BIC, paste0("X", active_p))){
    result_BIC[1] = 1
  }
  result_BIC[2] = sum(selected_BIC %in% paste0("X", active_p))
  result_BIC[3] = (length(selected_BIC) - sum(selected_BIC %in% paste0("X", active_p)))
  if (all(paste0("X", active_p) %in% selected_BIC)){
    result_BIC[4] = 1
  }
  result_BIC[5] = Acc_train
  return(list(
    result_BIC = result_BIC,
    selected_BIC = selected_BIC,
    beta_hat = BIC_model$beta,
    tau_hat = tau_hat
  ))
}

ordinalbayes_simulation_import = c("thresholds", "simulation_data_folder", "OB_model_folder", "n", "beta", "cor_X", "burn_in", "active_p")
ordinalbayes_simulation = function(seed_num){
  library(ordinalbayes)
  
  seed_num_data = seed_num
  seed_num_fit = seed_num
  K = length(thresholds) + 1
  if (!file.exists(simulation_data_folder)) {
    dir.create(simulation_data_folder, recursive = TRUE)
  }
  if (!file.exists(OB_model_folder)) {
    dir.create(OB_model_folder, recursive = TRUE)
  }
  if (file.exists(paste0(simulation_data_folder, seed_num_data, ".rds"))){
    simulation_data = readRDS(paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  else{
    simulation_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_data)
    saveRDS(simulation_data, file = paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  train_data = subset(simulation_data, select = -Y_star)
  train_data$Y = factor(train_data$Y, ordered=TRUE)
  p = length(train_data[1,]) - 1
  K = length(thresholds) + 1
  if (file.exists(paste0(OB_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))){
    OB_model = readRDS(paste0(OB_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))
  }
  else{
    OB_model=ordinalbayes(Y ~ 1,
                          data = train_data,
                          x = train_data[, 1:p],
                          center = FALSE,
                          scale = FALSE,
                          model = "regressvi",
                          gamma.ind = "fixed",
                          pi.fixed = 0.5,
                          seed = seed_num_fit,
                          burnInSteps = burn_in,
                          nChains = 1,
                          adaptSteps = 0,
                          numSavedSteps = 2000
    )
    saveRDS(OB_model, file = paste0(OB_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))
  }
  
  coefficients = coef(OB_model)
  selected = which(coefficients$gamma > 0.5)
  selected_bayes = paste0("X", selected)
  beta_hat = rep(0, p)
  beta_hat[selected] = coefficients$beta[selected]
  gamma_hat = rep(0, p)
  gamma_hat[selected] = 1
  tau_hat = coefficients$alpha
  result_bayes = rep(0, 5)
  if (identical(selected_bayes, paste0("X", active_p))){
    result_bayes[1] = 1
  }
  result_bayes[2] = sum(selected_bayes %in% paste0("X", active_p))
  result_bayes[3] = (length(selected_bayes) - sum(selected_bayes %in% paste0("X", active_p)))
  if (all(paste0("X", active_p) %in% selected_bayes)){
    result_bayes[4] = 1
  }
  
  n_train = length(train_data[, 1])
  Y_star_hat = as.matrix(train_data[, 1:p], nrow=n_train, byrow=TRUE) %*% beta_hat
  Y_star_hat = as.vector(Y_star_hat[,1])
  Y_hat = rep(0, n_train)
  for (l in 1:n_train){
    if (Y_star_hat[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_hat[l] < tau_hat))
    }
  }
  ACC_train = mean(Y_hat == train_data[, (p+1)])
  result_bayes[5] = ACC_train
  
  return(list(
    result_stat = result_bayes,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    tau_hat = tau_hat
  ))
}

ordinalbayes_test_import = c("seed_num_data", "seed_num_test_data", "n", "beta", "thresholds", "cor_X", "testing_data_folder", "OB_model_folder", "ordinal_data")
ordinalbayes_testing = function(seed_num){
  library(ordinalbayes)
  
  seed_num_data = seed_num
  seed_num_fit = seed_num
  K = length(thresholds) + 1
  if (!file.exists(testing_data_folder)) {
    dir.create(testing_data_folder, recursive = TRUE)
  }
  if (!file.exists(OB_model_folder)) {
    stop('OB_model_folder does not exist.')
  }
  if (file.exists(paste0(testing_data_folder, seed_num_test_data, ".rds"))){
    testing_data = readRDS(paste0(testing_data_folder, seed_num_test_data, ".rds"))
  }
  else{
    testing_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_test_data)
    saveRDS(testing_data, file = paste0(testing_data_folder, seed_num_test_data, ".rds"))
  }
  if (file.exists(paste0(OB_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))){
    OB_model = readRDS(paste0(OB_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))
  }
  else{
    stop('OB_model does not exist.')
  }
  coefficients = coef(OB_model)
  p = length(coefficients$gamma)
  selected = which(coefficients$gamma > 0.5)
  selected_bayes = paste0("X", selected)
  beta_hat = rep(0, p)
  beta_hat[selected] = coefficients$beta[selected]
  gamma_hat = rep(0, p)
  gamma_hat[selected] = 1
  tau_hat = coefficients$alpha
  testing_data = subset(testing_data, select = -Y_star)
  Y_star_pre = as.matrix(testing_data[,1:p], nrow=n, byrow=TRUE) %*% beta_hat
  Y_star_pre = as.vector(Y_star_pre[,1])
  Y_hat = rep(0, n)
  for (l in 1:n){
    if (Y_star_pre[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_pre[l] < tau_hat))
    }
  }
  ACC_test = mean(Y_hat == testing_data[, (p+1)])
  return(list(
    classification_table = table(testing_data[, (p+1)], Y_hat),
    ACC_test = ACC_test
  ))
}

bayes_import = c("n", "beta", "thresholds", "cor_X", "active_p", "a", "b", "burn_in", 
                 "simulation_data_folder", "bayes_model_folder", "ordinal_data", "bayes_ordinal_fitting")
bayes_method = function(seed_num){
  seed_num_data = seed_num
  seed_num_fit = seed_num
  K = length(thresholds) + 1
  if (!file.exists(simulation_data_folder)) {
    dir.create(simulation_data_folder, recursive = TRUE)
  }
  if (!file.exists(bayes_model_folder)) {
    dir.create(bayes_model_folder, recursive = TRUE)
  }
  if (file.exists(paste0(simulation_data_folder, seed_num_data, ".rds"))){
    simulation_data = readRDS(paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  else{
    simulation_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_data)
    saveRDS(simulation_data, file = paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  fitting_data = subset(simulation_data, select = -Y_star)
  save_bayes_model = paste0(bayes_model_folder, seed_num_data, "_", seed_num_fit, ".rds")
  iter = burn_in + 2000
  bayes_model = bayes_ordinal_fitting(fitting_data, iter=iter, a=a, b=b, theta=0.5, seed_num=seed_num_fit, save=save_bayes_model)
  p = length(bayes_model$gamma_record[, 1])
  gamma_sum = rep(0, p)
  beta_sum = rep(0, p)
  tau_sum = rep(0, (K-1))
  for (j in (burn_in + 1):iter){
    gamma_sum = gamma_sum + bayes_model$gamma_record[, j]
    beta_sum = beta_sum + bayes_model$beta_record[, j]
    tau_sum = tau_sum + bayes_model$tau_record[, j]
  }
  beta_hat = rep(0, p)
  gamma_hat = rep(0, p)
  tau_hat = tau_sum / (iter - burn_in)
  selected = which(gamma_sum >= (iter - burn_in) / 2)
  beta_hat[selected] = beta_sum[selected]/gamma_sum[selected]
  gamma_hat[selected] = 1
  Y_star_pre = as.matrix(fitting_data[,1:p], nrow=n, byrow=TRUE) %*% beta_hat
  Y_star_pre = as.vector(Y_star_pre[,1])
  Y_hat = rep(0, n)
  for (l in 1:n){
    if (Y_star_pre[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_pre[l] < tau_hat))
    }
  }
  Acc_train = mean(Y_hat == fitting_data[, (p+1)])
  
  selected_bayes = paste0("X", selected)
  result_bayes = rep(0, 5)
  if (identical(selected_bayes, paste0("X", active_p))){
    result_bayes[1] = 1
  }
  result_bayes[2] = sum(selected_bayes %in% paste0("X", active_p))
  result_bayes[3] = (length(selected_bayes) - sum(selected_bayes %in% paste0("X", active_p)))
  if (all(paste0("X", active_p) %in% selected_bayes)){
    result_bayes[4] = 1
  }
  result_bayes[5] = Acc_train
  return(list(
    result_stat = result_bayes,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    tau_hat = tau_hat
  ))
}

bayes_import_for_GRD = c("seed_num_data", "n", "beta", "thresholds", "cor_X", "active_p", "a", "b", "iter", "burn_in", 
                         "simulation_data_folder", "bayes_model_folder", "ordinal_data", "bayes_ordinal_fitting")
bayes_method_for_GRD = function(seed_num_fit){
  K = length(thresholds) + 1
  if (!file.exists(simulation_data_folder)) {
    dir.create(simulation_data_folder, recursive = TRUE)
  }
  if (!file.exists(bayes_model_folder)) {
    dir.create(bayes_model_folder, recursive = TRUE)
  }
  if (file.exists(paste0(simulation_data_folder, seed_num_data, ".rds"))){
    simulation_data = readRDS(paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  else{
    simulation_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_data)
    saveRDS(simulation_data, file = paste0(simulation_data_folder, seed_num_data, ".rds"))
  }
  fitting_data = subset(simulation_data, select = -Y_star)
  save_bayes_model = paste0(bayes_model_folder, seed_num_data, "_", seed_num_fit, ".rds")
  bayes_model = bayes_ordinal_fitting(fitting_data, iter=iter, a=a, b=b, theta=0.5, seed_num=seed_num_fit, save=save_bayes_model)
  p = length(bayes_model$gamma_record[, 1])
  gamma_sum = rep(0, p)
  beta_sum = rep(0, p)
  tau_sum = rep(0, (K-1))
  for (j in (burn_in + 1):iter){
    gamma_sum = gamma_sum + bayes_model$gamma_record[, j]
    beta_sum = beta_sum + bayes_model$beta_record[, j]
    tau_sum = tau_sum + bayes_model$tau_record[, j]
  }
  beta_hat = rep(0, p)
  gamma_hat = rep(0, p)
  tau_hat = tau_sum / (iter - burn_in)
  selected = which(gamma_sum >= (iter - burn_in) / 2)
  beta_hat[selected] = beta_sum[selected]/gamma_sum[selected]
  gamma_hat[selected] = 1
  Y_star_pre = as.matrix(fitting_data[,1:p], nrow=n, byrow=TRUE) %*% beta_hat
  Y_star_pre = as.vector(Y_star_pre[,1])
  Y_hat = rep(0, n)
  for (l in 1:n){
    if (Y_star_pre[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_pre[l] < tau_hat))
    }
  }
  Acc_train = mean(Y_hat == fitting_data[, (p+1)])
  
  selected_bayes = paste0("X", selected)
  result_bayes = rep(0, 5)
  if (identical(selected_bayes, paste0("X", active_p))){
    result_bayes[1] = 1
  }
  result_bayes[2] = sum(selected_bayes %in% paste0("X", active_p))
  result_bayes[3] = (length(selected_bayes) - sum(selected_bayes %in% paste0("X", active_p)))
  if (all(paste0("X", active_p) %in% selected_bayes)){
    result_bayes[4] = 1
  }
  result_bayes[5] = Acc_train
  return(list(
    result_stat = result_bayes,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    tau_hat = tau_hat
  ))
}

bayes_check_convergence = function(bayes_result, tau_threshold=1.1, gamma_threshold=0.06){
  if (!"mcmcse" %in% loadedNamespaces()){
    library(mcmcse)
  }
  
  num_MCMC = length(bayes_result)
  iter = length(bayes_result[[1]]$tau_record[1,])
  interval = 1000
  num_check = iter %/% interval
  
  # tau convergence
  num_tau = length(bayes_result[[1]]$tau_record[, 1])
  for (i in 1:num_MCMC){
    tau_matrix = matrix(0, iter, num_tau)
    for (j in 1:iter){
      tau_matrix[j,] <- bayes_result[[i]]$tau_record[, j]
    }
    if (i == 1){
      mcmc_tau = mcmc.list(mcmc(tau_matrix))
    }
    else{
      mcmc_tau = c(mcmc_tau, mcmc.list(mcmc(tau_matrix)))
    }
  }
  mc = mcmc.list(mcmc_tau)
  
  GRD = rep(0, num_check)
  for (j in 1:num_check){
    subset_mcmc_list = lapply(mc, function(chain){mcmc(chain[1:(interval*j), ])})
    GRD[j] = gelman.diag(subset_mcmc_list)$mpsrf
  }
  check_GRD = which(GRD < tau_threshold)
  check_GRD = check_GRD[check_GRD > 2]
  if (length(check_GRD) == 0){
    converge_tau = "tau does not converge!"
  }else{
    converge_tau = min(check_GRD) * interval
  }
  
  # gamma convergence
  converge_gamma_list = rep(0, num_MCMC)
  for (i in 1:num_MCMC){
    gamma_record = bayes_result[[i]]$gamma_record
    gamma_matrix = t(gamma_record)
    max_MCSEs = rep(0, num_check)
    
    for (k in 1:num_check){
      bn = ceiling((interval*k)^(1/3))
      an = ceiling((interval*k)/bn)
      max_MCSEs[k] = max(mcse.mat(gamma_matrix[1:(interval*k), ], size="cuberoot")[,2]) * qt(0.975, an-1)
    }
    
    check_gamma = which(max_MCSEs < gamma_threshold)
    if (length(check_gamma) == 0){
      converge_gamma_list[i] = "gamma does not converge!"
    }else{
      converge_gamma_list[i] = min(setdiff(check_gamma, c(1, 2))) * interval
    }
  }
  
  if ("gamma does not converge!" %in% converge_gamma_list){
    converge_gamma = "gamma does not converge!"
  }else{
    converge_gamma = max(converge_gamma_list)
  }
  
  if (is.numeric(converge_tau) && is.numeric(converge_gamma)){
    burn_in = max(converge_tau, converge_gamma)
  }else if (is.numeric(converge_tau)){
    burn_in = converge_gamma
  }else if (is.numeric(converge_gamma)){
    burn_in = converge_tau
  }else{
    burn_in = "tau and gamma does not converge!"
  }
  
  return(list(
    burn_in = burn_in,
    tau_GRD = GRD
  ))
}

bayes_test_import = c("seed_num_data", "seed_num_test_data", "n", "beta", "thresholds", "cor_X", "burn_in", 
                      "testing_data_folder", "bayes_model_folder", "ordinal_data")
bayes_testing = function(seed_num_fit){
  K = length(thresholds) + 1
  if (!file.exists(testing_data_folder)) {
    dir.create(testing_data_folder, recursive = TRUE)
  }
  if (!file.exists(bayes_model_folder)) {
    stop('bayes_model_folder does not exist.')
  }
  if (file.exists(paste0(testing_data_folder, seed_num_test_data, ".rds"))){
    testing_data = readRDS(paste0(testing_data_folder, seed_num_test_data, ".rds"))
  }
  else{
    testing_data = ordinal_data(n, beta, thresholds, cor_X, seed_num=seed_num_test_data)
    saveRDS(testing_data, file = paste0(testing_data_folder, seed_num_test_data, ".rds"))
  }
  if (file.exists(paste0(bayes_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))){
    bayes_model = readRDS(paste0(bayes_model_folder, seed_num_data, "_", seed_num_fit, ".rds"))
  }
  else{
    stop('bayes_model does not exist.')
  }
  p = length(bayes_model$gamma_record[, 1])
  iter = burn_in + 2000
  gamma_sum = rep(0, p)
  beta_sum = rep(0, p)
  tau_sum = rep(0, (K-1))
  for (j in (burn_in + 1):iter){
    gamma_sum = gamma_sum + bayes_model$gamma_record[, j]
    beta_sum = beta_sum + bayes_model$beta_record[, j]
    tau_sum = tau_sum + bayes_model$tau_record[, j]
  }
  beta_hat = rep(0, p)
  gamma_hat = rep(0, p)
  tau_hat = tau_sum / (iter - burn_in)
  selected = which(gamma_sum >= (iter - burn_in) / 2)
  beta_hat[selected] = beta_sum[selected]/gamma_sum[selected]
  gamma_hat[selected] = 1
  testing_data = subset(testing_data, select = -Y_star)
  Y_star_pre = as.matrix(testing_data[,1:p], nrow=n, byrow=TRUE) %*% beta_hat
  Y_star_pre = as.vector(Y_star_pre[,1])
  Y_hat = rep(0, n)
  for (l in 1:n){
    if (Y_star_pre[l] > tau_hat[K-1]){
      Y_hat[l] = K
    }
    else{
      Y_hat[l] = min(which(Y_star_pre[l] < tau_hat))
    }
  }
  ACC_test = mean(Y_hat == testing_data[, (p+1)])
  return(list(
    classification_table = table(testing_data[, (p+1)], Y_hat),
    ACC_test = ACC_test
  ))
}


var_select_performance = function(p, true_p, run1, run2, import_vector, var_select_func){
  if (!"parallel" %in% loadedNamespaces()){
    library(parallel)
  }
  
  # selected variables = true active variables
  sum_right_selection_rate = 0
  
  # AEIR: the active effect identified rate
  sum_AEIR = 0
  
  # IEIR: the inactive effect identified rate
  sum_IEIR = 0
  
  # selected variables cover all true active variables
  sum_coverage_rate = 0
  
  # collect all ACC_train
  ACC_train_collect = c()
    
  # cl = makeCluster(floor(detectCores() * (3/5)))  # 建立一個包含多核心的集群
  cl = makeCluster(5)
  clusterExport(cl, varlist = import_vector)
  result = parLapply(cl, run1:run2, var_select_func)  # 對每個元素進行並行計算
  stopCluster(cl)  # 停止集群
  run = run2 - run1 + 1
  
  for (i in 1:run){
    result_stat = result[[i]]$result_stat
    sum_right_selection_rate = sum_right_selection_rate + result_stat[1]
    sum_AEIR = sum_AEIR + result_stat[2]
    sum_IEIR = sum_IEIR + result_stat[3]
    sum_coverage_rate = sum_coverage_rate + result_stat[4]
    ACC_train_collect = c(ACC_train_collect, result_stat[5])
  }
  
  right_selection_rate = sum_right_selection_rate/run
  AEIR = sum_AEIR/(true_p*run)
  IEIR = sum_IEIR/((p - true_p)*run)
  coverage_rate = sum_coverage_rate/run
  
  performance = list(
    right_selection_rate = right_selection_rate,
    AEIR = AEIR,
    IEIR = IEIR,
    coverage_rate = coverage_rate,
    ACC_train_collect = ACC_train_collect
  )
  return(performance)
}
