ordinal_data = function(n, beta, thresholds, cor_X, seed_num=1){
  # -------------------------------------------------------------------------
  # PURPOSE:
  #
  # This function is used to generate data for an Ordinal Probit Regression  
  # Model, including explanatory variables X (n x p), latent variable 
  # Y_star (n x 1), and response variable Y (n x 1), where the 
  # explanatory variables X are correlated.
  # -------------------------------------------------------------------------
  # INPUTS:
  #
  # n:               Number of generating data.
  # beta:            p x 1 vector, which is the vector of coefficients of 
  #                  explanatory variables X.
  # thresholds:      (K-1) x 1 vector, which is the cutting points.
  #                  (the subsequent elements must be greater than the preceding 
  #                   elements.)
  # cor_X:           the correlations between two different variables.
  # seed_num:        Fix the seed number to reproduce the same output.  
  # -------------------------------------------------------------------------
  # OUTPUTS:
  #
  # The DataFrame of all data, including X, Y_star and Y.
  # -------------------------------------------------------------------------
  
  set.seed(seed_num)
  
  if (!"MASS" %in% loadedNamespaces()){
    library(MASS)
  }
  if (!"Matrix" %in% loadedNamespaces()){
    library(Matrix)
  }
  
  # check inputs
  if (!all(diff(thresholds) > 0)){
    stop('incorrect thresholds.')
  }
  if (cor_X > 1 | cor_X < 0){
    stop('only for positive cor_X.')
  }
  
  # Generate explanatory variables X (mean = 0, variance = 1), 
  # with pairwise correlation equal to cor_X.
  p = length(beta)
  X_samples = matrix(0, nrow = n, ncol = p)
  random_matrix = matrix(rnorm(n * p), nrow = n, ncol = p)
  cor_vec = rnorm(n)
  for (j in 1:p) {
    X_samples[, j] = sqrt(1-cor_X) * random_matrix[, j] + sqrt(cor_X) * cor_vec
  }
  
  # Calculate latent variable
  Y_star = X_samples %*% beta + rnorm(n, mean = 0, sd = 1)
  
  # Convert to response variable  
  Y = cut(Y_star, breaks = c(-Inf, thresholds, Inf), labels = 1:(length(thresholds) + 1), right = TRUE)
  
  # Combine data into a DataFrame
  colnames(X_samples) = paste0("X", 1:ncol(X_samples))
  simulation_data = data.frame(X_samples, Y_star = Y_star, Y = Y)
  
  return(simulation_data)
}

AIC_ordinal_fitting = function(fitting_data){
  # -------------------------------------------------------------------------
  # PURPOSE:
  #
  # This is an implementation of stepwise variable selection using AIC. 
  # The response model is assumed to be the Ordinal Probit Regression
  # Model and the variables can be correlated. Here the no-intercept  
  # term is defined in the model.
  # -------------------------------------------------------------------------
  # INPUTS:
  #
  # fitting_data:    An n x (p+1) DataFrame containing the candidate variables X1, 
  #                  X2,..., Xp and the ordinal response variable Y.
  #                  (There are K ordered levels, 1,2,...,K, in Y.)
  # -------------------------------------------------------------------------
  # OUTPUTS:
  #
  # The final selected model.
  # -------------------------------------------------------------------------
  
  if (!"ordinal" %in% loadedNamespaces()){
    library(ordinal)
  }
  
  # fit Ordinal Probit Regression Model 
  fit_model = clm(factor(Y, ordered = TRUE) ~ ., data = fitting_data, link = "probit")
  # stepwise variables selection using AIC
  AIC_stepwise_model = stepAIC(fit_model, direction="both", k=2, trace=0)
  
  return(AIC_stepwise_model)
}

BIC_ordinal_fitting = function(fitting_data){
  # -------------------------------------------------------------------------
  # PURPOSE:
  #
  # This is an implementation of stepwise variable selection using BIC. 
  # The response model is assumed to be the Ordinal Probit Regression
  # Model and the variables can be correlated. Here the no-intercept  
  # term is defined in the model.
  # -------------------------------------------------------------------------
  # INPUTS:
  #
  # fitting_data:    An n x (p+1) DataFrame containing the candidate variables X1, 
  #                  X2,..., Xp and the ordinal response variable Y.
  #                  (There are K ordered levels, 1,2,...,K, in Y.)
  # -------------------------------------------------------------------------
  # OUTPUTS:
  #
  # The final selected model.
  # -------------------------------------------------------------------------
  
  if (!"ordinal" %in% loadedNamespaces()){
    library(ordinal)
  }
  
  # fit Ordinal Probit Regression Model 
  fit_model = clm(factor(Y, ordered = TRUE) ~ ., data = fitting_data, link = "probit")
  # stepwise variables selection using BIC
  BIC_stepwise_model = stepAIC(fit_model, direction="both", k=log(nrow(fitting_data)), trace=0)
  
  return(BIC_stepwise_model)
}

bayes_ordinal_fitting = function(fitting_data, iter=3000, a=Inf, b, theta=0.5, seed_num=1, save=FALSE, group_list=list()){
  # -------------------------------------------------------------------------
  # PURPOSE:
  #
  # This is an implementation of the component-wise Gibbs sampler (CGS) for  
  # variable selection. The response model is assumed to be the Ordinal Probit 
  # Regression Model and the variables can be correlated. Here the no-intercept  
  # term is defined in the model.
  # (the prior of cutting points, "tau", is set to be informative.)
  # -------------------------------------------------------------------------
  # INPUTS:
  #
  # fitting_data:    An n x (p+1) DataFrame containing the candidate variables X1, 
  #                  X2,..., Xp and the ordinal response variable Y.
  #                  (There are K ordered levels, 1,2,...,K, in Y.)
  # iter:            The iteration number for the component-wise Gibbs 
  #                  sampler of variable and element level, iteration is usually set to 
  #                  be larger than 1000.
  #                  (default = 3000)
  # a:               A pre-specified standard deviation parameter in the prior of 
  #                  cutting points. "a" must be a positive real value. 
  #                  If a = Inf, then we sample "tau" from the closed full conditional
  #                  distribution.
  #                  else, we sample "tau" and "Y_star" from the joint full 
  #                  conditional distribution.
  # b:               A pre-specified standard deviation parameter in the prior of 
  #                  coefficient. "b" must be a positive real value by assuming 
  #                  the same "b" value for all variables. 
  # theta:           The parameter in the Bernoulli prior of "gamma", which controls     
  #                  the probability of that corresponding variable is active. 
  #                  "theta" must be a real value that belongs to (0,1).
  #                  (default = 0.5) 
  # seed_num:        Fix the seed number to reproduce the same output.  
  # save:            The path saving the output.  
  # group_list:      The list containing all variable groups.
  #                  (e.x. list(c(1, 2, 3), c(5, 6)))
  # -------------------------------------------------------------------------
  # OUTPUTS:
  #
  # A list containing the following fields:
  # gamma_record: The samples of the indicators of variables.
  # beta_record:  The samples of coefficients. 
  # tau_record:   The samples of cutting points.
  # -------------------------------------------------------------------------
  
  set.seed(seed_num)
  
  if (!"msm" %in% loadedNamespaces()){
    library(msm)
  }
  
  # Initial status of the selection procedure:
  # =========================================================================
  n_size = nrow(fitting_data)            # number of samples
  p_size = ncol(fitting_data)-1          # number of variables
  X_data = fitting_data[,1:p_size]
  Y_data = fitting_data[,(p_size+1)]
  Y_data = as.numeric(Y_data)
  K = max(Y_data)                        # number of ordered levels
  Y_level_count = as.numeric(table(Y_data))
  
  if (length(group_list) != 0){
    existing_numbers = unlist(group_list)
    missing_numbers = setdiff(1:p_size, existing_numbers)
    for (num in missing_numbers) {
      group_list = append(group_list, list(c(num)))
    }
    g_size = length(group_list)
  }
  beta = rep(0, p_size)
  gamma = rep(0, p_size)
  if (a == Inf){
    tau = sort(rnorm(K-1, 0, 3))
  }else{
    tau = sort(rnorm(K-1, 0, 3*a))
  }
  
  # Record the posterior samples
  # =========================================================================
  beta_record = array(0, dim = c(p_size, iter))
  gamma_record = array(0, dim = c(p_size, iter))
  tau_record = array(0, dim = c(K-1, iter))
  
  # Initial values of latent variables
  # =========================================================================
  Y_star_data = rep(0, n_size)
  for (i in 1:n_size){
    if (Y_data[i] == 1){
      Y_star_data[i] = -abs(rnorm(1)) + tau[1]
    }
    else if (Y_data[i] == K){
      Y_star_data[i] = abs(rnorm(1)) + tau[K-1]
    }
    else{
      lower = tau[Y_data[i]-1]
      upper = tau[Y_data[i]]
      Y_star_data[i] = runif(1, lower, upper)
    }
  }
  
  temp_start = 1
  if (save != FALSE){
    if (!file.exists(dirname(save))){
      dir.create(dirname(save), recursive = TRUE)
    }
    if (file.exists(save)){
      bayes_model = readRDS(save)
      temp_iter = length(bayes_model$gamma_record[1, ])
      if (temp_iter >= iter){
        return(bayes_model)
      }
      beta_record[, 1:temp_iter] = bayes_model$beta_record
      gamma_record[, 1:temp_iter] = bayes_model$gamma_record
      tau_record[, 1:temp_iter] = bayes_model$tau_record
      temp_start = temp_iter + 1
      
      beta = beta_record[, temp_iter]
      gamma = gamma_record[, temp_iter]
      tau = tau_record[, temp_iter]
      for (i in 1:n_size){
        mean_Y_star = sum(X_data[i,] * beta)
        if (Y_data[i] == 1){
          Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=-Inf, upper=tau[1])
        }
        else if (Y_data[i] == K){
          Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=tau[K-1], upper=Inf)
        }
        else{
          lower = tau[Y_data[i]-1]
          upper = tau[Y_data[i]]
          Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=lower, upper=upper)
        }
      }
    }
  }
  
  # Main code for Component-wise Gibbs Sampler 
  # =========================================================================
  for (counter in temp_start:iter) {
    # Draw gamma, beta
    if (length(group_list) == 0){
      for (j in 1:p_size){
        R_j = Y_star_data - as.matrix(X_data, nrow=n_size, byrow=TRUE) %*% beta + X_data[, j] * beta[j]  # current residual, n x 1 vector 
        mu_tilde_j = b^2 * sum(R_j * X_data[, j]) / (1 + b^2 * sum(X_data[, j] * X_data[, j]))
        sigma_2_tilde_j = b^2 / (1 + b^2 * sum(X_data[, j] * X_data[, j]))
        sigma_tilde_j = sqrt(sigma_2_tilde_j)
        G_j = sigma_tilde_j / b * exp(1/2 * mu_tilde_j^2 / sigma_2_tilde_j)
        # compute accepted probability
        if (G_j > 10^20){
          p_j = 1
        }
        else{
          p_j = theta * G_j / (theta * G_j + (1 - theta))
        }
        # decide whether a variable is active or inactive
        unif = runif(1)
        if (p_j < unif){
          beta[j] = 0
          gamma[j] = 0
        }
        else{
          beta[j] = rnorm(1, mu_tilde_j, sigma_tilde_j)
          gamma[j] = 1
        }
      }
    }else{
      for (j in 1:g_size){
        group_j = group_list[[j]]
        len_j = length(group_j)
        R_j = Y_star_data - as.matrix(X_data[, -group_j], nrow=n_size, byrow=TRUE) %*% beta[-group_j]  # current residual, n x 1 vector
        X_j = as.matrix(X_data[, group_j], nrow=n_size, byrow=TRUE)
        sigma_matrix_inverse_j = (t(X_j) %*% X_j) + 1/(b^2) * diag(len_j)
        sigma_matrix_j = solve(sigma_matrix_inverse_j)
        mu_tilde_j = sigma_matrix_j %*% t(X_j) %*% R_j
        G_j = 1/(b^len_j) * sqrt(determinant(sigma_matrix_j, logarithm=FALSE)$modulus) * exp((1/2)*t(mu_tilde_j)%*%sigma_matrix_inverse_j%*%mu_tilde_j)
        # compute accepted probability
        if (G_j > 10^20){
          p_j = 1
        }
        else{
          p_j = theta * G_j / (theta * G_j + (1 - theta))
        }
        # decide whether a variable is active or inactive
        unif = runif(1)
        if (p_j < unif){
          for (k in group_j){
            gamma[k] = 0
            beta[k] = 0
          }
        }
        else{
          for (k in group_j){
            gamma[k] = 1
            if (len_j == 1){
              beta[k] = rnorm(1, as.numeric(mu_tilde_j), sqrt(as.numeric(sigma_matrix_j)))
              break
            }
            R_k = Y_star_data - as.matrix(X_data, nrow=n_size, byrow=TRUE) %*% beta + X_data[, k] * beta[k]  # current residual, n x 1 vector 
            mu_tilde_k = b^2 * sum(R_k * X_data[, k]) / (1 + b^2 * sum(X_data[, k] * X_data[, k]))
            sigma_2_tilde_k = b^2 / (1 + b^2 * sum(X_data[, k] * X_data[, k]))
            sigma_tilde_k = sqrt(sigma_2_tilde_k)
            beta[k] = rnorm(1, mu_tilde_k, sigma_tilde_k)
          }
        }
      }
    }
    
    # Draw tau
    if (a == Inf){
      for (k in 1:(K-1)){
        lower = max(Y_star_data[Y_data == k])
        upper = min(Y_star_data[Y_data == k+1])
        tau[k] = runif(1, lower, upper)
      }
    }
    else{
      tau_star = rep(0, (K-1))
      for (k in 1:(K-1)){
        if (k == 1){
          lower = -10^20
          upper = tau[k+1]
        }
        else if (k == (K-1)){
          lower = tau_star[k-1]
          upper = 10^20
        }
        else{
          lower = tau_star[k-1]
          upper = tau[k+1]
        }
        tau_star[k] = rtnorm(1, mean=tau[k], sd=a, lower=lower, upper=upper)
      }
      r = exp(-1/(2*a^2) * (sum(tau_star * tau_star) - sum(tau * tau)))
      for (i in 1:n_size){
        mean_Y_star = sum(X_data[i,] * beta)
        if (Y_data[i] == 1){
          r = r * (pnorm(tau_star[Y_data[i]], mean_Y_star, 1)) / 
                  (pnorm(tau[Y_data[i]], mean_Y_star, 1))
        }
        else if (Y_data[i] == K){
          r = r * (1 - pnorm(tau_star[Y_data[i]-1], mean_Y_star, 1)) / 
                  (1 - pnorm(tau[Y_data[i]-1], mean_Y_star, 1))
        }
        else{
          r = r * (pnorm(tau_star[Y_data[i]], mean_Y_star, 1) - 
                   pnorm(tau_star[Y_data[i]-1], mean_Y_star, 1)) / 
                  (pnorm(tau[Y_data[i]], mean_Y_star, 1) - 
                   pnorm(tau[Y_data[i]-1], mean_Y_star, 1))
        }
        if (is.na(r)){
          break
        }
        if (is.infinite(r)){
          break
        }
        if (r == 0){
          break
        }
      }
      if (!is.na(r)){
        r = min(1, r)
        unif = runif(1)
        if (r >= unif) {
          tau = tau_star
        }
      }
    }
    
    # Draw Y_star
    for (i in 1:n_size){
      mean_Y_star = sum(X_data[i,] * beta)
      if (Y_data[i] == 1){
        Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=-Inf, upper=tau[1])
      }
      else if (Y_data[i] == K){
        Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=tau[K-1], upper=Inf)
      }
      else{
        lower = tau[Y_data[i]-1]
        upper = tau[Y_data[i]]
        Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=lower, upper=upper)
      }
    }
    
    # record
    beta_record[, counter] = beta
    gamma_record[, counter] = gamma
    tau_record[, counter] = tau
    
    if (save != FALSE){
      if (counter %% 2000 == 0) {
        result = list(
          gamma_record = gamma_record[, 1:counter],
          beta_record = beta_record[, 1:counter],
          tau_record = tau_record[, 1:counter]
        )
        saveRDS(result, file = save)
      }
    }
  }
  
  result = list(
    gamma_record = gamma_record,
    beta_record = beta_record,
    tau_record = tau_record
  )
  if (save != FALSE){
    saveRDS(result, file = save)
  }
  return(result)
}