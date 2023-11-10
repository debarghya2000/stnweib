#' Implementing a lag-1 stochastic process model particularly for those time-series dataset which has tie case in two consecutive time points.
#' @export
#' @param x numeric variable
#' @description
#' This package will implement a new proposed discrete time and continuous state space stochastic process model based on Weibull distribution.It works very effectively when there are ties  in two consecutive time points in the user-given time series data.The proposed stochastic process here is a lag-1 process.So,it will simply ask for time-series data from the user and then it will provide the test statistic for the model and it will also provide the corresponding p-value to understand how good the proposed model fits with the data.




stnwep <- function(x)
{
  # for lamda_0 = lamda_1
  I_1 <- character()
  for(i in 1:(length(x)-1))
  {
    I_1[i] <- x[i] < x[i+1]
  }
  
  size_1 <- length(which(I_1 == "TRUE"))
  
  
  
  I_2 <- character()
  for(i in 1:(length(x)-1))
  {
    I_2[i] <- x[i] > x[i+1]
  }
  
  size_2 <- length(which(I_2 == "TRUE"))
  
  
  
  I_0 <- character()
  for(i in 1:(length(x)-1))
  {
    I_0[i] <- x[i] == x[i+1]
  }
  
  size_0 <- length(which(I_0 == "TRUE"))
  
  
  
  
  neg_log_likelihood <- function(alpha,data)
  {
    I_1_size <- character()
    for(i in 1:(length(data)-1))
    {
      I_1_size[i] <- data[i] < data[i+1]
    }
    
    n_1 <- length(which(I_1_size == "TRUE"))
    
    
    
    I_2_size <- character()
    for(i in 1:(length(data)-1))
    {
      I_2_size[i] <- data[i] > data[i+1]
    }
    
    n_2 <- length(which(I_2_size == "TRUE"))
    
    
    
    I_0_size <- character()
    for(i in 1:(length(data)-1))
    {
      I_0_size[i] <- data[i] == data[i+1]
    }
    
    n_0 <- length(which(I_0_size == "TRUE"))
    
    
    
    
    a <- log(data[which(I_1_size == "TRUE")])
    b <- log(data[which(I_2_size == "TRUE")])
    L <- c(a,b)
    k_1 <- (data[which(I_1_size == "TRUE")])^(alpha)
    k_2 <- 2*(data[which(I_1_size == "TRUE")+1])^(alpha)
    k_3 <- 2*(data[which(I_2_size == "TRUE")])^(alpha)
    k_4 <- (data[which(I_2_size == "TRUE") + 1])^(alpha)
    k_5 <- 3*(data[which(I_0_size == "TRUE")])^(alpha)
    
    k_6 <- numeric()
    for(i in 2:(length(data)-1))
    {
      k_6[i-1] <-  2*((data[i])^(alpha))
    }
    
    G_new <- sum(k_6)
    g_1_D <- sum(k_1) + sum(k_2) + sum(k_3) + sum(k_4) + sum(k_5) - G_new
    h_alpha <- ((n_1 + n_2 + 1)* log(alpha)) + ((n_1 + n_2 + 1)*(log(n_1 + n_2 + 1) - log(g_1_D) - 1 )) + ((alpha-1)*(log(data[1]) + sum(L)))
    output <- -h_alpha
    return(output)
    
  }
  
  
  results = optim(par = 2,
                  data = x,
                  fn = neg_log_likelihood,
                  method = "L-BFGS-B",
                  lower = c(0),
                  upper = c(Inf),
                  hessian = TRUE)
  
  
  
  
  k_1_dat <- (x[which(I_1 == "TRUE")])^(results$par)
  k_2_dat <- 2*(x[which(I_1 == "TRUE")+1])^(results$par)
  k_3_dat <- 2*(x[which(I_2 == "TRUE")])^(results$par)
  k_4_dat <- (x[which(I_2 == "TRUE") + 1])^(results$par)
  k_5_dat <- 3*(x[which(I_0 == "TRUE")])^(results$par)
  
  k_6_dat <- numeric()
  for(i in 2:(length(x)-1))
  {
    k_6_dat[i-1] <-  2*((x[i])^(results$par))
  }
  
  G <- sum(k_6_dat)
  g_1_dat <- sum(k_1_dat) + sum(k_2_dat) + sum(k_3_dat) + sum(k_4_dat) + sum(k_5_dat) - G
  
  mle_alpha_new <- results$par
  mle_lambda <- (size_1 + size_2 + 1)/g_1_dat
  
  
  # for lambda_0 != lambda_1
  
  
  # beta <- (0.15/0.04)^(1/3)
  
  I_1_beta <- character()
  for(i in 1:(length(x)-1))
  {
    I_1_beta[i] <- x[i] < x[i+1]
  }
  
  size_1_beta <- length(which(I_1_beta == "TRUE"))
  
  
  
  I_2_beta <- character()
  for(i in 1:(length(x)-1))
  {
    I_2_beta[i] <- x[i] > x[i+1]
  }
  
  size_2_beta <- length(which(I_2_beta == "TRUE"))
  
  
  
  I_0_beta <- character()
  for(i in 1:(length(x)-1))
  {
    I_0_beta[i] <- x[i] == x[i+1]
  }
  
  size_0_beta <- length(which(I_0_beta == "TRUE"))
  
  
  
  
  
  
  
  
  
  gamma_beta <- 0.15/0.04
  
  
  g_2_fixed_alpha_gamma <- sum((x[which(I_1_beta == TRUE)])^3) + (1+gamma_beta)*(sum((x[which(I_1_beta == TRUE) + 1])^(3)) + sum((x[which(I_2_beta == TRUE)])^(3))) + (gamma_beta)*sum((x[which(I_2_beta == TRUE) + 1])^(3)) + ((gamma_beta)^(2) + gamma_beta + 1)*sum((x[which(I_0_beta == TRUE)])^(3)) - (1 + gamma_beta)*sum((x[seq(from = 2 , to = 74)])^(3))
  
  
  lambda_1_cap_fixed_gamma_alpha <- (size_1_beta + size_2_beta + 1)/g_2_fixed_alpha_gamma
  
  
  neg_lik_new <- function(parameter , data)
  {
    gamma <- parameter[1]
    alpha <- parameter[2]
    # beta <- (0.15/0.04)^(1/3)
    I_1_beta_data <- character()
    for(i in 1:(length(data)-1))
    {
      I_1_beta_data[i] <- data[i] < data[i+1]
    }
    
    size_1_beta_data <- length(which(I_1_beta_data == "TRUE"))
    
    
    
    I_2_beta_data <- character()
    for(i in 1:(length(data)-1))
    {
      I_2_beta_data[i] <- data[i] > data[i+1]
    }
    
    size_2_beta_data <- length(which(I_2_beta_data == "TRUE"))
    
    
    
    I_0_beta_data <- character()
    for(i in 1:(length(data)-1))
    {
      I_0_beta_data[i] <- data[i] == data[i+1]
    }
    
    size_0_beta_data <- length(which(I_0_beta_data == "TRUE"))
    
    p <- numeric()
    for(i in 2:(length(data) - 1))
    {
      p[i-1] <- (data[i])^(alpha)
    }
    Z <- sum(p)
    
    
    G_2_alpha_gamma <- sum((data[which(I_1_beta_data == TRUE)])^alpha) + (1+gamma)*(sum((data[which(I_1_beta_data == TRUE) + 1])^(alpha)) + sum((data[which(I_2_beta_data == TRUE)])^(alpha))) + (gamma)*sum((data[which(I_2_beta_data == TRUE) + 1])^(alpha)) + ((gamma)^(2) + gamma + 1)*sum((data[which(I_0_beta_data == TRUE)])^(alpha)) - (1 + gamma)*Z
    H <- log(data[1]) + sum(log(c(data[which(I_1_beta_data == TRUE) + 1],data[which(I_2_beta_data == TRUE) + 1])))
    l_gamma_alpha <- ((size_1_beta_data + size_2_beta_data + 1) * (log(lambda_1_cap_fixed_gamma_alpha) + log(alpha))) - (lambda_1_cap_fixed_gamma_alpha*(G_2_alpha_gamma) ) - (size_0_beta_data - 1)*log(1 + gamma) + ((size_0_beta_data -(1/alpha)) + size_2_beta_data)*(log(gamma)) + ((alpha - 1)*H)
    save_new <- -(l_gamma_alpha)
    return(save_new)
  }
  
  
  
  
  
  results_data_2 = optim(par = c(gamma =  2  , alpha = 2 ),
                         data = x,
                         fn = neg_lik_new,
                         method = "L-BFGS-B",
                         lower = c(0,0),
                         upper = c(Inf,Inf),
                         hessian = TRUE)$par
  
  mle_gamma <- results_data_2[[1]]
  mle_alpha_new_1 <- results_data_2[[2]]
  p_estimated <- numeric()
  for(i in 2:(length(x) - 1))
  {
    p_estimated[i-1] <- (x[i])^(mle_alpha_new_1)
  }
  Z_estimated <- sum(p_estimated)
  g_2_estimated <- sum((x[which(I_1_beta == TRUE)])^(mle_alpha_new_1)) + (1+mle_gamma)*(sum((x[which(I_1_beta == TRUE) + 1])^(mle_alpha_new_1)) + sum((x[which(I_2_beta == TRUE)])^(mle_alpha_new_1))) + (mle_gamma)*sum((x[which(I_2_beta == TRUE) + 1])^(mle_alpha_new_1)) + ((mle_gamma)^(2) + mle_gamma + 1)*sum((x[which(I_0_beta == TRUE)])^(mle_alpha_new_1)) - (1 + mle_gamma)*Z_estimated
  
  mle_lamda_1 <- (size_1_beta + size_2_beta + 1)/g_2_estimated
  
  mle_lamda_0 <- mle_gamma*mle_lamda_1
  
  
  all <- c(mle_alpha_new , mle_lambda , mle_alpha_new_1 , mle_lamda_0 , mle_lamda_1)
  weibull <- function(n,lambda_not,lambda_one,alpha)
  {
    U <- runif(n+1 , 0 , 1)
    d <- numeric()
    for(i in 1:n)
    {
      d[i] <- min((((-1/lambda_not)*(log(U[i+1])))^(1/alpha)) , (((-1/lambda_one)*(log(U[i])))^(1/alpha)))
    }
    
    return(d)
  }
  
  
  samp1 <- list()
  for(i in 1:1000)
  {
    samp1[[i]] <- sort(weibull(length(x),all[2],all[2],all[1]))
  }
  
  est1 <- matrix(0 , nrow = 1000 , ncol = length(x))
  for(i in 1:1000)
  {
    for(j in 1:length(x))
    {
      est1[i,j] <- samp1[[i]][j]
    }
  }
  est1_new <- numeric()
  for(k in 1:length(x))
  {
    est1_new[k] <- mean(est1[,k])
  }
  teststat1 <-numeric()
  for(l in 1:1000)
  {
    teststat1[l] <- max(abs(samp1[[l]]- est1_new))
  }
  
  ord_test <- sort(teststat1)
  cn <- ord_test[990]
  teststat_org <- max(abs(sort(weibull(length(x),all[2],all[2],all[1])) - est1_new))
  p1 <- length(which(teststat1 > teststat_org))/length(teststat1)
  out <- c(teststat_org , p1)
 
  
  # if(teststat_org > cn)
  #  {
  # print("Data is coming from WEP where lamda_0 != lambda_1")
  #}else
  #{
  #  print("Data is coming from WEP where lamda_0!= lambda_1")
  #}
  
  samp2 <- list()
  for(m in 1:1000)
  {
    samp2[[m]] <- sort(weibull(length(x),all[4],all[5],all[3]))
  }
  
  est2 <- matrix(0 , nrow = 1000 , ncol = length(x))
  for(n in 1:1000)
  {
    for(o in 1:length(x))
    {
      est2[n,o] <- samp2[[n]][o]
    }
  }
  est2_new <- numeric()
  for(p in 1:length(x))
  {
    est2_new[p] <- mean(est2[,p])
  }
  teststat2 <-numeric()
  for(q in 1:1000)
  {
    teststat2[q] <- max(abs(samp2[[q]]- est2_new))
  }
  
  ord_test2 <- sort(teststat2)
  cn2 <- ord_test2[990]
  teststat_org2 <- max(abs(samp2[[1]] - est2_new))
  p2 <- length(which(teststat2 > teststat_org2))/length(teststat2)
  out2 <- c(teststat_org2 , p2)
  print(paste("Test Statistic for Model WEP(alpha,lambda,lambda) =", out[1] , " and P-value =",out[2])  )
  print(paste("Test Statistic for Model WEP(alpha,lambda_0,lambda_1) =", out2[1] , " and P-value =",out2[2])  )
  par(mfrow = c(1,2))
  hist(teststat1 , col = "red" , xlab = "Test Statistic" , main = expression(paste("WEP (", alpha,",",lambda,",",lambda ,") Model")) )
  hist(teststat2 , col = "blue" , xlab = "Test Statistic" , main = expression(paste("WEP (" ,alpha,",",lambda[0],",",lambda[1] ,") Model")))
  print(out)
}

