if (!require("mvtnorm", quietly = TRUE)) install.packages("mvtnorm")
if (!require("Gmedian", quietly = TRUE)) install.packages("Gmedian")
if (!require("flare", quietly = TRUE)) install.packages("flare")
if (!require("extraDistr", quietly = TRUE)) install.packages("extraDistr")

library(mvtnorm)
library(Gmedian)
library(flare)
library(extraDistr)

#' @param dat a list of data matrices, with each element from each group
#' @param K the number of groups
#' @param N a vector contains the number of observation in each group
#' @param p dimensions of the data
#' @param n_boot the number of bootstrap iterations
MANOVA_Median <- function (dat, K, N, p, n_boot=400) {

  median_list <- vector("list", K)
  for (k in 1:K) {
    median_list[[k]] <- drop(Weiszfeld(dat[[k]])$median)
  }

  dif_matrix <- matrix(0, K, K)
  for (k1 in 1:(K-1)) {
    for (k2 in (k1+1):K) {
      tmp1 <- abs(median_list[[k1]]-median_list[[k2]])
      dif_matrix[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp1)
    }
  }
  test_stat_median <- max(dif_matrix)

  ##standardized median
  median_sd_est <- vector("list", K)
  for (k in 1:K) {
    Xk_central <- dat[[k]]-rep(1,N[k])%o%median_list[[k]]
    norm_Xk <- sqrt(diag(Xk_central%*%t(Xk_central)))
    c_Xk <- mean(1/norm_Xk)
    B_Xk <- (t(Xk_central)%*%diag(1/norm_Xk))%*%t(t(Xk_central)%*%diag(1/norm_Xk))/N[k]
    median_sd_est[[k]] <- B_Xk/c_Xk^2
  }

  dif_matrix_std <- matrix(0, K, K)
  for (k1 in 1:(K-1)) {
    for (k2 in (k1+1):K) {
      Bxy <- (N[k2]*median_sd_est[[k1]]+N[k1]*median_sd_est[[k2]])/(N[k1]+N[k2])
      bxy <- sqrt(diag(Bxy))
      tmp2 <- abs(median_list[[k1]]-median_list[[k2]])/bxy
      dif_matrix_std[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp2)
    }
  }
  test_stat_median_std <- max(dif_matrix_std)

  #### bootstrap
  res_median_boot <- rep(0, n_boot)
  res_median_std_boot <- rep(0, n_boot)
  res_median_std_boot_old <- rep(0, n_boot)

  for(b in 1:n_boot){
    dat_boot <- vector("list", K)
    median_boot_list <- vector("list", K)
    for (k in 1:K) {
      Z_Xk <- rsign(N[k])%o%rep(1,p)
      dat_boot[[k]] <- (dat[[k]]-rep(1,N[k])%o%median_list[[k]])*Z_Xk
      median_boot_list[[k]] <- drop(Weiszfeld(dat_boot[[k]])$median)
    }

    dif_matrix_boot <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        tmp3 <- abs(median_boot_list[[k1]]-median_boot_list[[k2]])
        dif_matrix_boot[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp3)
      }
    }
    res_median_boot[b] <- max(dif_matrix_boot)

    median_sd_est_boot <- vector("list", K)
    for (k in 1:K) {
      Yk_central <- dat_boot[[k]]-rep(1,N[k])%o%median_boot_list[[k]]
      norm_Yk <- sqrt(diag(Yk_central%*%t(Yk_central)))
      c_Yk <- mean(1/norm_Yk)
      B_Yk <- (t(Yk_central)%*%diag(1/norm_Yk))%*%t(t(Yk_central)%*%diag(1/norm_Yk))/N[k]
      median_sd_est_boot[[k]] <- B_Yk/c_Yk^2
    }

    dif_matrix_std_boot <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        Bxy_boot <- (N[k2]*median_sd_est_boot[[k1]]+N[k1]*median_sd_est_boot[[k2]])/(N[k1]+N[k2])
        bxy_boot <- sqrt(diag(Bxy_boot))
        tmp4 <- abs(median_boot_list[[k1]]-median_boot_list[[k2]])/bxy_boot
        dif_matrix_std_boot[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp4)
      }
    }
    res_median_std_boot[b] <- max(dif_matrix_std_boot)

    dif_matrix_std_boot_old <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        Bxy <- (N[k2]*median_sd_est[[k1]]+N[k1]*median_sd_est[[k2]])/(N[k1]+N[k2])
        bxy <- sqrt(diag(Bxy))
        tmp5 <- abs(median_boot_list[[k1]]-median_boot_list[[k2]])/bxy
        dif_matrix_std_boot_old[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp5)
      }
    }
    res_median_std_boot_old[b] <- max(dif_matrix_std_boot_old)
  }

  return(list("test_stat_median"=test_stat_median,
              "test_stat_median_std"=test_stat_median_std,
              "res_median_boot"=res_median_boot,
              "res_median_std_boot"=res_median_std_boot,
              "res_median_std_boot_old"=res_median_std_boot_old))
}


MANOVA_Mean <- function (dat, K, N, p, n_boot=400) {
  mean_list <- vector("list", K)
  for (k in 1:K) {
    mean_list[[k]] <- apply(dat[[k]], 2, mean)
  }

  mean_dif_matrix <- matrix(0, K, K)
  for (k1 in 1:(K-1)) {
    for (k2 in (k1+1):K) {
      tmp1 <- abs(mean_list[[k1]]-mean_list[[k2]])
      mean_dif_matrix[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp1)
    }
  }
  test_stat_mean <- max(mean_dif_matrix)

  ##standardized mean
  mean_sd_est <- vector("list", K)
  for (k in 1:K) {
    mean_sd_est[[k]] <- cov(dat[[k]])*(N[k]-1)/N[k]
  }

  mean_dif_matrix_std <- matrix(0, K, K)
  for (k1 in 1:(K-1)) {
    for (k2 in (k1+1):K) {
      Sxy <- (N[k2]*mean_sd_est[[k1]]+N[k1]*mean_sd_est[[k2]])/(N[k1]+N[k2])
      sxy <- sqrt(diag(Sxy))
      tmp2 <- abs(mean_list[[k1]]-mean_list[[k2]])/sxy
      mean_dif_matrix_std[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp2)
    }
  }
  test_stat_mean_std <- max(mean_dif_matrix_std)

  #### bootstrap
  res_mean_boot <- rep(0, n_boot)
  res_mean_std_boot <- rep(0, n_boot)
  res_mean_std_boot_old <- rep(0, n_boot)

  for(b in 1:n_boot){
    dat_boot <- vector("list", K)
    mean_boot_list <- vector("list", K)
    for (k in 1:K) {
      e1 <- rsign(N[k])%o%rep(1,p)
      dat_boot[[k]] <- (dat[[k]]-rep(1,N[k])%o%mean_list[[k]])*e1
      mean_boot_list[[k]] <- apply(dat_boot[[k]], 2, mean)
    }

    mean_dif_matrix_boot <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        tmp3 <- abs(mean_boot_list[[k1]]-mean_boot_list[[k2]])
        mean_dif_matrix_boot[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp3)
      }
    }
    res_mean_boot[b] <- max(mean_dif_matrix_boot)

    ##standardized mean
    mean_sd_est_boot <- vector("list", K)
    for (k in 1:K) {
      mean_sd_est_boot[[k]] <- cov(dat_boot[[k]])*(N[k]-1)/N[k]
    }

    mean_dif_matrix_std_boot <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        Sxy_boot <- (N[k2]*mean_sd_est_boot[[k1]]+N[k1]*mean_sd_est_boot[[k2]])/(N[k1]+N[k2])
        sxy_boot <- sqrt(diag(Sxy_boot))
        tmp4 <- abs(mean_boot_list[[k1]]-mean_boot_list[[k2]])/sxy_boot
        mean_dif_matrix_std_boot[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp4)
      }
    }
    res_mean_std_boot[b] <- max(mean_dif_matrix_std_boot)

    mean_dif_matrix_std_boot_old <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        Sxy <- (N[k2]*mean_sd_est[[k1]]+N[k1]*mean_sd_est[[k2]])/(N[k1]+N[k2])
        sxy <- sqrt(diag(Sxy))
        tmp5 <- abs(mean_boot_list[[k1]]-mean_boot_list[[k2]])/sxy
        mean_dif_matrix_std_boot_old[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp5)
      }
    }
    res_mean_std_boot_old[b] <- max(mean_dif_matrix_std_boot_old)
  }

  return(list("test_stat_mean"=test_stat_mean,
              "test_stat_mean_std"=test_stat_mean_std,
              "res_mean_boot"=res_mean_boot,
              "res_mean_std_boot"=res_mean_std_boot,
              "res_mean_std_boot_old"=res_mean_std_boot_old))
}


MANOVA_Mean_Lin2021 <- function (dat, K, N, p, n_boot=400) {
  mean_list <- vector("list", K)
  for (k in 1:K) {
    mean_list[[k]] <- apply(dat[[k]], 2, mean)
  }

  mean_dif_matrix <- matrix(0, K, K)
  for (k1 in 1:(K-1)) {
    for (k2 in (k1+1):K) {
      tmp1 <- abs(mean_list[[k1]]-mean_list[[k2]])
      mean_dif_matrix[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp1)
    }
  }
  test_stat_mean <- max(mean_dif_matrix)

  ##standardized mean
  mean_sd_est <- vector("list", K)
  for (k in 1:K) {
    mean_sd_est[[k]] <- cov(dat[[k]])*(N[k]-1)/N[k]
  }

  mean_dif_matrix_std <- matrix(0, K, K)
  for (k1 in 1:(K-1)) {
    for (k2 in (k1+1):K) {
      Sxy <- (N[k2]*mean_sd_est[[k1]]+N[k1]*mean_sd_est[[k2]])/(N[k1]+N[k2])
      sxy <- sqrt(diag(Sxy))
      tmp2 <- abs(mean_list[[k1]]-mean_list[[k2]])/sxy
      mean_dif_matrix_std[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp2)
    }
  }
  test_stat_mean_std <- max(mean_dif_matrix_std)

  #### bootstrap
  res_mean_boot <- rep(0, n_boot)
  res_mean_std_boot <- rep(0, n_boot)

  for(b in 1:n_boot){
    dat_boot <- vector("list", K)
    mean_boot_list <- vector("list", K)
    for (k in 1:K) {
      dat_boot[[k]] <- rmvnorm(N[k], rep(0, p), mean_sd_est[[k]])
      mean_boot_list[[k]] <- apply(dat_boot[[k]], 2, mean)
    }

    mean_dif_matrix_boot <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        tmp3 <- abs(mean_boot_list[[k1]]-mean_boot_list[[k2]])
        mean_dif_matrix_boot[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp3)
      }
    }
    res_mean_boot[b] <- max(mean_dif_matrix_boot)

    mean_dif_matrix_std_boot <- matrix(0, K, K)
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {
        Sxy <- (N[k2]*mean_sd_est[[k1]]+N[k1]*mean_sd_est[[k2]])/(N[k1]+N[k2])
        sxy <- sqrt(diag(Sxy))
        tmp5 <- abs(mean_boot_list[[k1]]-mean_boot_list[[k2]])/sxy
        mean_dif_matrix_std_boot[k1, k2] <- max(sqrt(N[k1]*N[k2]/(N[k1]+N[k2]))*tmp5)
      }
    }
    res_mean_std_boot[b] <- max(mean_dif_matrix_std_boot)
  }

  return(list("test_stat_mean"=test_stat_mean,
              "test_stat_mean_std"=test_stat_mean_std,
              "res_mean_boot"=res_mean_boot,
              "res_mean_std_boot"=res_mean_std_boot))
}


MANOVA_CaiXia <- function (dat, K, N, p, alpha=0.05) {

  # CovO=array(0:0,dim=c(p,p))
  # for (j in 1:K){
  #   z = dat[[j]]
  #   CovO=CovO+cov(z)*(N[j]-1)
  # }
  # A = CovO/(sum(N)-K)

  # t1=array(0:0:0,dim=c(K,p,p))
  # for (j in 1:K){
  #   e = array(1:1,dim=c(1,N[j]))
  #   ayy = dat[[j]]
  #   ayybar = (e%*%ayy)/N[j];
  #   ayycov = cov(ayy)*(N[j]-1)/N[j];
  #   ayyy = ayy-((t(e))%*%ayybar);
  #   c1 = (t(ayyy^2))%*%(ayyy^2)/N[j];
  #   c2 = ((t(ayyy))%*%(ayyy)/N[j])*(ayycov);
  #   c3 = ayycov^2;
  #   t1[j,,] = c1-2*c2+c3;
  # }
  # AA = sqrt(K*sum(N)*(A^2)/(t1[1,,]+t1[2,,]+t1[3,,]));
  #
  # lambda = rep(0,50)
  # dif = rep(0,50)
  # for (j in 1:50){
  #   lambda[j] = 2*j/50*sqrt(log(p))
  #   dif[j] = norm(A*(abs(AA)>=lambda[j])-A,"F")
  # }
  # thr = max(lambda[which(dif==min(dif))])
  # O1hat = solve(A*(abs(AA)>=thr*10))

  # dat0 <- matrix(0, sum(N), p)
  # dat0[1:N[1], ] <- dat[[1]]
  # for (j in 2:K) {
  #   dat0[(sum(N[1:(j-1)])+1):(sum(N[1:j])), ] <- dat[[j]]
  # }

  # O1hat = sugm(dat0, lambda = NULL, nlambda = 10, lambda.min.ratio = NULL,
  #              rho = NULL, method = "clime", sym = "or", shrink=NULL,
  #              prec = 1e-4, max.ite = 1e4, standardize = FALSE,
  #              perturb = TRUE, verbose = TRUE)
  O1hat = diag(1, p, p) #O1hat$icov[[10]]

  estx = array(0:0,dim=c(K,p))
  Cov = 0
  for (j in 1:K){
    estz = t(O1hat %*% t(dat[[j]]))
    estx[j,] = colMeans(estz)
    Cov = Cov+cov(estz)*(N[j]-1)
  }
  Cov = Cov/(sum(N)-K)


  est = 0
  for (l in 1:(K-1)){
    for (t in (l+1):K){
      est = est + (N[l]*N[t]/(N[l]+N[t]))*(estx[l,]-estx[t,])^2/diag(Cov)
    }
  }
  test_stat_clime <- max(est)

  Cov_Ip = matrix(0, p, p)
  for (j in 1:K){
    Cov_Ip = Cov_Ip + cov(dat[[j]])*(N[j]-1)
  }
  Cov_Ip = Cov_Ip/(sum(N)-K)

  mean_list <- vector("list", K)
  for (k in 1:K) {
    mean_list[[k]] <- apply(dat[[k]], 2, mean)
  }

  est_Ip = 0
  for (l in 1:(K-1)){
    for (t in (l+1):K){
      est_Ip = est_Ip + (N[l]*N[t]/(N[l]+N[t]))*(mean_list[[l]]-mean_list[[t]])^2/diag(Cov_Ip)
    }
  }
  test_stat_Ip <- max(est_Ip)

  crinew = -K*log(-gamma((K-1)/2)*log(1-alpha))
  crinew = K*log(p)+K*(K-3)/2*log(log(p))+crinew

  return(c(test_stat_clime, test_stat_Ip, crinew,
           test_stat_clime>crinew, test_stat_Ip>crinew))
}
