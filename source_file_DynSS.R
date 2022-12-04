
## source file for DynSS method

################################################################################ 
################################################################################

rm(list = ls())

library(dplyr); library(data.table); library(spTimer); library(LaplacesDemon); 
library(extraDistr); library(pbapply); library(MASS); 
library(ggplot2); library(grid); library(gridExtra); library(ggridges); 

## some source functions 

gen_data_arm <- function(n=c(300),rand_pr=NULL,n_arm=NULL){
  ##
  ## n is the sample in one interim 
  ##
  dat <- data.table::CJ(
    id = 1:n,
    id_arm = NA_integer_
  )
  ##
  if(is.null(rand_pr)){
    stop("provide value for rand_pr")
  }
  ##
  if(is.null(n_arm)){
    n_arm <- length(rand_pr)
    dat$id_arm <- sample(1:n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  }
  else{
    dat$id_arm <- sample(n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  }
  ##  
  dat
}
##
gen_data_binodist <- function(N, theta=c(0.5,0.5,0.25), rand_pr=c(0.33,0.33,0.33)){
  ##
  ## theta = prob for all arms (control and treatment arms)
  ## N = total (including the control and treatment arms)
  ##
  input <- gen_data_arm(n=N, rand_pr = rand_pr, n_arm = NULL)
  out <- data.frame(input, y = NA)
  ##
  ss <- table(input$id_arm)
  arm <- unique(out$id_arm)
  for(i in arm){
    out[out$id_arm%in%i,"y"] <- rbinom(ss[names(ss)%in%i],size=1,prob=theta[i])
  }
  ##
  out
}
##
trial_binodist <- function(y = NULL, beta_prior=c(0.5,0.5), # beta(a,b)
                           sim=FALSE, nRep=5000) 
{
  ##
  ## beta-binomial distribution
  a = (beta_prior[1]+y[2]); b = (beta_prior[2]+y[1]) 
  ##
  if(isTRUE(sim)){
    ## sims 
    tmp <- rbeta(nRep, a, b)
    return(list(mc=tmp))
  }
  else{
    ## math expression 
    ## mean and sd of beta dist 
    post_theta_mn <- a/(a+b); post_theta_sd <- sqrt((a*b)/((a+b)^2*(a+b+1)))
    return(list(post_theta_mn=post_theta_mn,post_theta_sd=post_theta_sd))
  }
  ##
}
##
EDSS_binodist_opt <- function(dat=NULL, rand_assignment=c(1,1),
                              nRep=5000, true_theta=c(0.5,0.5),
                              prior_interest=list(c(1,4),c(1,1)), # first one is control, use beta-dist(alpha,beta)
                              prior_reference=list(c(1,1),c(1,1)),  # first one is control, use beta-dist(alpha,beta) 
                              n_grid=seq(100,200,1),
                              EDSS_cost="kld", use_arm=1){
  ##
  ## prior_interest = informed prior 
  ## prior_reference = reference/objective/baseline prior (may be vague)
  ## n_ref = sample size for the baseline/reference prior
  ## n_grid = grid sample size to check (bit different compared to the 2019 ECSS paper, where they used max(n_grid) = n_ref)
  ##  
  if(is.null(dat)){
    rand_pr <- rand_assignment/sum(rand_assignment)
    gen_data_binodist_fixed <- function(N, theta, rand_pr){
      ##
      ## theta = prob for all arms (control and treatment arms)
      ## N = total (including the control and treatment arms)
      ##
      input <- gen_data_arm(n=N, rand_pr = rand_pr, n_arm = NULL)
      input$id_arm <- rep(1:length(rand_pr),N/length(rand_pr))
      out <- data.frame(input, y = NA)
      ##
      ss <- table(input$id_arm)
      arm <- unique(out$id_arm)
      for(i in arm){
        out[out$id_arm%in%i,"y"] <- rbinom(ss[names(ss)%in%i],size=1,prob=theta[i])
      }
      ##
      out
    }
    dat <- gen_data_binodist_fixed(N=max(c(n_grid)), theta=true_theta, rand_pr=rand_pr)
  }
  else{
    dat <- dat
  }
  ##
  prior_reference <- as.matrix(data.frame(prior_reference))
  dimnames(prior_reference) <- NULL
  ##
  pp = use_arm # only for the control arm (or defined by use_arm)
  res_ref <- list()
  ddat <- table(dat$id_arm,dat$y)
  res_ref[[pp]] = trial_binodist(y = ddat[pp,], beta_prior=c(prior_reference[,pp]))
  ##
  bint_ref <- rnorm(nRep,res_ref[[pp]]$post_theta_mn,res_ref[[pp]]$post_theta_sd)
  b_mse <- mean((bint_ref-true_theta[1])^2)
  b_kld <- KLD(px=c(bint_ref),py=rep(true_theta[1],each=nRep))$mean.sum.KLD
  ##
  #store_valStat <- matrix(NA,length(n_grid),5)
  store_valStat <- matrix(NA,length(n_grid),7) 
  store_res <- list()
  prior_interest <- as.matrix(data.frame(prior_interest))
  dimnames(prior_interest) <- NULL
  for(j in 1:length(n_grid)){
    ddat <- table(dat[1:n_grid[j],]$id_arm,dat[1:n_grid[j],]$y)
    res <- list()
    res[[pp]] = trial_binodist(y = ddat[pp,], beta_prior=c(prior_interest[,pp]))
    bint <- rnorm(nRep,res[[1]]$post_theta_mn,res[[1]]$post_theta_sd)
    p_mse <- mean((bint-true_theta[1])^2)
    p_kld <- KLD(px=c(bint),py=rep(true_theta[1],each=nRep))$mean.sum.KLD
    ##
    tmp <- array(unlist(res),dim=c(2,1))
    dimnames(tmp)[[1]] <- c("mean","se"); dimnames(tmp)[[2]] <- paste0("arm",use_arm)
    store_res[[j]] <- tmp
    ##
    store_valStat[j,] <- c(sum(ddat[1,]),tmp[,1],p_mse,b_mse,p_kld,b_kld)
    ##
  }
  store_valStat <- data.frame(store_valStat)
  names(store_valStat) <- c("control_ss","mn","se","MSE","MSE_b","KL_mean","KL_mean_b") 
  store_valStat <- store_valStat %>% group_by(control_ss) %>% 
    summarize(mn=mean(mn),se=mean(se),MSE=mean(MSE),MSE_b=mean(MSE_b),
              KL_mean=mean(KL_mean),KL_mean_b=mean(KL_mean_b))
  store_valStat <- data.frame(store_valStat)
  ##
  store_valStat$kl_diff <- abs(store_valStat$KL_mean-store_valStat$KL_mean_b)
  store_valStat$mse_diff <- abs(store_valStat$MSE-store_valStat$MSE_b)
  ##
  if(EDSS_cost%in%"mse"){
    optmSS <- store_valStat[which(store_valStat$mse_diff==min(store_valStat$mse_diff)),]
  }
  else if(EDSS_cost%in%"kld"){
    optmSS <- store_valStat[which(store_valStat$kl_diff==min(store_valStat$kl_diff)),]
  }
  else{
    stop("define EDSS_cost correctly, can take arguments mse or kld ")
  }  
  ##
  return(list(optmSS=optmSS,store_valStat=store_valStat))
  ##
}
##
trail_interim_binoEDSS <- function(
    trial_type = "noninf", # takes "noninf" and "sup"
    margin_of_error = 0.1, # percentage 
    n_sample=c(100,50,50), 
    rand_assignment=c(1,1),
    nRep=5000, 
    prior_interest=list(c(1,1),c(1,1)), # first one is control, use beta-dist(alpha,beta)
    prior_reference=list(c(1,1),c(1,1)),  # first one is control, use beta-dist(alpha,beta) 
    method=c("EDSS"), # can take EDSS "ECSS","NONE" 
    prior_intrest_arm=c("ctr"), # can take "ctr" or c("ctr","trt")
    true_theta=c(0.9,0.5), # first one is the success rate for control
    cost = "kld", # kld or mse
    drop_arm_type=c("futility"), # can take either "futility", "efficacy", "none 
    D0=0.10, D1=0.90, sim=FALSE, mixture=TRUE) 
{
  ##
  ## non-inferiority trial: H1: theta_trt > theta_crt - margin_of_error => (theta_trt - theta_crt) > - margin_of_error
  ## superiority trial: H1: theta_trt > theta_crt => (theta_trt - theta_crt) > 0
  ##
  pr_change <- function(old, new){
    (new-old)/old
  }
  ##
  K <- length(true_theta)
  n_j <- length(n_sample)
  dat <- list()
  res <- list()
  optm <- list()
  optmSS_ctr <- c()
  store_pc <- c()
  n_dyn_rand_pr <- matrix(NA,K,n_j)
  pr_decision <- matrix(NA,K-1,n_j)
  n_interimSS <- matrix(NA,K,n_j)
  dyn_prior_interest <- list()
  drop_arm <- c()
  e = 0.0
  ##
  for(j in 1:n_j){
    ##
    if(j==1){
      ##
      rand_pr <- rand_assignment/sum(rand_assignment)
      n_dyn_rand_pr[,j] <- rand_pr
      dat[[j]] = gen_data_binodist(N=n_sample[j], theta=true_theta, rand_pr=rand_pr)
      ddat <- table(dat[[j]]$id_arm,dat[[j]]$y)
      ##
      ## for conventional model
      if(method%in%"NONE"){
        ##
        res_tmp <- matrix(NA,K,2)
        res_sim <- list()
        ##
        for(arm in 1:K){
          if(isTRUE(sim)){
            res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]]),sim=TRUE)$mc
            res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
          }
          else{
            res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]])))
          }
        }
        dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
        dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
        res[[j]] = res_tmp
        dyn_rand_pr <- rand_pr
        dyn_prior_interest[[j]] <- as.list(prior_interest)
        dyn_prior_interest[[j+1]] <- as.list(prior_interest)
        ##
      }
      ##
      else{
        ##
        if(method%in%"EDSS"){
          ##
          n_grid <- seq(round(nrow(dat[[j]])/2),nrow(dat[[j]]),1)
          optm[[j]] = EDSS_binodist_opt(dat=dat[[j]], 
                                        nRep=nRep, true_theta=true_theta,
                                        prior_interest=prior_interest, # first one is control, use beta-dist(alpha,beta)
                                        prior_reference=prior_reference,  # first one is control, use beta-dist(alpha,beta) 
                                        n_grid=n_grid,
                                        EDSS_cost=cost)$optmSS
          optmSS_ctr[j] <- unlist(optm[[j]][1,1])
          pc <- pr_change(old=rowSums(ddat)[1],new=optmSS_ctr[j])
          store_pc[j] <- unlist(pc)
          res_tmp <- matrix(NA,K,2)
          res_sim <- list()
          ##
          if(isTRUE(mixture)){
            if(pc < (-e) & method%in%"EDSS"){
              ##
              for(arm in 1:K){
                if(isTRUE(sim)){
                  res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]]),sim=TRUE)$mc
                  res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
                }
                else{
                  res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]])))
                }
              }
              dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
              dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
              res[[j]] = res_tmp
              tmp <- rand_pr[1]*pc+rand_pr[1]
              dyn_rand_pr <- unlist(c(tmp,rep((1-tmp)/(K-1),K-1))) # equal allocation of randomisation prob to treatment arms
              dyn_prior_interest[[j]] <- as.list(prior_interest)
              ##
              ## convert the bino-probs into beta shape parameters: concept of betaReg 
              ## a = x+1 and b = n-x+1
              ##
              dddat <- matrix(NA,K,2)
              for(arm in 1:K){
                tmp = gen_data_binodist(N=rowSums(ddat)[arm], theta=res_tmp[arm,1], rand_pr=1) # using postmean
                dddat[arm,] = table(tmp$id_arm,tmp$y)
              }
              a = dddat[,2]+1; b = dddat[,1]+1; ab = t(cbind(a,b)); ab = data.frame(ab)
              names(ab) <- paste0("arm_",1:K); 
              dyn_prior_interest[[j+1]] <- as.list(ab)
              ##
            }
            else{
              ##
              res_tmp <- matrix(NA,K,2)
              res_sim <- list()
              ##
              for(arm in 1:K){
                if(isTRUE(sim)){
                  res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_reference[[arm]]),sim=TRUE)$mc
                  res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
                }
                else{
                  res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_reference[[arm]])))
                }
              }
              dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
              dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
              res[[j]] = res_tmp
              dyn_rand_pr <- rand_pr
              dyn_prior_interest[[j]] <- as.list(prior_interest)
              dyn_prior_interest[[j+1]] <- as.list(prior_interest)
              ##
            }
          }
          else{
            ##
            res_tmp <- matrix(NA,K,2)
            res_sim <- list()
            ##
            for(arm in 1:K){
              if(isTRUE(sim)){
                res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]]),sim=TRUE)$mc
                res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
              }
              else{
                res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]])))
              }
            }
            dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
            dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
            res[[j]] = res_tmp
            tmp <- rand_pr[1]*pc+rand_pr[1]
            dyn_rand_pr <- unlist(c(tmp,rep((1-tmp)/(K-1),K-1))) # equal allocation of randomisation prob to treatment arms
            dyn_prior_interest[[j]] <- as.list(prior_interest)
            ##
            ## convert the bino-probs into beta shape parameters: concept of betaReg 
            ## a = x+1 and b = n-x+1
            ##
            dddat <- matrix(NA,K,2)
            for(arm in 1:K){
              tmp = gen_data_binodist(N=rowSums(ddat)[arm], theta=res_tmp[arm,1], rand_pr=1) # using postmean
              dddat[arm,] = table(tmp$id_arm,tmp$y)
            }
            a = dddat[,2]+1; b = dddat[,1]+1; ab = t(cbind(a,b)); ab = data.frame(ab)
            names(ab) <- paste0("arm_",1:K); 
            dyn_prior_interest[[j+1]] <- as.list(ab)
            ##
          }
          ##
        }
        else if(method%in%"ECSS"){
          ##
          res_tmp <- matrix(NA,K,2)
          res_sim <- list()
          ##
          for(arm in 1:K){
            if(isTRUE(sim)){
              res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]]),sim=TRUE)$mc
              res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
            }
            else{
              res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]])))
            }
          }
          dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
          dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
          res[[j]] = res_tmp
          tmp <- rand_pr[1]*pc+rand_pr[1]
          dyn_rand_pr <- unlist(c(tmp,rep((1-tmp)/(K-1),K-1)))
          dyn_prior_interest[[j]] <- as.list(prior_interest)
          dyn_prior_interest[[j+1]] <- as.list(prior_interest)
        }
        else{
          ##
          stop("Define argument method correctly, can take EDSS, ECSS and NONE")
          ##
        }
        ##
      }
      ##
      ## superiority trial: H1: theta_trt > theta_crt
      ##                    => (theta_crt - theta_trt) < 0
      ##
      ## non-inferiority trial: H1: theta_trt > theta_crt - margin_of_error 
      ##                        => (theta_trt - theta_crt) > - margin_of_error
      ##                        => (theta_crt - theta_trt) < margin_of_error
      ##
      ## Approximation method: considering cov(x,y) = 0 we get Var(x-y) = var(x)+Var(y) - 2*0 
      ##
      if(trial_type%in%"noninf"){
        for(arm in 2:K){
          if(isTRUE(sim)){
            pr_decision[arm-1,j] <- mean(res_sim[[arm]] > res_sim[[1]] - margin_of_error)
          }
          else{
            pr_decision[arm-1,j] <- pnorm(margin_of_error, res_tmp[1,1]-res_tmp[arm,1], res_tmp[1,2]+res_tmp[arm,2]) # Pr(b < delta)
            #1-pnorm(0, mu, sd) # Pr(b > delta)
          }
        }
      }
      else if(trial_type%in%"sup"){
        for(arm in 2:K){
          if(isTRUE(sim)){
            pr_decision[arm-1,j] <- mean(res_sim[[arm]] > res_sim[[1]])
          }
          else{
            pr_decision[arm-1,j] <- pnorm(0, res_tmp[1,1]-res_tmp[arm,1], res_tmp[1,2]+res_tmp[arm,2]) # Pr(b < 0)
          }
        }
      }
      else{
        stop("define trial_type correctly")
      }
      ##
      n_interimSS[,j] <- table(dat[[j]]$id_arm)
      names(dyn_prior_interest[[j]]) = c(paste0("arm_",1:K))
      ##
      n_dyn_rand_pr[,j+1] <- dyn_rand_pr
      ##
    }
    ##
    else{
      ## For no sample size cap: if one arm is dropped then other arms 
      ## may receive more participants up to the nMax[sum(n_sample)] participants 
      if("futility"%in%drop_arm_type & "efficacy"%in%drop_arm_type){
        drop_arm <- c(drop_arm,c(which(pr_decision[,j-1] < D0) + 1,which(pr_decision[,j-1] > D1) + 1))
        dyn_rand_pr[drop_arm] <- 0
        dyn_rand_pr <- dyn_rand_pr/sum(dyn_rand_pr)
        n_dyn_rand_pr[,j] <- dyn_rand_pr
      }
      else if(("futility"%in%drop_arm_type) & !("efficacy"%in%drop_arm_type)){
        drop_arm <- c(drop_arm,c(which(pr_decision[,j-1] < D0) + 1))
        dyn_rand_pr[drop_arm] <- 0
        dyn_rand_pr <- dyn_rand_pr/sum(dyn_rand_pr)
        n_dyn_rand_pr[,j] <- dyn_rand_pr
      }
      else if(!("futility"%in%drop_arm_type) & ("efficacy"%in%drop_arm_type)){
        drop_arm <- c(drop_arm,c(which(pr_decision[,j-1] > D1) + 1))
        dyn_rand_pr[drop_arm] <- 0
        dyn_rand_pr <- dyn_rand_pr/sum(dyn_rand_pr)
        n_dyn_rand_pr[,j] <- dyn_rand_pr
      }
      else if("none"%in%drop_arm_type){
        n_dyn_rand_pr[,j] <- dyn_rand_pr
      }
      else{
        stop("define drop_arm_type correctly")
      }
      ##
      new_dat <- gen_data_binodist(N=n_sample[j], theta=true_theta, rand_pr=dyn_rand_pr)
      dat[[j]] = data.frame(mapply(c,dat[[j-1]],new_dat,SIMPLIFY=FALSE))
      dat[[j]]$id <- 1:nrow(dat[[j]])
      n_interimSS[,j] <- table(dat[[j]]$id_arm)
      ddat <- table(dat[[j]]$id_arm,dat[[j]]$y)
      ##
      if(method%in%"NONE"){
        ##
        optmSS_ctr[j] <- unlist(rowSums(ddat)[1])
        pc <- pr_change(old=rowSums(ddat)[1],new=optmSS_ctr[j])
        store_pc[j] <- unlist(pc)
        res_tmp <- matrix(NA,K,2)
        for(arm in 1:K){
          if(isTRUE(sim)){
            res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]]),sim=TRUE)$mc
            res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
          }
          else{
            res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]])))
          }
        }
        dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
        dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
        res[[j]] = res_tmp
        dyn_prior_interest[[j]] <- as.list(prior_interest)
        dyn_prior_interest[[j+1]] <- as.list(prior_interest)
        ##
      }
      ##
      else{
        ##
        if(method%in%"EDSS"){
          #if(!isTRUE(unique(unlist(prior_interest[[1]])%in%unlist(prior_reference[[1]])))){
          n_grid <- seq(nrow(dat[[j-1]])+1,nrow(dat[[j]]),1)
          optm[[j]] = EDSS_binodist_opt(dat=dat[[j]], #dat[[j]][(nrow(dat[[j-1]])+1):nrow(dat[[j]]),], #new_dat, # 
                                        nRep=nRep, true_theta=true_theta,
                                        prior_interest=prior_interest, # first one is control, use beta-dist(alpha,beta)
                                        prior_reference=prior_reference,  # first one is control, use beta-dist(alpha,beta) 
                                        n_grid=n_grid,
                                        EDSS_cost=cost)$optmSS
          optmSS_ctr[j] <- unlist(optm[[j]][1,1])
          #}
          #else{
          #  optmSS_ctr[j] <- unlist(rowSums(ddat)[1])
          #}
          pc <- pr_change(old=rowSums(ddat)[1],new=optmSS_ctr[j])
          store_pc[j] <- unlist(pc)
          res_tmp <- matrix(NA,K,2)
          ##
          if(isTRUE(mixture)){
            if(pc < (-e) & method%in%"EDSS"){
              ##
              for(arm in 1:K){
                if(isTRUE(sim)){
                  res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(dyn_prior_interest[[1]][[arm]]),sim=TRUE)$mc
                  res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
                }
                else{
                  res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(dyn_prior_interest[[1]][[arm]])))
                }
              }
              dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
              dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
              res[[j]] = res_tmp
              tmp <- dyn_rand_pr[1]*pc+dyn_rand_pr[1]
              dyn_rand_pr[1] <- tmp
              dyn_rand_pr <- dyn_rand_pr/sum(dyn_rand_pr)
              ##
              ## convert the bino-probs into beta shape parameters: concept of betaReg 
              ## a = x+1 and b = n-x+1
              ##
              dddat <- matrix(NA,K,2)
              for(arm in 1:K){
                tmp = gen_data_binodist(N=rowSums(ddat)[arm], theta=res_tmp[arm,1], rand_pr=1) # using postmean
                dddat[arm,] = table(tmp$id_arm,tmp$y)
              }
              a = dddat[,2]+1; b = dddat[,1]+1; ab = t(cbind(a,b)); ab = data.frame(ab)
              names(ab) <- paste0("arm_",1:K); 
              dyn_prior_interest[[j+1]] <- as.list(ab)
              ##
            }
            else{
              ##
              res_tmp <- matrix(NA,K,2)
              res_sim <- list()
              ##
              for(arm in 1:K){
                if(isTRUE(sim)){
                  res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_reference[[arm]]),sim=TRUE)$mc
                  res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
                }
                else{
                  res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_reference[[arm]])))
                }
              }
              dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
              dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
              res[[j]] = res_tmp
              dyn_rand_pr <- rand_pr
              dyn_prior_interest[[j]] <- as.list(prior_interest)
              dyn_prior_interest[[j+1]] <- as.list(prior_interest)
              ##
            }
          }
          else{
            ##
            for(arm in 1:K){
              if(isTRUE(sim)){
                res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(dyn_prior_interest[[1]][[arm]]),sim=TRUE)$mc
                res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
              }
              else{
                res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(dyn_prior_interest[[1]][[arm]])))
              }
            }
            dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
            dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
            res[[j]] = res_tmp
            tmp <- dyn_rand_pr[1]*pc+dyn_rand_pr[1]
            dyn_rand_pr[1] <- tmp
            dyn_rand_pr <- dyn_rand_pr/sum(dyn_rand_pr)
            ##
            ## convert the bino-probs into beta shape parameters: concept of betaReg 
            ## a = x+1 and b = n-x+1
            ##
            dddat <- matrix(NA,K,2)
            for(arm in 1:K){
              tmp = gen_data_binodist(N=rowSums(ddat)[arm], theta=res_tmp[arm,1], rand_pr=1) # using postmean
              dddat[arm,] = table(tmp$id_arm,tmp$y)
            }
            a = dddat[,2]+1; b = dddat[,1]+1; ab = t(cbind(a,b)); ab = data.frame(ab)
            names(ab) <- paste0("arm_",1:K); 
            dyn_prior_interest[[j+1]] <- as.list(ab)
            ##
          }
        }
        ##
        ## for ECSS
        else if(method%in%"ECSS"){
          ##
          for(arm in 1:K){
            if(isTRUE(sim)){
              res_sim[[arm]] = trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]]),sim=TRUE)$mc
              res_tmp[arm,] = c(mean(res_sim[[arm]]),sd(res_sim[[arm]]))
            }
            else{
              res_tmp[arm,] = unlist(trial_binodist(y = ddat[arm,], beta_prior=c(prior_interest[[arm]])))
            }
          }
          dimnames(res_tmp)[[2]] = c("post_theta_mn","post_theta_sd")
          dimnames(res_tmp)[[1]] = c(paste0("arm_",1:K))
          res[[j]] = res_tmp
          tmp <- dyn_rand_pr[1]*pc+dyn_rand_pr[1]
          dyn_rand_pr[1] <- tmp
          dyn_rand_pr <- dyn_rand_pr/sum(dyn_rand_pr)
          dyn_prior_interest[[j]] <- as.list(prior_interest)
          dyn_prior_interest[[j+1]] <- as.list(prior_interest)
          ##
        }
        ##
        else{
          ##
          stop("Define argument method correctly, can take EDSS, ECSS and NONE")
          ##
        }
        ##
      }
      ##
      ##
      ## non-inferiority trial: H1: theta_trt > theta_crt - margin_of_error 
      ##                        => (theta_trt - theta_crt) > - margin_of_error
      ##                        => (theta_crt - theta_trt) < margin_of_error
      ##
      ## Approximation method: considering cov(x,y) = 0 we get Var(x-y) = var(x)+Var(y) - 2*0 
      ##
      if(trial_type%in%"noninf"){
        for(arm in 2:K){
          if(isTRUE(sim)){
            pr_decision[arm-1,j] <- mean(res_sim[[arm]] > res_sim[[1]] - margin_of_error)
          }
          else{
            pr_decision[arm-1,j] <- pnorm(margin_of_error, res_tmp[1,1]-res_tmp[arm,1], res_tmp[1,2]+res_tmp[arm,2]) # Pr(b < delta)
          }
        }
      }
      else if(trial_type%in%"sup"){
        for(arm in 2:K){
          if(isTRUE(sim)){
            pr_decision[arm-1,j] <- mean(res_sim[[arm]] > res_sim[[1]])
          }
          else{
            pr_decision[arm-1,j] <- pnorm(0, res_tmp[1,1]-res_tmp[arm,1], res_tmp[1,2]+res_tmp[arm,2]) # Pr(b < 0)
          }
        }
      }
      else{
        stop("define trial_type correctly")
      }
      ##
    }
  }
  ##
  dimnames(n_dyn_rand_pr)[[1]] = c("control",paste0("treat",1:(length(true_theta)-1)))
  dimnames(n_dyn_rand_pr)[[2]] = c(paste0("interim",1:n_j))
  ##
  dimnames(pr_decision)[[1]] = c(paste0("treat",1:(length(true_theta)-1)))
  dimnames(pr_decision)[[2]] = c(paste0("interim",1:n_j))
  ##
  dimnames(n_interimSS)[[1]] = c("control",paste0("treat",1:(length(true_theta)-1)))
  dimnames(n_interimSS)[[2]] = c(paste0("interim",1:n_j))
  ##
  post_theta <- t(sapply(res, function(x) x))
  post_theta <- array(post_theta,dim=c(n_j,K,2))
  dimnames(post_theta)[[1]] = c(paste0("interim",1:n_j))
  dimnames(post_theta)[[2]] = c("control",paste0("treat",1:(length(true_theta)-1)))
  dimnames(post_theta)[[3]] = c("mean","sd")
  ##
  dyn_prior <- array(unlist(dyn_prior_interest),dim=c(2,K,n_j))
  dimnames(dyn_prior)[[2]] = c("control",paste0("treat",1:(length(true_theta)-1)))
  dimnames(dyn_prior)[[3]] = c(paste0("interim",1:n_j))
  ##
  return(list(pr_decision=pr_decision, post_theta=post_theta, 
              n_dyn_rand_pr=n_dyn_rand_pr, n_interimSS=n_interimSS, optmSS_ctr=optmSS_ctr,
              per_change_ctr=store_pc*100, dyn_prior=dyn_prior)) #, optm=optm, res=res, dat=dat))
  ##
}
##
trail_interim_binoEDSS_simulator <- function(nSim=10,
                                             trial_type = "noninf", # takes "noninf" and "sup"
                                             margin_of_error = 0.1, # percentage 
                                             n_sample=c(100,rep(100,2)), 
                                             rand_assignment=c(1,1,1),
                                             nRep=5000, 
                                             prior_interest=list(c(2,1),c(1,1),c(1,1)), # first one is control, use beta-dist(alpha,beta)
                                             prior_reference=list(c(1,1),c(1,1),c(1,1)),  # first one is control, use beta-dist(alpha,beta) 
                                             method=c("EDSS"), # can take EDSS "ECSS","NONE" 
                                             prior_intrest_arm=c("ctr"), # can take "ctr" or c("ctr","trt")
                                             true_theta=c(0.5,0.8,0.1), # first one is the success rate for control
                                             cost = "kld", # kld or mse
                                             drop_arm_type=c("none"),
                                             D0=0.10, D1=0.90, 
                                             plot=TRUE, sim=FALSE, mixture=FALSE){
  ##
  options(warn=-1)
  ##
  library(pbapply)
  res <- pblapply(1:nSim, 
                  function(x) trail_interim_binoEDSS(
                    trial_type = trial_type, # takes "noninf" and "sup"
                    margin_of_error = margin_of_error, # percentage 
                    n_sample = n_sample, 
                    rand_assignment = rand_assignment,
                    nRep = nRep, 
                    prior_interest = prior_interest, # first one is control, use beta-dist(alpha,beta)
                    prior_reference = prior_reference,  # first one is control, use beta-dist(alpha,beta) 
                    method=method, # can take EDSS "ECSS","NONE" 
                    prior_intrest_arm=prior_intrest_arm, # can take "ctr" or c("ctr","trt")
                    true_theta = true_theta, # first one is the success rate for control
                    cost = cost, # kld or mse
                    drop_arm_type = drop_arm_type,
                    D0 = D0, D1 = D1, sim=sim, mixture=mixture))
  ##
  K <- length(true_theta)
  n_j <- length(n_sample)
  ##
  post_theta <- sapply(res, function(x) x$post_theta)
  post_theta <- array(c(post_theta),dim=c(n_j,K,2,nSim))
  dimnames(post_theta)[[1]] = c(paste0("interim",1:n_j))
  dimnames(post_theta)[[2]] = c("control",paste0("treat",1:(K-1)))
  dimnames(post_theta)[[3]] = c("mean","sd")
  ##
  pr_decision_criterion <- sapply(res, function(x) x$pr_decision)
  pr_decision_criterion <- array(c(pr_decision_criterion),dim=c(K-1,n_j,nSim))
  dimnames(pr_decision_criterion)[[1]] <- c(paste0("treat",1:(K-1)))
  dimnames(pr_decision_criterion)[[2]] <- c(paste0("interim",1:n_j))
  ##
  n_dyn_rand_pr <- sapply(res, function(x) x$n_dyn_rand_pr)
  n_dyn_rand_pr <- array(c(n_dyn_rand_pr),dim=c(K,n_j,nSim))
  dimnames(n_dyn_rand_pr)[[1]] = c("control",paste0("treat",1:(K-1)))
  dimnames(n_dyn_rand_pr)[[2]] = c(paste0("interim",1:n_j))
  ##
  sample_size <- sapply(res, function(x) x$n_interimSS)
  sample_size <- array(c(sample_size),dim=c(K,n_j,nSim))
  dimnames(sample_size)[[1]] <- c("control",paste0("treat",1:(K-1)))
  dimnames(sample_size)[[2]] <- c(paste0("interim",1:n_j))
  ##
  optmSS_ctr <- sapply(res, function(x) x$optmSS_ctr)
  dimnames(optmSS_ctr)[[1]] <- c(paste0("interim",1:n_j))
  ##
  per_change_ctr <- sapply(res, function(x) x$per_change_ctr)
  dimnames(per_change_ctr)[[1]] <- c(paste0("interim",1:n_j))
  ##
  dyn_prior <- sapply(res, function(x) x$dyn_prior)
  dyn_prior <- apply(dyn_prior,1,mean)
  dyn_prior <- array(dyn_prior,dim=c(2,K,n_j))
  dimnames(dyn_prior)[[2]] = c("control",paste0("treat",1:(length(true_theta)-1)))
  dimnames(dyn_prior)[[3]] = c(paste0("interim",1:n_j))
  ##
  pr_decisionD0 <- NULL
  pr_decisionD1 <- NULL
  fnc_replace_to_one <- function(x){
    ## x is a vector 
    if(length(which(x==1))>0){
      x[which(x==1)[1]:length(x)] <- 1
    }
    x
  }
  for(i in 1:(K-1)){
    tmp <- (pr_decision_criterion[i,,]<D0)
    pr_decisionD0 <- rbind(pr_decisionD0,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
    tmp <- (pr_decision_criterion[i,,]>D1)
    pr_decisionD1 <- rbind(pr_decisionD1,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
  }  
  pr_decision <- cbind(rep(1:length(n_sample),(K-1)),c(t(pr_decisionD0)),c(t(pr_decisionD1)))
  dimnames(pr_decision)[[2]] <- c("interim",paste0("D0_",D0),paste0("D1_",D1))
  pr_decision <- data.frame(pr_decision)
  pr_decision$arm <- rep(2:K,each=length(n_sample))
  pr_decision$tot_sample_size <- cumsum(n_sample)
  pr_decision$dec_pr <- NA
  pr_decision$dec_pr_sd <- NA
  for(i in 1:(K-1)){
    pr_decision[pr_decision$arm%in%(i+1),"dec_pr"] <- apply(pr_decision_criterion[i,,],1,mean)
    pr_decision[pr_decision$arm%in%(i+1),"dec_pr_sd"] <- apply(pr_decision_criterion[i,,],1,sd)
  }
  ##
  dat <- data.table::CJ(
    interim = 1:length(n_sample),
    arm = 1:K,
    true_theta = NA_real_,
    est_theta_mn = NA_real_,
    est_theta_sd = NA_real_,
    est_dyn_rand_pr_mn = NA_real_,
    est_dyn_rand_pr_sd = NA_real_,
    est_sample_mn = NA_real_,
    est_sample_sd = NA_real_
  )
  for(i in 1:K){
    dat[dat$arm%in%i,"true_theta"] <- true_theta[i]
    dat[dat$arm%in%i,"est_theta_mn"] <- apply(post_theta[,i,1,],1,mean)
    dat[dat$arm%in%i,"est_theta_sd"] <- apply(post_theta[,i,2,],1,mean)
    dat[dat$arm%in%i,"est_dyn_rand_pr_mn"] <- apply(n_dyn_rand_pr[i,,],1,mean)
    dat[dat$arm%in%i,"est_dyn_rand_pr_sd"] <- apply(n_dyn_rand_pr[i,,],1,sd)
    dat[dat$arm%in%i,"est_sample_mn"] <- apply(sample_size[i,,],1,mean)
    dat[dat$arm%in%i,"est_sample_sd"] <- apply(sample_size[i,,],1,sd)
  }
  dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
  dat$tot_sample_size <- rep(cumsum(n_sample),each=K)
  dat$arm <- as.factor(dat$arm)
  pr_decision$arm <- as.factor(pr_decision$arm)
  names(pr_decision)[2:3] <- c("x1","x2") # ??
  ##
  if(isTRUE(plot)){
    ##
    plot_EDSS(dat=dat,pr_decision=pr_decision,true_theta=true_theta,
              method=method,trial_type=trial_type,D0=D0,D1=D1)
    ##
  }
  ##
  names(pr_decision)[2:3] <- c(paste0("D0_",D0),paste0("D1_",D1))
  ##
  return(list(stat=dat,dyn_prior=dyn_prior,ref_prior=prior_reference,
              pr_decision=pr_decision, method=method, trial_type=trial_type,
              D0=D0, D1=D1))
  ##
}
##
EDSS_simulation <- function(nSim=10, dat=NULL, rand_assignment=c(1,1), nRep=5000,
                            true_theta=c(0.5,0.5),
                            prior_interest=list(c(10,10),c(1,1)), # list function for all arms 
                            prior_reference=list(c(1,1),c(1,1)), # list function for all arms 
                            n_grid=seq(100,200,1), cost ="kld", use_arm=1,
                            plot=TRUE){
  ##
  library(pbapply)
  out <- pblapply(1:nSim, 
                  function(x) EDSS_binodist_opt(dat=dat,
                                                rand_assignment=rand_assignment,nRep=nRep,
                                                true_theta=true_theta,prior_interest=prior_interest, # first one is control, use beta-dist(alpha,beta)
                                                prior_reference=prior_reference,  # first one is control, use beta-dist(alpha,beta) 
                                                n_grid=n_grid, EDSS_cost=cost, use_arm=use_arm)$store_valStat) #[,c(1,6)])
  ##
  valStat <- sapply(out, function(x) as.matrix(x[]))
  mn <- rowMeans(valStat)
  valStat <- data.frame(matrix(c(mn),dim(out[[1]])[1],dim(out[[1]])[2]))
  names(valStat) <- dimnames(out[[1]])[[2]]
  ##
  if(isTRUE(plot)){
    library(ggplot2)
    if(cost=="kld"){
      pp = ggplot(valStat,aes(x=control_ss,y=kl_diff)) + geom_smooth() +
        labs(x="Sample size",y="Kullback–Leibler divergence (KLD) differences")
      print(pp)
    }
    else if(cost=="mse"){
      pp = ggplot(valStat,aes(x=control_ss,y=mse_diff)) + geom_smooth() +
        labs(x="Sample size",y="Mean squared error (MSE) differences")
      print(pp)
    }
    else{
      print("Define cost function using kld or mse\n")
    }
  }
  ##
  return(list(valStat=valStat))
  ##
}
##
plot_EDSS <- function(dat,pr_decision,true_theta,method,trial_type,D0,D1){
  ##
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  #
  p1 <- ggplot(dat,aes(x=tot_sample_size,y=est_theta_mn,color=arm,shape=arm)) +
    geom_point(size=2.5,position=position_dodge(0.1)) + geom_hline(yintercept = true_theta,lty=2) +
    geom_errorbar(aes(ymin=est_theta_mn-1.96*est_theta_sd,ymax=est_theta_mn+1.96*est_theta_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1.1) +
    labs(y="Estimated success",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  ##
  p2 <- ggplot(dat,aes(x=tot_sample_size,y=est_sample_mn,color=arm,shape=arm)) +
    geom_point(size=2,position=position_dodge(0.1)) + 
    geom_errorbar(aes(ymin=est_sample_mn-1.96*est_sample_sd,ymax=est_sample_mn+1.96*est_sample_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1) +
    labs(y="Estimated sample size",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  ##
  df = dat[,c("interim","arm","est_dyn_rand_pr_mn","est_dyn_rand_pr_sd","tot_sample_size")]
  df$id = 1:nrow(df)
  myrnorm <- function(mn,sd,n){
    rnorm(n,mn,sd)
  }
  tmp = matrix(unlist(tapply(df$est_dyn_rand_pr_mn,list(df$id),myrnorm,df$est_dyn_rand_pr_sd,100)),100,nrow(df))
  df = df[rep(seq_len(nrow(df)), each = 100), ]
  df$est_dyn_rand_pr_mn = c(tmp)
  p3 <- ggplot() +
    geom_point(data=dat,aes(x=tot_sample_size,y=est_dyn_rand_pr_mn,color=arm,shape=arm),size=2) +
    geom_smooth(data=df,aes(x=tot_sample_size,y=est_dyn_rand_pr_mn,color=arm)) + ylim(0,1) +
    labs(y="Randomisation probability",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette) 
  ##
  if(trial_type%in%"noninf"){
    yylab <- "Pr(theta > theta_s - delta)"
  }
  else if(trial_type%in%"sup"){
    yylab <- "Pr(theta > theta_s)"
  }
  else{
    yylab <- ""
  }
  p4 <- ggplot(pr_decision,aes(x=tot_sample_size,y=dec_pr,color=arm,shape=arm)) +
    geom_point(size=2) + geom_smooth() + ylim(0,1) +
    labs(y=paste0(yylab,"\n D_0: ",D0,", D_1: ",D1),x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette[-1])
  ##    
  p5 <- ggplot(pr_decision,aes(x=tot_sample_size,y=x1,color=arm,shape=arm)) + 
    geom_point(size=2) + geom_smooth() + ylim(0,1) + 
    labs(y=paste0("Pr(Dropout for futility)\n D_0: ",D0),x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette[-1])
  ##
  p6 <- ggplot(pr_decision,aes(x=tot_sample_size,y=x2,color=arm,shape=arm)) + 
    geom_point(size=2) + geom_smooth() + ylim(0,1) +
    labs(y=paste0("Pr(Staying in the trial)\n D_1: ",D1),x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette[-1])
  ##
  grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,
               top = textGrob(paste(method),gp=gpar(fontsize=15,font=3)))
  ##
}
##
prior_post_plots <- function(x, nRep=5000){
  ##
  dyn_prior <- x$dyn_prior # beta(a,b) x arms x interim
  post <- x$stat[,c("interim","arm","est_theta_mn","est_theta_sd","true_theta")]
  d <- data.table::CJ(
    interim = 1:dim(dyn_prior)[3],
    arm = 1:dim(dyn_prior)[2],
    nRep = 1:nRep,
    prior = NA_real_,
    post = NA_real_,
    true = NA_real_,
    ref_prior = NA_real_
  )
  for(i in 1:dim(dyn_prior)[3]){ # interim
    for(j in 1:dim(dyn_prior)[2]){ # arm
      d[d$interim==i & d$arm==j,"prior"] = rbeta(nRep,dyn_prior[1,j,i],dyn_prior[2,j,i])
      d[d$interim==i & d$arm==j,"post"] = rnorm(nRep,unlist(post[post$interim==i & post$arm==j,"est_theta_mn"]),
                                                unlist(post[post$interim==i & post$arm==j,"est_theta_sd"]))
      d[d$interim==i & d$arm==j,"true"] = post[post$interim==i & post$arm==j,"true_theta"]
    }
  }
  d$interim <- as.factor(d$interim); d$arm <- as.factor(d$arm)
  levels(d$interim) <- c(paste0("Interim",1:dim(dyn_prior)[3]))
  levels(d$arm) <- c("control",paste0("treat",1:(dim(dyn_prior)[2]-1)))
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  p1 <- ggplot(d, aes(x = prior, y = arm, group = arm, fill=arm)) + 
    geom_density_ridges(rel_min_height = 0.01) + xlim(0,1) +
    facet_wrap(~interim) + labs(x="Dynamic prior distribution",y="Trial arms") +
    scale_fill_manual(values=safe_colorblind_palette) +
    theme(legend.position = "none")
  p2 <- ggplot(d, aes(x = post, y = arm, group = arm, fill=arm)) + 
    geom_density_ridges(rel_min_height = 0.01) + xlim(0,1) +
    geom_vline(aes(xintercept=true),d,linetype="dotdash") +
    facet_wrap(~interim) + labs(x="Posterior distribution",y="Trial arms") +
    scale_fill_manual(values=safe_colorblind_palette) +
    theme(legend.position = "none")
  grid.arrange(p1,p2)
}
##
plot_dynRand <- function(dat){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  df = dat[,c("interim","arm","est_dyn_rand_pr_mn","est_dyn_rand_pr_sd","tot_sample_size")]
  df$id = 1:nrow(df)
  myrnorm <- function(mn,sd,n){
    rnorm(n,mn,sd)
  }
  tmp = matrix(unlist(tapply(df$est_dyn_rand_pr_mn,list(df$id),myrnorm,df$est_dyn_rand_pr_sd,100)),100,nrow(df))
  df = df[rep(seq_len(nrow(df)), each = 100), ]
  df$est_dyn_rand_pr_mn = c(tmp)
  p <- ggplot() +
    geom_point(data=dat,aes(x=tot_sample_size,y=est_dyn_rand_pr_mn,color=arm,shape=arm),size=2) +
    geom_smooth(data=df,aes(x=tot_sample_size,y=est_dyn_rand_pr_mn,color=arm)) + ylim(0,1) +
    labs(y="Randomisation probability",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette) 
  p
}
##
plot_dynSS <- function(dat){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  p <- ggplot(dat,aes(x=tot_sample_size,y=est_sample_mn,color=arm,shape=arm)) +
    geom_point(size=2,position=position_dodge(0.1)) + 
    geom_errorbar(aes(ymin=est_sample_mn-1.96*est_sample_sd,ymax=est_sample_mn+1.96*est_sample_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1) +
    labs(y="Estimated sample size",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p
}
##
plot_dynSSRand <- function(dat,true_theta){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  df = dat[,c("interim","arm","est_dyn_rand_pr_mn","est_dyn_rand_pr_sd","tot_sample_size")]
  df$id = 1:nrow(df)
  myrnorm <- function(mn,sd,n){
    rnorm(n,mn,sd)
  }
  tmp = matrix(unlist(tapply(df$est_dyn_rand_pr_mn,list(df$id),myrnorm,df$est_dyn_rand_pr_sd,100)),100,nrow(df))
  df = df[rep(seq_len(nrow(df)), each = 100), ]
  df$est_dyn_rand_pr_mn = c(tmp)
  p1 <- ggplot(dat,aes(x=tot_sample_size,y=est_theta_mn,color=arm,shape=arm)) +
    geom_point(size=2.5,position=position_dodge(0.1)) + geom_hline(yintercept = true_theta,lty=2) +
    geom_errorbar(aes(ymin=est_theta_mn-1.96*est_theta_sd,ymax=est_theta_mn+1.96*est_theta_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1.1) +
    labs(y="Estimated success",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p2 <- ggplot() +
    geom_point(data=dat,aes(x=tot_sample_size,y=est_dyn_rand_pr_mn,color=arm,shape=arm),size=2) +
    geom_smooth(data=df,aes(x=tot_sample_size,y=est_dyn_rand_pr_mn,color=arm)) + ylim(0,1) +
    labs(y="Randomisation probability",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette) 
  p3 <- ggplot(dat,aes(x=tot_sample_size,y=est_sample_mn,color=arm,shape=arm)) +
    geom_point(size=2,position=position_dodge(0.1)) + 
    geom_errorbar(aes(ymin=est_sample_mn-1.96*est_sample_sd,ymax=est_sample_mn+1.96*est_sample_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1) +
    labs(y="Estimated sample size",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  grid.arrange(p1,p2,p3,nrow=1,ncol=3)
}
##
fig_1 <- function(){
  out1 <- EDSS_simulation(nSim=50,prior_interest=list(c(1,1),c(1,1)),plot=FALSE)
  out2 <- EDSS_simulation(nSim=50,prior_interest=list(c(5,5),c(1,1)),plot=FALSE)
  out3 <- EDSS_simulation(nSim=50,prior_interest=list(c(10,10),c(1,1)),plot=FALSE)
  out <- rbind(data.frame(out1$valStat,prior="Beta(1,1)"),
               data.frame(out2$valStat,prior="Beta(5,5)"),
               data.frame(out3$valStat,prior="Beta(10,10)"))
  out$prior <- as.factor(out$prior)
  out$prior <- ordered(out$prior, levels =c("Beta(1,1)", "Beta(5,5)", "Beta(10,10)"))
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  p = ggplot(out,aes(x=control_ss,y=kl_diff,colour=prior)) + geom_smooth() + 
    labs(x="Sample size",y="Kullback-Leibler divergence (KLD) differences", colour="") +
    theme(legend.position="bottom") +
    scale_color_manual(values=safe_colorblind_palette)
  print(p)
}
##
fig_2 <- function(){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  true_theta = c(0.5,0.5)
  out1 <- trail_interim_binoEDSS_simulator(nSim=50, method="NONE", n_sample=c(100,100,100),
                                           rand_assignment=c(1,1),
                                           true_theta = true_theta,
                                           prior_interest=list(c(1,1),c(1,1)), 
                                           prior_reference=list(c(1,1),c(1,1)),plot=FALSE)
  p1 <- ggplot(out1$stat,aes(x=tot_sample_size,y=est_sample_mn,color=arm,shape=arm)) +
    geom_point(size=2,position=position_dodge(0.1)) + 
    geom_errorbar(aes(ymin=est_sample_mn-1.96*est_sample_sd,ymax=est_sample_mn+1.96*est_sample_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1) +
    labs(y="Estimated sample size",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p1 <- p1 + annotate(geom="text",x=115, y=100, label="Beta(1,1)")
  p11 <- ggplot(out1$stat,aes(x=tot_sample_size,y=est_theta_mn,color=arm,shape=arm)) +
    geom_point(size=2.5,position=position_dodge(0.1)) + geom_hline(yintercept = true_theta,lty=2) +
    geom_errorbar(aes(ymin=est_theta_mn-1.96*est_theta_sd,ymax=est_theta_mn+1.96*est_theta_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1.1) +
    labs(y="Estimated success",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p11 <- p11 + annotate(geom="text",x=125, y=0.6, label="Beta(1,1)")
  out2 <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                           rand_assignment=c(1,1),
                                           true_theta = true_theta,
                                           prior_interest=list(c(26,26),c(1,1)), 
                                           prior_reference=list(c(1,1),c(1,1)),plot=FALSE)
  p2 <- ggplot(out2$stat,aes(x=tot_sample_size,y=est_sample_mn,color=arm,shape=arm)) +
    geom_point(size=2,position=position_dodge(0.1)) + 
    geom_errorbar(aes(ymin=est_sample_mn-1.96*est_sample_sd,ymax=est_sample_mn+1.96*est_sample_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1) +
    labs(y="Estimated sample size",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p2 <- p2 + annotate(geom="text",x=115, y=100, label="Beta(26,26)")
  p22 <- ggplot(out2$stat,aes(x=tot_sample_size,y=est_theta_mn,color=arm,shape=arm)) +
    geom_point(size=2.5,position=position_dodge(0.1)) + geom_hline(yintercept = true_theta,lty=2) +
    geom_errorbar(aes(ymin=est_theta_mn-1.96*est_theta_sd,ymax=est_theta_mn+1.96*est_theta_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1.1) +
    labs(y="Estimated success",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p22 <- p22 + annotate(geom="text",x=125, y=0.6, label="Beta(26,26)")
  out3 <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                           rand_assignment=c(1,1),
                                           true_theta = true_theta,
                                           prior_interest=list(c(42,10),c(1,1)), 
                                           prior_reference=list(c(1,1),c(1,1)),plot=FALSE)
  p3 <- ggplot(out3$stat,aes(x=tot_sample_size,y=est_sample_mn,color=arm,shape=arm)) +
    geom_point(size=2,position=position_dodge(0.1)) + 
    geom_errorbar(aes(ymin=est_sample_mn-1.96*est_sample_sd,ymax=est_sample_mn+1.96*est_sample_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1) +
    labs(y="Estimated sample size",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p3 <- p3 + annotate(geom="text",x=115, y=100, label="Beta(42,10)")
  p33 <- ggplot(out3$stat,aes(x=tot_sample_size,y=est_theta_mn,color=arm,shape=arm)) +
    geom_point(size=2.5,position=position_dodge(0.1)) + geom_hline(yintercept = true_theta,lty=2) +
    geom_errorbar(aes(ymin=est_theta_mn-1.96*est_theta_sd,ymax=est_theta_mn+1.96*est_theta_sd),
                  colour="black", width=2, position=position_dodge(0.1)) + geom_line(size=1.1) +
    labs(y="Estimated success",x="Total number of enrolments",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p33 <- p33 + annotate(geom="text",x=125, y=0.65, label="Beta(42,10)")
  grid.arrange(p1,p11,p2,p22,p3,p33,nrow=3,ncol=2)
}
##
fig_3 <- function(){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  true_theta = c(0.5,0.5)
  out_r1 <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                             rand_assignment=c(1,1),
                                             true_theta = true_theta,
                                             prior_interest=list(c(26,26),c(1,1)), 
                                             prior_reference=list(c(1,1),c(1,1)),plot=FALSE)
  p_r1 <- ggplot(out_r1$stat,aes(x=interim,y=est_dyn_rand_pr_mn,color=arm,shape=arm)) +
    geom_point(size=2) + geom_smooth() + ylim(0,1) + 
    labs(y="Randomisation probability",x="Stage (t)",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p_r1 <- p_r1 + annotate(geom="text",x=1.5, y=0.25, label="Beta(26,26)") 
  ##
  out_r2 <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                             rand_assignment=c(1,1),
                                             true_theta = true_theta,
                                             prior_interest=list(c(42,10),c(1,1)), 
                                             prior_reference=list(c(1,1),c(1,1)),plot=FALSE)
  p_r2 <- ggplot(out_r2$stat,aes(x=interim,y=est_dyn_rand_pr_mn,color=arm,shape=arm)) +
    geom_point(size=2) + geom_smooth() + ylim(0,1) +
    labs(y="Randomisation probability",x="Stage (t)",color="Arm",shape="Arm") +
    scale_color_manual(values=safe_colorblind_palette)
  p_r2 <- p_r2 + annotate(geom="text",x=1.5, y=0.25, label="Beta(42,10)")
  grid.arrange(p_r1,p_r2,nrow=1,ncol=2)
}
##
fig_4 <- function(D=list(c(0.05,0.95))){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  n_sample <- c(100,100,100)
  true_theta <- list(c(0.5,0.5),c(0.5,0.6),c(0.5,0.7),c(0.5,0.8),c(0.5,0.9))
  rand_assignment <- c(1,1)
  trial_type <- c("sup")
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method="NONE", 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(26,26),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE,mixture=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p1=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(superiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(26,26)")
  p2=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point() + ylim(49,153) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(26,26)")
  ##
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(26,26),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE,mixture=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p11=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(superiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(26,26)")
  p22=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point() + ylim(49,153) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(26,26)")
  ##
  #grid.arrange(p1,p2,p11,p22,nrow=2,ncol=2)
  grid.arrange(arrangeGrob(p11, top = 'Dynamic method'),arrangeGrob(p22, top = 'Dynamic method'),
               arrangeGrob(p1, top = 'Conventional method'),arrangeGrob(p2, top = 'Conventional method'),
               nrow=2,ncol=2)
}
##
fig_S3 <- function(D=list(c(0.05,0.99))){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  n_sample <- c(100,100,100)
  true_theta <- list(c(0.5,0.5),c(0.5,0.6),c(0.5,0.7),c(0.5,0.8),c(0.5,0.9))
  rand_assignment <- c(1,1)
  trial_type <- c("sup")
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method="NONE", 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(42,10),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE,mixture=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p1=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(superiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(42,10)")
  p2=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point() + ylim(49,153) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(42,10)")
  ##
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(42,10),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE,mixture=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p11=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(superiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(42,10)")
  p22=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point() + ylim(49,153) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(42,10)")
  ##
  #grid.arrange(p1,p2,p11,p22,nrow=2,ncol=2)
  grid.arrange(arrangeGrob(p11, top = 'Dynamic method'),arrangeGrob(p22, top = 'Dynamic method'),
               arrangeGrob(p1, top = 'Conventional method'),arrangeGrob(p2, top = 'Conventional method'),
               nrow=2,ncol=2)
}
##
fig_5 <- function(D=list(c(0.050,0.950)),margin_of_error=0.1){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  n_sample <- c(100,100,100)
  true_theta <- list(c(0.5,0.5),c(0.5,0.6),c(0.5,0.7),c(0.5,0.8),c(0.5,0.9))
  rand_assignment <- c(1,1)
  trial_type <- c("noninf")
  method <- c("EDSS")
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method=method, 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(26,26),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       margin_of_error = margin_of_error,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p1=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(non-inferiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(26,26)")
  p2=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+ylim(49,153) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(26,26)")
  ##
  trial_type <- c("noninf")
  method <- c("NONE")
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method=method, 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(26,26),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       margin_of_error = margin_of_error,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p11=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(non-inferiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(26,26)")
  p22=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+ylim(49,151) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(26,26)")
  #grid.arrange(p1,p2,p11,p22,nrow=2,ncol=2)
  grid.arrange(arrangeGrob(p1, top = 'Dynamic method'),arrangeGrob(p2, top = 'Dynamic method'),
               arrangeGrob(p11, top = 'Conventional method'),arrangeGrob(p22, top = 'Conventional method'),
               nrow=2,ncol=2)
  
}
##
fig_S4 <- function(D=list(c(0.050,0.950)),margin_of_error=0.05){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  n_sample <- c(100,100,100)
  true_theta <- list(c(0.5,0.5),c(0.5,0.6),c(0.5,0.7),c(0.5,0.8),c(0.5,0.9))
  rand_assignment <- c(1,1)
  trial_type <- c("noninf")
  method <- c("EDSS")
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method=method, 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(42,10),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       margin_of_error = margin_of_error,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p1=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(non-inferiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(42,10)")
  p2=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+ylim(47,153) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(42,10)")
  ##
  trial_type <- c("noninf")
  method <- c("NONE")
  out_sup_store <- list()
  store <- NULL
  for(i in 1:length(true_theta)){
    out_sup <- list()
    for(j in 1:length(D)){
      out_sup[[j]] <- trail_interim_binoEDSS_simulator(nSim=50, method=method, 
                                                       n_sample = n_sample,
                                                       rand_assignment = rand_assignment,
                                                       true_theta = true_theta[[i]],
                                                       prior_interest = list(c(42,10),c(1,1)), 
                                                       prior_reference = list(c(1,1),c(1,1)),
                                                       trial_type = trial_type,
                                                       margin_of_error = margin_of_error,
                                                       drop_arm_type=c("none"),
                                                       D0=D[[j]][1],D1=D[[j]][2], 
                                                       plot=FALSE)
      tmp <- cbind(out_sup[[j]]$stat,true_theta[[i]][2]-true_theta[[i]][1],
                   true_theta[[i]][1],true_theta[[i]][2],D[[j]][1],D[[j]][2])
      names(tmp)[10:11] <- c("D0","D1")
      store <- rbind(store, tmp)
    }
    out_sup_store[[i]] = out_sup
  }
  ##
  p11=ggplot(store[store$arm==2,],aes(x=interim,y=D1)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+
    labs(y="Pr(non-inferiority)",x="Stage (t)",color="Effect size") +
    annotate(geom="text",x=1.5, y=0.25, label="Beta(42,10)")
  p22=ggplot(store[store$arm==1,],aes(x=interim,y=est_sample_mn)) +
    geom_line(aes(color=as.factor(V2),shape=as.factor(V2))) +
    geom_point()+ylim(47,151) +
    labs(y="Sample size for control arm",x="Stage (t)",color="Effect size")+
    annotate(geom="text",x=1.5, y=125, label="Beta(42,10)")
  #grid.arrange(p1,p2,p11,p22,nrow=2,ncol=2)
  grid.arrange(arrangeGrob(p1, top = 'Dynamic method'),arrangeGrob(p2, top = 'Dynamic method'),
               arrangeGrob(p11, top = 'Conventional method'),arrangeGrob(p22, top = 'Conventional method'),
               nrow=2,ncol=2)
  
}
##
##
fig_6 <- function(drop_arm_type=c("futility"),D0=0.025,D1=0.975){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  true_theta <- c(0.5,0.5,0.3,0.7)
  out_f_sup <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                                rand_assignment=c(1,1,1,1),
                                                true_theta = true_theta,
                                                prior_interest=list(c(26,26),c(1,1),c(1,1),c(1,1)), 
                                                prior_reference=list(c(1,1),c(1,1),c(1,1),c(1,1)),
                                                trial_type = "sup",
                                                drop_arm_type=drop_arm_type,
                                                D0=D0,D1=D1,plot=FALSE)
  p1 = plot_dynSS(dat=out_f_sup$stat) + annotate(geom="text",x=115, y=80, label="Beta(26,26)")
  p2 = plot_dynRand(dat=out_f_sup$stat) + annotate(geom="text",x=115, y=0.6, label="Beta(26,26)")
  out_f_sup <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                                rand_assignment=c(1,1,1,1),
                                                true_theta = true_theta,
                                                prior_interest=list(c(42,10),c(1,1),c(1,1),c(1,1)), 
                                                prior_reference=list(c(1,1),c(1,1),c(1,1),c(1,1)),
                                                trial_type = "sup",
                                                drop_arm_type=drop_arm_type,
                                                D0=D0,D1=D1,plot=FALSE)
  p3 = plot_dynSS(dat=out_f_sup$stat) + annotate(geom="text",x=115, y=80, label="Beta(42,10)")
  p4 = plot_dynRand(dat=out_f_sup$stat) + annotate(geom="text",x=115, y=0.6, label="Beta(42,10)")
  #grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,top = textGrob("Superiority trial (arm dropping for futility)",gp=gpar(fontsize=12,font=3)))
  grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
}
##
##
fig_7 <- function(drop_arm_type=c("futility"),D0=0.025,D1=0.975,margin_of_error=0.10){
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  true_theta <- c(0.5,0.5,0.3,0.7)
  out_f_noninf <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                                   rand_assignment=c(1,1,1,1),
                                                   true_theta = true_theta,
                                                   prior_interest=list(c(26,26),c(1,1),c(1,1),c(1,1)), 
                                                   prior_reference=list(c(1,1),c(1,1),c(1,1),c(1,1)),
                                                   trial_type = "noninf",
                                                   margin_of_error = margin_of_error,
                                                   drop_arm_type=drop_arm_type,
                                                   D0=D0,D1=D1,plot=FALSE)
  p1 = plot_dynSS(dat=out_f_noninf$stat) + annotate(geom="text",x=115, y=80, label="Beta(26,26)")
  p2 = plot_dynRand(dat=out_f_noninf$stat) + annotate(geom="text",x=115, y=0.6, label="Beta(26,26)")
  out_f_noninf <- trail_interim_binoEDSS_simulator(nSim=50, method="EDSS", n_sample=c(100,100,100),
                                                   rand_assignment=c(1,1,1,1),
                                                   true_theta = true_theta,
                                                   prior_interest=list(c(42,10),c(1,1),c(1,1),c(1,1)), 
                                                   prior_reference=list(c(1,1),c(1,1),c(1,1),c(1,1)),
                                                   trial_type = "noninf",
                                                   margin_of_error = margin_of_error,
                                                   drop_arm_type=drop_arm_type,
                                                   D0=D0,D1=D1,plot=FALSE)
  p3 = plot_dynSS(dat=out_f_noninf$stat) + annotate(geom="text",x=115, y=80, label="Beta(42,10)")
  p4 = plot_dynRand(dat=out_f_noninf$stat) + annotate(geom="text",x=115, y=0.6, label="Beta(42,10)")
  #grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,top = textGrob("Non-inferiority trial (arm dropping for futility)",gp=gpar(fontsize=12,font=3)))
  grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
}
##
##
################################################################################ 
################################################################################ 
