

rm(list = ls())

## prior effective dynamic sample size - code 

library(dplyr); library(data.table); library(spTimer); library(LaplacesDemon); 
library(extraDistr); library(pbapply); library(MASS);
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)  

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
gen_data_norm <- function(N, bta, sig, rand_pr=c(0.25,0.25,0.25,0.25)){
  ##
  input <- gen_data_arm(n=N, rand_pr = rand_pr, n_arm = NULL)
  ##
  xx <- model.matrix(~factor(input$id_arm, levels=c(1:length(rand_pr))))
  y <- xx%*%c(bta) + rnorm(N,0,sig)
  out <- data.frame(input, y = y, x=as.factor(input$id_arm))
  out
}
##
model <- stan_model(file="lm_beta.stan")
trial_LaplaceApprox <- function(data, 
                                model,
                                mu_beta = c(0,0,0,0), # with intercept -
                                sig_beta = c(5,5,5,5), # with intercept - 
                                a=0.5, b=2){
  ##
  ##
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)  
  X <- model.matrix( ~ as.factor(x), data = data)
  y <- data$y
  K <- ncol(X) # with interpect
  N <- nrow(X)
  if(!length(mu_beta)%in%K){
    stop("check mu_beta")
  }
  if(!length(sig_beta)%in%K){
    stop("check sig_beta")
  }
  ##
  lmdat <- list(
    N = N,
    K = K, # with interpect
    x = X,
    y = y,
    mu_beta = mu_beta, sig_beta=sig_beta,
    a = a, b = b
  )
  ##
  lmLA <- optimizing(model, data=lmdat, hessian = TRUE)
  # Covariance is negative inverse Hessian
  opt_cov <- solve(-lmLA$hessian)
  ##
  post.beta.mn <- lmLA$par[1:K]
  post.beta.sd <- sqrt(diag(opt_cov[1:K,1:K]))
  post.sigma2.mn <- lmLA$par[K+1]
  post.sigma2.sd <- sqrt(opt_cov[K+1,K+1])
  ##
  rm(lmLA)
  ##
  sample_size <- unique(sapply(data,length))
  ##
  return(list(post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              post.sigma2.mn=post.sigma2.mn, post.sigma2.sd=post.sigma2.sd,
              sample_size=sample_size))
}
##
##
EDSS_ <- function(dat=NULL, nRep=5000, 
                  prior_interest=list(c(1,0.5),c(1,0.5)), # list function for all arms 
                  prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                  n_ref=50, n_grid=seq(2,50,1), bta=c(0.5,1), sig=0.5,
                  rand_pr=c(0.5,0.5), model=model, onlyValStat=TRUE){
  ##
  ## prior_baseline = reference prior (may not be vauge)
  ## n_ref = sample size for the baseline/reference prior
  ## n_grid = grid sample size to check (bit different compared to the 2019 ECSS paper, where they used max(n_grid) = n_ref)
  ##
  if(is.null(dat)){
    dat <- gen_data_norm(N=max(c(n_ref,n_grid)), bta=bta, sig=sig, rand_pr=rand_pr)
  }
  else{
    dat <- dat
  }
  ##
  prior_baseline <- as.matrix(data.frame(prior_baseline))
  dimnames(prior_baseline) <- NULL
  ##
  res_ref <- trial_LaplaceApprox(data=dat[1:n_ref,], #dat_ref, 
                                 model=model, 
                                 mu_beta = c(prior_baseline[1,]), # with intercept -
                                 sig_beta = sqrt(c(prior_baseline[2,])), # with intercept - 
                                 a=3, b=1) # for sig prior
  ##
  s <- matrix(0,length(bta),length(bta))
  diag(s) <- res_ref$post.beta.sd
  bint_ref <- mvrnorm(nRep,res_ref$post.beta.mn,s)
  b_mse <- mean(rowSums((t(bint_ref)-bta)^2)/nRep)
  b_kld <- KLD(px=c(bint_ref),py=rep(bta,each=nRep))$mean.sum.KLD
  ##
  ## store MSE_b, MSE_p and KL-divergence (mean),
  store_valStat <- matrix(NA,length(n_grid),6) 
  store_res <- list()
  prior_interest <- as.matrix(data.frame(prior_interest))
  dimnames(prior_interest) <- NULL
  for(j in 1:length(n_grid)){
    #if(j==1){
    #  dyn_prior_interest <- prior_interest
    #}
    res <- trial_LaplaceApprox(data=dat[1:n_grid[j],], #dat[[j]], 
                               model=model, 
                               mu_beta = c(prior_interest[1,]), # with intercept -
                               #mu_beta = c(dyn_prior_interest[1,]), # with intercept -
                               sig_beta = sqrt(c(prior_interest[2,])), # with intercept - 
                               a=3, b=1) # for sig prior
    dyn_prior_interest <- rbind(c(res$post.beta.mn),c(res$post.beta.sd))
    s <- matrix(0,length(bta),length(bta))
    diag(s) <- res$post.beta.sd
    bint <- mvrnorm(nRep,res$post.beta.mn,s)
    p_mse <- mean(rowSums((t(bint)-bta)^2)/nRep)
    p_kld <- KLD(px=c(bint),py=rep(bta,each=nRep))$mean.sum.KLD
    bp_kld <- KLD(px=c(bint),py=c(bint_ref))$mean.sum.KLD
    store_valStat[j,] <- c(n_grid[j],p_mse,b_mse,p_kld,b_kld,bp_kld)
    store_res[[j]] <- res
    rm(res)
  }
  store_valStat <- data.frame(store_valStat)
  names(store_valStat) <- c("Sample_size","MSE","MSE_b","KL_mean","KL_mean_b","KL_mean_bp")
  if(isTRUE(onlyValStat)){
    return(list(valStat=store_valStat,dyn_prior_interest=dyn_prior_interest))
  }
  else{
    return(list(valStat=store_valStat,result=store_res,res_ref=res_ref,n_grid=n_grid,
                dyn_prior_interest=dyn_prior_interest))
  }
}
##
##
EDSS_simulation <- function(nSim=10, dat=NULL, nRep=5000, 
                            prior_interest=list(c(1,0.5),c(1,0.5)), # list function for all arms 
                            prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                            n_ref=50, n_grid=seq(30,50,1), bta=c(0.5,1), sig=0.5,
                            rand_pr=c(0.5,0.5), model=model){
  ##
  library(pbapply)
  out <- pblapply(1:nSim, 
                  function(x) EDSS_(dat=dat , nRep=nRep, 
                                    prior_interest=prior_interest, # only for the treatment arm
                                    prior_baseline=prior_baseline, # only for the treatment arm 
                                    n_ref=n_ref, n_grid=n_grid, bta=bta, sig=sig,
                                    rand_pr=rand_pr, model=model))
  ##
  valStat <- sapply(out, function(x) as.matrix(x$valStat))
  mn <- rowMeans(valStat)
  se <- apply(valStat,1,sd)
  valStat <- array(c(mn,se),dim=c(dim(out[[1]]$valStat),2))
  dimnames(valStat)[[2]] <- names(out[[1]]$valStat)
  dimnames(valStat)[[3]] <- c("mean","se")
  ##
  dyn_prior_interest <- sapply(out, function(x) as.matrix(x$dyn_prior_interest))
  dyn_prior_interest <- matrix(rowMeans(dyn_prior_interest),2,length(bta))
  dimnames(dyn_prior_interest)[[1]] <- c("mean","sd")
  ##
  return(list(valStat=valStat,n_ref=n_ref,n_grid=n_grid,
              dyn_prior_interest=dyn_prior_interest))
  ##
}
##
## null model 
#out <- EDSS_simulation(nSim=100,model=model,
#                       prior_interest=list(c(0,1000),c(0,1)),
#                       prior_baseline=list(c(0,1000),c(0,1000)),
#                       bta=c(0,0),n_ref=100, n_grid=seq(50,100,1))
#mn <- data.frame(out$valStat[,,1])
##
#ggplot(mn,aes(x=Sample_size,y=MSE)) +
#  geom_point() + geom_smooth(span=0.3) +
#  geom_line(aes(x=Sample_size,y=MSE_b))
#ggplot(mn,aes(x=Sample_size,y=KL_mean)) +
#  geom_point() + geom_smooth(span=0.3) +
#  geom_line(aes(x=Sample_size,y=KL_mean_b))
##
#mn$mse_diff <- abs(mn$MSE-mn$MSE_b)
#ggplot(mn,aes(x=Sample_size,y=MSE)) +
#  geom_point() + geom_smooth(span=0.3) +
#  geom_line(aes(x=Sample_size,y=0))
#mn[which(mn$mse_diff==min(mn$mse_diff)),]
#mn$kl_diff <- abs(mn$KL_mean-mn$KL_mean_b)
#ggplot(mn,aes(x=Sample_size,y=kl_diff)) +
#  geom_point() + geom_smooth(span=0.3) +
#  geom_line(aes(x=Sample_size,y=0))
#mn[which(mn$kl_diff==min(mn$kl_diff)),]

##
## null model - interim 1
out1 <- EDSS_simulation(nSim=500,model=model,
                       prior_interest=list(c(0,1000),c(0,10)),
                       prior_baseline=list(c(0,1000),c(0,1000)),
                       bta=c(0,0),n_ref=100, n_grid=seq(50,100,1))
mn1 <- data.frame(out1$valStat[,,1])
mn1$kl_diff <- abs(mn1$KL_mean-mn1$KL_mean_b)
ggplot(mn1,aes(x=Sample_size,y=kl_diff)) +
  geom_smooth(span=0.3) + geom_line(aes(x=Sample_size,y=0)) +
  labs(y="Kullback–Leibler divergence (KLD) differences")
optmSS <- mn1[which(mn1$kl_diff==min(mn1$kl_diff)),]
dyn_prior_interest <- out1$dyn_prior_interest
## null model - interim 2
out2 <- EDSS_simulation(nSim=500,model=model,
                        prior_interest=list(c(0,1000),c(0,10)), # ECSS
                        #prior_interest=list(c(dyn_prior_interest[,1]),c(dyn_prior_interest[,2])),  # EDSS
                        prior_baseline=list(c(0,1000),c(0,1000)),
                        bta=c(0,0),n_ref=optmSS[1,1]+100, n_grid=seq(optmSS[1,1],optmSS[1,1]+100,1))
mn2 <- data.frame(out2$valStat[,,1])
mn2$kl_diff <- abs(mn2$KL_mean-mn2$KL_mean_b)
ggplot(mn2,aes(x=Sample_size,y=kl_diff)) +
  geom_smooth(span=0.3) + geom_line(aes(x=Sample_size,y=0)) +
  labs(y="Kullback–Leibler divergence (KLD) differences")
optmSS <- rbind(optmSS,mn2[which(mn2$kl_diff==min(mn2$kl_diff)),])
dyn_prior_interest <- out2$dyn_prior_interest
## null model - interim 3
out3 <- EDSS_simulation(nSim=500,model=model,
                        #prior_interest=list(c(0,1000),c(0,10)), # ECSS
                        prior_interest=list(c(dyn_prior_interest[,1]),c(dyn_prior_interest[,2])),  # EDSS
                        prior_baseline=list(c(0,1000),c(0,1000)),
                        bta=c(0,0),n_ref=optmSS[2,1]+100, n_grid=seq(optmSS[2,1],optmSS[2,1]+100,1))
mn3 <- data.frame(out3$valStat[,,1])
mn3$kl_diff <- abs(mn3$KL_mean-mn3$KL_mean_b)
ggplot(mn3,aes(x=Sample_size,y=kl_diff)) +
  geom_smooth(span=0.3) + geom_line(aes(x=Sample_size,y=0)) +
  labs(y="Kullback–Leibler divergence (KLD) differences")
optmSS <- rbind(optmSS,mn3[which(mn3$kl_diff==min(mn3$kl_diff)),])
##
optmSS_EDSS <- optmSS
rm(optmSS)
##
#save.image(file="prior_edss.RData")
##

ggplot() +
  geom_smooth(data=mn1,aes(x=Sample_size,y=kl_diff,color="Interim 1"),span=0.3) + 
  geom_smooth(data=mn2,aes(x=Sample_size,y=kl_diff,color="Interim 2"),span=0.3) + 
  geom_smooth(data=mn3,aes(x=Sample_size,y=kl_diff,color="Interim 3"),span=0.3) + 
  scale_color_manual(name=' ',
                     breaks=c('Interim 1', 'Interim 2', 'Interim 3'),
                     values=c('Interim 1'='red', 'Interim 2'='blue', 'Interim 3'='green')) +
  geom_line(data=mn1, aes(x=Sample_size,y=0)) +
  geom_line(data=mn2, aes(x=Sample_size,y=0)) +
  geom_line(data=mn3, aes(x=Sample_size,y=0)) +
  labs(y="Kullback–Leibler divergence (KLD) differences",x="Sample size")


ggplot() +
  geom_smooth(data=mn1,aes(x=Sample_size,y=kl_diff),span=0.3,col=2) + #geom_point(data=mn1,aes(x=Sample_size,y=kl_diff),col=2,pch=2) +
  geom_smooth(data=mn2,aes(x=Sample_size,y=kl_diff),span=0.3,col=3) + #geom_point(data=mn2,aes(x=Sample_size,y=kl_diff),col=3,pch=3) +
  geom_smooth(data=mn3,aes(x=Sample_size,y=kl_diff),span=0.3,col=4) + #geom_point(data=mn3,aes(x=Sample_size,y=kl_diff),col=4,pch=4) +
  geom_line(data=mn1, aes(x=Sample_size,y=0)) +
  geom_line(data=mn2, aes(x=Sample_size,y=0)) +
  geom_line(data=mn3, aes(x=Sample_size,y=0)) +
  labs(y="Kullback–Leibler divergence (KLD) differences")

##  
load(file="prior_edss.RData")
##

##
## code for effective dynamic sample size (EDSS) for adaptive desing

trail_interim_EDSS <- function(model,n_sample=c(100,50,50),
                               rand_assignment=c(1,1),
                               nRep=5000, 
                               prior_interest=list(c(0,1000),c(0,1)), # list function for all arms 
                               prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                               true_beta = c(0.5,1), # first input is intercept
                               true_sig = c(1),
                               EDSS = "mse", # kld or mse
                               capType="noCap",
                               D0=0.10, D1=0.90){
  ##
  K <- length(true_beta)
  n_j <- length(n_sample)
  dat <- list()
  res <- list()
  optm <- list()
  optmSS <- c()
  n_rand_pr <- matrix(NA,K,n_j)
  pr_decision <- matrix(NA,K-1,n_j)
  n_interimSS <- matrix(NA,K,n_j)
  dyn_prior_interest <- list()
  ##
  for(j in 1:n_j){
    ##
    if(j==1){
      ##
      rand_pr <- rand_assignment/sum(rand_assignment)
      dat[[j]] = gen_data_norm(N=n_sample[j], bta=true_beta, sig=true_sig, rand_pr=rand_pr)
      n_ref <- nrow(dat[[j]]); n_grid <- seq(round(nrow(dat[[j]])/2),nrow(dat[[j]]),1)
      optm[[j]] = EDSS_(dat=dat[[j]], nRep=nRep, 
                         prior_interest=prior_interest, # only for the treatment arm
                         prior_baseline=prior_baseline, # only for the treatment arm 
                         n_ref=n_ref, n_grid=n_grid, bta=true_beta, sig=true_sig,
                         rand_pr=rand_pr, model=model)
      mn <- data.frame(optm[[j]]$valStat)
      mn$mse_diff <- abs(mn$MSE-mn$MSE_b)
      mn$kl_diff <- abs(mn$KL_mean-mn$KL_mean_b)
      if(EDSS%in%"mse"){
        opt_m <- mn[which(mn$mse_diff==min(mn$mse_diff)),1]
      }
      else if(EDSS%in%"kld"){
        opt_m <- mn[which(mn$kl_diff==min(mn$kl_diff)),1]
      }
      else{
        stop("define EDSS correctly, can take arguments mse or kld ")
      }
      optmSS[j] <- opt_m
      #dat[[j]] <- dat[[j]][1:opt_m,]
      #dat[[j]]$id <- 1:nrow(dat[[j]])
      dyn_prior_interest[[j]] = as.matrix(data.frame(prior_interest))
      dimnames(dyn_prior_interest) <- NULL
      res[[j]] = trial_LaplaceApprox(data=dat[[j]][1:opt_m,], model=model, #data=dat[[j]], model=model, 
                                     mu_beta = c(dyn_prior_interest[[j]][1,]), # with intercept -
                                     sig_beta = c(dyn_prior_interest[[j]][2,]), # with intercept - 
                                     a=3, b=1) # for sig prior
      dyn_prior_interest[[j]] <- rbind(c(res[[j]]$post.beta.mn),c(res[[j]]$post.beta.sd))
      #pr_decision[,j] <- pnorm(0, res$post.beta.mn[2], res$post.beta.sd[2]) # Pr(b < 0)
      pr_decision[,j] <- 1-pnorm(0, res[[j]]$post.beta.mn[2], res[[j]]$post.beta.sd[2]) # Pr(b > 0)
      n_rand_pr[,j] <- rand_pr
      n_interimSS[,j] <- table(dat[[j]][1:opt_m,]$x)
      ##
    }
    else{
      ## For no sample size cap: if one arm is dropped then other arms 
      ## may receive more participants up to the nMax[sum(n_sample)] participants 
      #if(capType=="noCap"){
      #  ## need ot work on the code below 
      #  drop_arm <- c(which(pr.decision[,j-1] < D0) + 1,which(pr.decision[,j-1] > D1) + 1)
      #  rand_pr[drop_arm] <- 0
      #  rand_pr <- rand_pr/sum(rand_pr)
      #}
      #else{
      rand_pr <- rand_assignment/sum(rand_assignment)
      #}
      ##
      gen_sample <- optmSS[j-1]+n_sample[j]-nrow(dat[[j-1]])
      dat[[j]] = data.frame(mapply(c,dat[[j-1]],gen_data_norm(N=gen_sample, bta=true_beta, sig=true_sig, rand_pr=rand_pr),SIMPLIFY=FALSE))
      dat[[j]]$id <- 1:nrow(dat[[j]])
      n_ref <- n_sample[j]; n_grid <- seq(round(n_sample[j]/2),n_sample[j],1)
      #n_ref <- nrow(dat[[j]]); n_grid <- seq(optmSS[j-1]+1,nrow(dat[[j]]),1)
      #n_ref <- nrow(dat[[j]]); n_grid <- seq(max(c(round(nrow(dat[[j]])-cumsum(n_sample)[j]/2),nrow(dat[[j-1]])+1)),nrow(dat[[j]]),1)
      optm[[j]] = EDSS_(dat=dat[[j]][(optmSS[j-1]+1):nrow(dat[[j]]),], nRep=nRep, 
                         prior_interest=prior_interest, # only for the treatment arm
                         prior_baseline=prior_baseline, # only for the treatment arm 
                         n_ref=n_ref, n_grid=n_grid, bta=true_beta, sig=true_sig,
                         rand_pr=rand_pr, model=model)
      mn <- data.frame(optm[[j]]$valStat)
      mn$Sample_size <- mn$Sample_size + optmSS[j-1]
      mn$mse_diff <- abs(mn$MSE-mn$MSE_b)
      mn$kl_diff <- abs(mn$KL_mean-mn$KL_mean_b)
      if(EDSS%in%"mse"){
        opt_m <- mn[which(mn$mse_diff==min(mn$mse_diff)),1]
      }
      else if(EDSS%in%"kld"){
        opt_m <- mn[which(mn$kl_diff==min(mn$kl_diff)),1]
      }
      else{
        stop("define EDSS correctly, can take arguments mse or kld ")
      }
      optmSS[j] <- opt_m
      #dat[[j]] <- dat[[j]][1:opt_m,]
      #dat[[j]]$id <- 1:nrow(dat[[j]])
      res[[j]] = trial_LaplaceApprox(data=dat[[j]][1:opt_m,], model=model, 
                                     mu_beta = c(dyn_prior_interest[[j-1]][1,]), # with intercept -
                                     sig_beta = c(dyn_prior_interest[[j-1]][2,]), # with intercept - 
                                     a=3, b=1) # for sig prior
      dyn_prior_interest[[j]] <- rbind(c(res[[j]]$post.beta.mn),c(res[[j]]$post.beta.sd))
      #pr_decision[,j] <- pnorm(0, res$post.beta.mn[2], res$post.beta.sd[2]) # Pr(b < 0)
      pr_decision[,j] <- 1-pnorm(0, res[[j]]$post.beta.mn[2], res[[j]]$post.beta.sd[2]) # Pr(b > 0)
      n_rand_pr[,j] <- rand_pr
      n_interimSS[,j] <- table(dat[[j]][1:opt_m,]$x)      
      ##
    }
  }
  ##
  dimnames(n_rand_pr)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_rand_pr)[[2]] = c(paste0("interim",1:n_j))
  dimnames(pr_decision)[[1]] = c(paste0("treat",1:(length(true_beta)-1)))
  ##
  dimnames(n_interimSS)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_interimSS)[[2]] = c(paste0("interim",1:n_j))
  ##
  post.beta.mn <- sapply(res, function(x) x$post.beta.mn)
  post.beta.sd <- sapply(res, function(x) x$post.beta.sd)
  post.sigma2.mn <- sapply(res, function(x) x$post.sigma2.mn)
  post.sigma2.sd <- sapply(res, function(x) x$post.sigma2.sd)
  dimnames(post.beta.mn)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(post.beta.sd)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  ##
  dyn_prior_interest <- array(sapply(dyn_prior_interest, function(x) x),dim=c(2,K,n_j))
  ##
  return(list(pr_decision=pr_decision,
              post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              post.sigma2.mn=post.sigma2.mn, post.sigma2.sd=post.sigma2.sd,
              sample_size=n_interimSS, rand_pr=n_rand_pr, 
              optmSS=optmSS, dyn_prior_interest=dyn_prior_interest)) #,optm=optm))
  ##
}
##
#out <- trail_interim_EDSS(model=model,
#                          n_sample=c(100,50,50),
#                          rand_assignment=c(1,1),   
#                          true_beta = c(0.5,1), # first input is intercept
#                          true_sig = c(1),
#                          prior_interest=list(c(0,1000),c(0,1000)), # list function for all arms 
#                          prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
#                          EDSS = "mse",
#                          capType="noCap",
#                          D0=0.10, D1=0.90)
##
#out$pr_decision
#out$sample_size
#out$post.beta.mn
#out$optmSS
##
## simulator 

trail_interim_EDSS_simulator <- function(nSim,model,n_sample=c(100,100,100),
                               rand_assignment=c(1,1), nRep=5000, 
                               prior_interest=list(c(0,1000),c(0,1)), # list function for all arms 
                               prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                               true_beta = c(0,0), # first input is intercept
                               true_sig = c(0.5),
                               EDSS = "kld", # kld or mse
                               capType="noCap",
                               D0=0.10, D1=0.90){
  ##
  library(pbapply)
  res <- pblapply(1:nSim, 
                  function(x) trail_interim_EDSS(model=model,
                                                 n_sample=n_sample,
                                                 rand_assignment=rand_assignment,   
                                                 true_beta = true_beta, # first input is intercept
                                                 true_sig = true_sig,
                                                 prior_interest=prior_interest, # list function for all arms 
                                                 prior_baseline=prior_baseline, # list function for all arms 
                                                 EDSS = EDSS,
                                                 capType=capType,
                                                 D0=D0, D1=D1))
  ##
  K <- length(true_beta)
  n_j <- length(n_sample)
  ##
  post.beta.mn <- sapply(res, function(x) x$post.beta.mn)
  post.beta.mn <- array(c(post.beta.mn),dim=c(K,n_j,nSim))
  dimnames(post.beta.mn)[[1]] <- c("control",paste0("treat",1:(K-1)))
  dimnames(post.beta.mn)[[2]] <- c(paste0("interim",1:n_j))
  post.beta.sd <- sapply(res, function(x) x$post.beta.sd)
  post.beta.sd <- array(c(post.beta.sd),dim=c(K,n_j,nSim))
  dimnames(post.beta.sd)[[1]] <- c("control",paste0("treat",1:(K-1)))
  dimnames(post.beta.sd)[[2]] <- c(paste0("interim",1:n_j))
  post.sigma2.mn <- sapply(res, function(x) x$post.sigma2.mn)
  dimnames(post.sigma2.mn)[[1]] <- c(paste0("interim",1:n_j))
  post.sigma2.sd <- sapply(res, function(x) x$post.sigma2.sd)
  dimnames(post.sigma2.sd)[[1]] <- c(paste0("interim",1:n_j))
  ##
  pr_decision <- sapply(res, function(x) x$pr_decision)
  pr_decision <- array(c(pr_decision),dim=c(K-1,n_j,nSim))
  dimnames(pr_decision)[[1]] <- c(paste0("treat",1:(K-1)))
  dimnames(pr_decision)[[2]] <- c(paste0("interim",1:n_j))
  ##
  sample_size <- sapply(res, function(x) x$sample_size)
  sample_size <- array(c(sample_size),dim=c(K,n_j,nSim))
  dimnames(sample_size)[[1]] <- c("control",paste0("treat",1:(K-1)))
  dimnames(sample_size)[[2]] <- c(paste0("interim",1:n_j))
  ##
  rand_pr <- sapply(res, function(x) x$rand_pr)
  rand_pr <- array(c(rand_pr),dim=c(K,n_j,nSim))
  dimnames(rand_pr)[[1]] <- c("control",paste0("treat",1:(K-1)))
  dimnames(rand_pr)[[2]] <- c(paste0("interim",1:n_j))
  ##
  optm_SS <- sapply(res, function(x) x$optmSS)
  dimnames(optm_SS)[[1]] <- c(paste0("interim",1:n_j))
  #apply(optm_SS,1,mean)
  ##
  dyn_prior_interest <- array(sapply(res, function(x) x$dyn_prior_interest),dim=c(2,K,n_j,nSim))
  dimnames(dyn_prior_interest)[[1]] = c("mean","sd")
  dimnames(dyn_prior_interest)[[2]] = c("control",paste0("treat",1:(K-1)))
  dimnames(dyn_prior_interest)[[3]] = c(paste0("interim",1:n_j))
  ##
  return(list(pr_decision=pr_decision,sample_size=sample_size,optm_SS=optm_SS,
              rand_pr=rand_pr, post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              post.sigma2.mn=post.sigma2.mn, post.sigma2.sd=post.sigma2.sd, 
              dyn_prior_interest=dyn_prior_interest))
  ##
}

##

fnc_replace_to_one <- function(x){
  ## x is a vector 
  if(length(which(x==1))>0){
    x[which(x==1)[1]:length(x)] <- 1
  }
  x
}

## null model

res <- trail_interim_EDSS_simulator(nSim=50, model=model, n_sample=c(100,100,100),
                                    prior_interest=list(c(0,1000),c(0,10)), # list function for all arms 
                                    prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                                    true_beta = c(0,0), # first input is intercept
                                    true_sig = c(0.5))
res2 <- trail_interim_EDSS_simulator(nSim=50, model=model, n_sample=c(100,100,100),
                                    prior_interest=list(c(0,1000),c(5,10)), # list function for all arms 
                                    prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                                    true_beta = c(0,0), # first input is intercept
                                    true_sig = c(0.5))

##
#save.image(file="prior_edss.RData")
##
##  
load(file="prior_edss.RData")
##
names(res)
apply(res$optm_SS,1,mean)
fnc_replace_to_one <- function(x){
  ## x is a vector 
  if(length(which(x==1))>0){
    x[which(x==1)[1]:length(x)] <- 1
  }
  x
}
D0=0.1; D1=0.9
tmp <- res$pr_decision[1,,]<D0
rowMeans(apply(tmp,2,fnc_replace_to_one))
tmp <- res$pr_decision[1,,]>D1
rowMeans(apply(tmp,2,fnc_replace_to_one))
rm(tmp)
##




##
## summary stat
##
dat <- data.table::CJ(
  interim = 1:length(n_sample),
  arm = 1:length(rand_assignment),
  id_scenario = NA_integer_,
  per_missing = NA_real_,
  per_missing_sd = NA_real_,
  true_beta = NA_real_,
  est_beta = NA_real_,
  est_beta_sd = NA_real_,
  true_sig = NA_real_,
  est_sig = NA_real_,
  est_sig_sd = NA_real_,
  est_sample =  NA_real_,
  est_sd =  NA_real_
)
#id_scenario <- id_scenario + 1
#dat$id_scenario <-  id_scenario
dat$id_scenario <-  1
dat$true_beta = rep(c(true_beta),length(n_sample))
dat$true_sig = true_sig
for(j in 1:length(n_sample)){
  #
  dat[dat$interim==j,"est_beta"] <- rowMeans(res$post.beta.mn[,j,])
  dat[dat$interim==j,"est_beta_sd"] <- rowMeans(res$post.beta.sd[,j,])
  #
  dat[dat$interim==j,"est_sig"] <-  mean(res$post.sigma.mn[,j])
  dat[dat$interim==j,"est_sig_sd"] <- mean(res$post.sigma.sd[,j])
  #
  dat[dat$interim==j,"est_sample"] <- rowMeans(res$sample_size[,j,])
  dat[dat$interim==j,"est_sd"] <- apply(res$sample_size[,j,],1,sd)
  #
}
## dim(res$pr.decision) # beta, interim, nSim
pr_decisionD0 <- NULL
pr_decisionD1 <- NULL
for(i in 1:(length(rand_assignment)-1)){
  tmp <- (res$pr.decision[i,,]<D0)
  pr_decisionD0 <- rbind(pr_decisionD0,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
  tmp <- (res$pr.decision[i,,]>D1)
  pr_decisionD1 <- rbind(pr_decisionD1,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
}
pr_decision <- cbind(rep(1:length(n_sample),(length(rand_assignment)-1)),c(pr_decisionD0),c(pr_decisionD1))
dimnames(pr_decision)[[2]] <- c("interim","D0_0.1","D1_0.975")
pr_decision <- data.frame(pr_decision)
pr_decision$arm <- rep(2:length(rand_assignment),each=length(n_sample))
##
dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
##
store_res <- rbind(store_res,dat)




########################################################################
### back up
########################################################################


trail_interim_EDSSsim <- function(model,n_sample=c(100,50,50),
                               rand_assignment=c(1,1),
                               nSim=1000, nRep=5000, 
                               prior_interest=list(c(1,0.5),c(1,0.5)), # list function for all arms 
                               prior_baseline=list(c(0,1000),c(0,1000)), # list function for all arms 
                               true_beta = c(0.5,1), # first input is intercept
                               true_sig = c(1),
                               EDSS = "mse", # kld or mse
                               capType="noCap",
                               D0=0.10, D1=0.90){
  ##
  K <- length(true_beta)
  n_j <- length(n_sample)
  dat <- list()
  res <- list()
  optm <- list()
  n_rand_pr <- matrix(NA,K,n_j)
  pr_decision <- matrix(NA,K-1,n_j)
  n_interimSS <- matrix(NA,K,n_j)
  ##
  for(j in 1:n_j){
    ##
    if(j==1){
      ##
      rand_pr <- rand_assignment/sum(rand_assignment)
      dat[[j]] = gen_data_norm(N=n_sample[j], bta=true_beta, sig=true_sig, rand_pr=rand_pr)
      n_ref <- nrow(dat[[j]]); n_grid <- seq(round(nrow(dat[[j]])/2),nrow(dat[[j]]),1)
      optm[[j]] <- EDSS_simulation(dat=dat[[j]], nSim=nSim, nRep=nRep, n_ref=n_ref, n_grid=n_grid, 
                                   bta=true_beta, sig=true_sig, model=model, rand_pr=rand_pr,
                                   prior_interest=prior_interest,  prior_baseline=prior_baseline)
      mn <- data.frame(optm[[j]]$valStat[,,1])
      mn$mse_diff <- abs(mn$MSE-mn$MSE_b)
      mn$kl_diff <- abs(mn$KL_mean-mn$KL_mean_b)
      if(EDSS%in%"mse"){
        opt_m <- mn[which(mn$mse_diff==min(mn$mse_diff)),1]
      }
      else if(EDSS%in%"kld"){
        opt_m <- mn[which(mn$kl_diff==min(mn$kl_diff)),1]
      }
      else{
        stop("define EDSS correctly, can take arguments mse or kld ")
      }
      dat[[j]] <- dat[[j]][1:opt_m,]
      dat[[j]]$id <- 1:nrow(dat[[j]])
      dyn_prior_interest <- as.matrix(data.frame(prior_interest))
      dimnames(dyn_prior_interest) <- NULL
      res[[j]] = trial_LaplaceApprox(data=dat[[j]], model=model, #data=dat[[j]], model=model, 
                                     mu_beta = c(dyn_prior_interest[1,]), # with intercept -
                                     sig_beta = sqrt(c(dyn_prior_interest[2,])), # with intercept - 
                                     a=3, b=1) # for sig prior
      dyn_prior_interest <- rbind(c(res[[j]]$post.beta.mn),c(res[[j]]$post.beta.sd))
      #pr_decision[,j] <- pnorm(0, res$post.beta.mn[2], res$post.beta.sd[2]) # Pr(b < 0)
      pr_decision[,j] <- 1-pnorm(0, res[[j]]$post.beta.mn[2], res[[j]]$post.beta.sd[2]) # Pr(b > 0)
      n_rand_pr[,j] <- rand_pr
      n_interimSS[,j] <- table(dat[[j]]$x)
      ##
    }
    else{
      ## For no sample size cap: if one arm is dropped then other arms 
      ## may receive more participants up to the nMax[sum(n_sample)] participants 
      #if(capType=="noCap"){
      #  ## need ot work on the code below 
      #  drop_arm <- c(which(pr.decision[,j-1] < D0) + 1,which(pr.decision[,j-1] > D1) + 1)
      #  rand_pr[drop_arm] <- 0
      #  rand_pr <- rand_pr/sum(rand_pr)
      #}
      #else{
      rand_pr <- rand_assignment/sum(rand_assignment)
      #}
      ##
      dat[[j]] = data.frame(mapply(c,dat[[j-1]],gen_data_norm(N=n_sample[j], bta=true_beta, sig=true_sig, rand_pr=rand_pr),SIMPLIFY=FALSE))
      dat[[j]]$id <- 1:nrow(dat[[j]])
      n_ref <- nrow(dat[[j]]); n_grid <- seq(max(c(round(nrow(dat[[j]])/2),nrow(dat[[j-1]]))),nrow(dat[[j]]),1)
      optm[[j]] <- EDSS_simulation(dat=dat[[j]], nSim=nSim, nRep=nRep, n_ref=n_ref, n_grid=n_grid, 
                                   bta=true_beta, sig=true_sig, model=model, rand_pr=rand_pr,
                                   prior_interest=prior_interest,  prior_baseline=prior_baseline)
      mn <- data.frame(optm[[j]]$valStat[,,1])
      mn$mse_diff <- abs(mn$MSE-mn$MSE_b)
      mn$kl_diff <- abs(mn$KL_mean-mn$KL_mean_b)
      if(EDSS%in%"mse"){
        opt_m <- mn[which(mn$mse_diff==min(mn$mse_diff)),1]
      }
      else if(EDSS%in%"kld"){
        opt_m <- mn[which(mn$kl_diff==min(mn$kl_diff)),1]
      }
      else{
        stop("define EDSS correctly, can take arguments mse or kld ")
      }
      dat[[j]] <- dat[[j]][1:opt_m,]
      res[[j]] = trial_LaplaceApprox(data=dat[[j]], model=model, 
                                     mu_beta = c(dyn_prior_interest[1,]), # with intercept -
                                     sig_beta = sqrt(c(dyn_prior_interest[2,])), # with intercept - 
                                     a=3, b=1) # for sig prior
      dyn_prior_interest <- rbind(c(res[[j]]$post.beta.mn),c(res[[j]]$post.beta.sd))
      #pr_decision[,j] <- pnorm(0, res$post.beta.mn[2], res$post.beta.sd[2]) # Pr(b < 0)
      pr_decision[,j] <- 1-pnorm(0, res[[j]]$post.beta.mn[2], res[[j]]$post.beta.sd[2]) # Pr(b > 0)
      n_rand_pr[,j] <- rand_pr
      n_interimSS[,j] <- table(dat[[j]]$x)      
      ##
    }
  }
  ##
  dimnames(n_rand_pr)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_rand_pr)[[2]] = c(paste0("interim",1:n_j))
  dimnames(pr_decision)[[1]] = c(paste0("treat",1:(length(true_beta)-1)))
  ##
  dimnames(n_interimSS)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_interimSS)[[2]] = c(paste0("interim",1:n_j))
  ##
  post.beta.mn <- sapply(res, function(x) x$post.beta.mn)
  post.beta.sd <- sapply(res, function(x) x$post.beta.sd)
  post.sigma.mn <- sapply(res, function(x) x$post.sigma.mn)
  post.sigma.sd <- sapply(res, function(x) x$post.sigma.sd)
  dimnames(post.beta.mn)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(post.beta.sd)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  ##
  return(list(pr_decision=pr_decision,
              post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              post.sigma.mn=post.sigma.mn, post.sigma.sd=post.sigma.sd,
              sample_size=n_interimSS, rand_pr=n_rand_pr, optm=optm))
  ##
}
##