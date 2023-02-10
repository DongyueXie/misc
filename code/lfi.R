
## generate data
gen_data = function(n,p,pi0,var0,var1,sigma2,XI = FALSE){
  beta0 = rnorm(p,0,sqrt(var0))
  beta1 = rnorm(p,0,sqrt(var1))
  gamma = rbinom(p,1,1-pi0)
  beta = beta1 * gamma + (1-gamma) * beta0
  if(XI){
    y = beta + rnorm(n,0,sqrt(sigma2))
  }else{
    X = matrix(rnorm(n*p),nrow=n,ncol=p)
    y = X%*%beta + rnorm(n,0,sqrt(sigma2))
  }
  return(list(beta=beta,y=y,X=X))
}

ssvs = function(X,y,pi0,var0,var1,sigma2,fix_sigma2=T,ig_a=0.01,ig_b=0.01,n_burnin=1000,n_post=5000,printevery = 10){
  n = nrow(X)
  p = ncol(X)
  XtX = crossprod(X,X)
  Xty = crossprod(X,y)
  beta_draws = matrix(nrow=n_burnin + n_post,ncol=p)
  gamma_draws = matrix(nrow=n_burnin + n_post,ncol=p)
  sigma2_draws = c()
  gamma = rep(0,p)
  if(!fix_sigma2){
    sigma2 = var(y)
  }
  for(i in 1:(n_burnin + n_post)){
    if(i%%printevery==0){
      print(paste('drawing sample',i))
    }
    d_inv = 1/(gamma*var1 + (1-gamma)*var0)
    A = solve(XtX/sigma2 + diag(d_inv))
    beta = c(mvnfast::rmvn(n=1,mu=A%*%Xty/sigma2,sigma=A))
    if(!fix_sigma2){
      sigma2 = invgamma::rinvgamma(n/2+ig_a,scale=sum((y-X%*%beta)^2)/2+ig_b)
    }
    for(j in 1:p){
      d1 = (gamma*var1 + (1-gamma)*var0)
      d0 = d1
      d1[j] = var1
      d0[j] = var0
      p1 = log(1-pi0) + sum(dnorm(beta,rep(0,p),sqrt(d1),log = T))
      p0 = log(pi0) + sum(dnorm(beta,rep(0,p),sqrt(d0),log = T))
      p_vec = c(p0,p1)
      p_vec = p_vec-max(p_vec)
      p_vec = exp(p_vec)/sum(exp(p_vec))
      gamma[j] = rbinom(1,1,p_vec[2])
    }
    beta_draws[i,] = beta
    gamma_draws[i,] = gamma
    sigma2_draws[i] = sigma2
  }
  return(list(beta_draws=beta_draws,gamma_draws=gamma_draws,sigma2_draws=sigma2_draws))
}

library(BBSSL)
source('code/sslasso_code/lasso_inference.r')

simu_study = function(n_simu = 100,
                      n = 50,
                      p = 50,
                      pi0 = 0.95,
                      var0 = 0.01,
                      var1 = 25,
                      sigma2 = 1){

  res = list()
  covered = matrix(nrow=n_simu,ncol=3)
  for(i in 1:n_simu){
    print(paste('Running simulation:',i))
    datax = gen_data(n,p,pi0,var0,var1,sigma2)
    fit_ssvs = ssvs(datax$X,datax$y,pi0,var0,var1,sigma2,n_burnin = 500,n_post = 2000,printevery = 500)


    lambda0 = 7;
    lambda1 = 0.15;
    lambda = c(lambda0, lambda1)
    a = 1;
    b = ncol(datax$X)

    fit_bb = BB_SSL(datax$y, datax$X, method = 3, lambda=lambda, NSample=2000, a=a, b=b, maxiter=100, eps = 1e-3, burn.in = T,
                  length.out = 50, discard = FALSE, alpha = 3, sigma = sqrt(sigma2), initial.beta = rep(0,p),
                  penalty = "adaptive", theta=0.5)

    fit_debiased = SSLasso(datax$X, datax$y, alpha = 0.05, lambda = NULL, mu = NULL, intercept = F,
                   resol = 1.3, maxiter = 50, threshold = 1e-2, verbose = TRUE)
    res[[i]] = list(data=datax,fit_ssvs=fit_ssvs,fit_bb=fit_bb,fit_debiased=fit_debiased)
    covered[i,1] = mean((datax$beta<c(apply(fit_ssvs$beta_draws,2,function(x){quantile(x,0.975)})))*
                          (datax$beta>c(apply(fit_ssvs$beta_draws,2,function(x){quantile(x,0.025)}))))
    covered[i,2] = mean((datax$beta<c(apply(fit_bb$beta,2,function(x){quantile(x,0.975)})))*
                          (datax$beta>c(apply(fit_bb$beta,2,function(x){quantile(x,0.025)}))))
    covered[i,3] = mean((datax$beta<fit_debiased$up.lim)*((datax$beta>fit_debiased$low.lim)))
  }
  colnames(covered) = c('ssvs','bb_ssl','debiased')
  return(list(res=res,covered=covered))
}

n = 100
p = 100
pi0 = 0.95
var0 = 0.01
var1 = 25
sigma2 = 1
output = simu_study(n=n,p=p,pi0=pi0,var0=var0,var1=var1,sigma2=sigma2)
saveRDS(output,file='output/lfi/n100.rds')
