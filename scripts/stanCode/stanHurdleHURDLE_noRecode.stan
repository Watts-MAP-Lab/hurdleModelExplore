data {
  int<lower=0> N; // number of subjects
  int<lower=0> p; // number of items -- equivalent for 2PL and GRM items
  int y[N,p]; // response vector 
  int k; // number of categories including 0
  vector[2] Zero; // theta means
}

parameters {
  vector[2] theta[N]; // theta matrix
  cholesky_factor_corr[2] Lcorr; // cholesky factor of the corr matrix
  //vector<lower=0>[2] sigma; // factor cor value
  vector<lower=0>[p] alpha_pl; // discrim 2pl
  vector[p] beta_pl; // difficulty 2pl
  real mu_beta_pl; // mean difficulty
  real<lower=0> sigma_alpha_pl; // std discrim 2pl
  real<lower=0> sigma_beta_pl; // std difficulty 2pl
  real<lower=0> sigma_alpha_grm; // std discrim grm
  // now do the grm parameters down here
  vector<lower=0>[p] alpha_grm; // discrim grm items
  ordered[k - 2] kappa[p];
  real mu_kappa; //mean of the prior distribution of category difficulty
  real<lower=0> sigma_kappa; //sd of the prior distribution of category difficulty
}

## Now do the transformed paramters here
## if 0 --> 2PL probability
## if > 0 --> 2PL * GRM probability
model {
  // priors
  theta ~ multi_normal_cholesky(Zero, Lcorr);
  Lcorr ~ lkj_corr_cholesky(1);
  
  alpha_pl ~ normal(1,sigma_alpha_pl);
  sigma_alpha_pl ~ cauchy(0,5);
  alpha_grm ~ normal(1, sigma_alpha_grm);
  sigma_alpha_grm ~ cauchy(0,5);
  
  // difficulty priors
  mu_beta_pl ~ normal(1,.4);
  sigma_beta_pl ~ cauchy(0,3);
  beta_pl ~ normal(mu_beta_pl,sigma_beta_pl);
  for (i in 1:p){
    for (z in 1:(k-2)){
      kappa[i,z] ~ normal(mu_kappa,sigma_kappa);
    }
  }
  mu_kappa ~ normal(0,1); // grm difficulty hyperprior mean
  sigma_kappa ~ cauchy(0,3); // grm difficulty hyperprior std
  
  // model
  for(i in 1:N){
    for(j in 1:p){
        real eta_z = alpha_pl[j] * (theta[i,1] - beta_pl[j]);
        if (y[i, j] == 0) {
          target += bernoulli_logit_lpmf(0 | eta_z);
        } else {
          int r = y[i, j];  // 1..M
          target += bernoulli_logit_lpmf(1 | eta_z);
          target += ordered_logistic_lpmf(r | alpha_grm[j] * theta[i, 2], alpha_grm[j] * kappa[j]);
      }
    }
  }
}

