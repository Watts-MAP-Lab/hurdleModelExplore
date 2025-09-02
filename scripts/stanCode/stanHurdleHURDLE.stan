data {
  int<lower=0> N; // number of subjects
  int<lower=0> p; // number of items -- equivalent for 2PL and GRM items
  int Y[N,p]; // response vector -- binary
  // Now do the GRM data down here -- responses should be vectorized to accomodate missing data
  int k; // number of grm categories excluding 0
  int<lower=1> NG; // number of grm responses
  int<lower=1,upper=NG> jj[NG];  // respondent for observation n
  int<lower=1,upper=NG> kk[NG];  // question for observation n
  int<lower=0,upper=k> y[NG];   // correctness for observation n
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[N] theta1; // 2PL theta
  vector[N] theta2; // GRM theta
  vector<lower=0>[p] alpha_pl; // discrim 2pl
  vector[p] beta_pl; // difficulty 2pl
  real mu_beta_pl; // mean difficulty
  real<lower=0> sigma_alpha_pl; // std discrim 2pl
  real<lower=0> sigma_beta_pl; // std difficulty 2pl
  // now do the grm parameters down here
  vector<lower=0>[p] alpha_grm; // discrim grm items
  ordered[k-1] kappa[p]; //category difficulty
  real mu_kappa; //mean of the prior distribution of category difficulty
  real<lower=0> sigma_kappa; //sd of the prior distribution of category difficulty
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // assign priors
  // Theta priors
  theta1 ~ normal(0,1);
  theta2 ~ normal(0,1);
  // discrim priors
  alpha_pl ~ normal(0,sigma_alpha_pl);
  sigma_alpha_pl ~ cauchy(0,5);
  alpha_grm ~ normal(1,.4);

  // difficulty priors
  mu_beta_pl ~ normal(0,5);
  sigma_beta_pl ~ cauchy(0,5);
  beta_pl ~ normal(mu_beta_pl,sigma_beta_pl);
  for (i in 1:p){
    for (z in 1:(k-1)){
      kappa[i,z] ~ normal(mu_kappa,sigma_kappa);
    }
  }
  mu_kappa ~ normal(0,5); // grm difficulty hyperprior mean
  sigma_kappa ~ cauchy(0,5); // grm difficulty hyperprior std
  // now estimate the 2PL lik
  for(i in 1:N){
    for (j in 1:p){
       	real pq; //create a local variable within the loop to make Stan code more readable
       	pq= inv_logit(alpha_pl[j]*(theta1[i] - beta_pl[j]));
       	Y[i,j] ~ bernoulli(pq);
      }
  }
  // now esitmate GRM portion here
  for(i in 1:NG) {
    y[i] ~ ordered_logistic(theta2[jj[i]]*alpha_grm[kk[i]],kappa[kk[i]]);
  }
}