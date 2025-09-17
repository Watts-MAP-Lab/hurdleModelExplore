data {
int<lower=0> n;
int<lower=0> p;
int<lower=0,upper=1> Y[n,p];
}
parameters {
vector[n] theta;
vector<lower=0> [p] alpha;
vector[p] delta;
real mu_delta;
real<lower=0> sigma_alpha;
real<lower=0> sigma_delta;
}
transformed parameters{
vector<lower=0,upper=1>[p] prob[n];
for(i in 1:n){
 for (j in 1:p){
  prob[i,j] = inv_logit(alpha[j]*(theta[i] - delta[j]));
 }
}
}
model {
theta ~ normal(0,1);
delta ~ normal(mu_delta,sigma_delta);
mu_delta ~ normal(0,5);
sigma_delta ~ cauchy(0,5);
alpha ~ lognormal(0,sigma_alpha);
sigma_alpha ~ cauchy(0,5);
for(i in 1:n){
 for (j in 1:p){
  Y[i,j] ~ bernoulli(prob[i,j]);
 }
}
}
generated quantities {
vector[p] loglik_y[n];
vector[p] Y_rep[n];
for (i in 1: n){
 for (j in 1: p){
   loglik_y[i,j] = bernoulli_lpmf(Y[i,j] | prob[i,j]);
   Y_rep[i,j] = bernoulli_rng(prob[i,j]); 
 }
}
}
