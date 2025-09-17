data{
  int<lower=2, upper=8> K; //number of categories
  int <lower=0> n_student;
  int <lower=1> n_item;
  int<lower=1,upper=K> Y[n_student,n_item];
}
parameters {
vector[n_student] theta;
real<lower=0> alpha [n_item];
ordered[K-1] kappa[n_item]; //category difficulty
real mu_kappa; //mean of the prior distribution of category difficulty
real<lower=0> sigma_kappa; //sd of the prior distribution of category difficulty
}
model{
alpha ~ cauchy(0,5);
theta ~ normal(0,1);
for (i in 1: n_item){
  for (k in 1:(K-1)){
    kappa[i,k] ~ normal(mu_kappa,sigma_kappa);
  }
}
mu_kappa ~ normal(0,5);
sigma_kappa ~ cauchy(0,5);
for (i in 1:n_student){
  for (j in 1:n_item){
    Y[i,j] ~ ordered_logistic(theta[i]*alpha[j],kappa[j]);
  }
}
}
generated quantities {
vector[n_item] log_lik[n_student];
for (i in 1: n_student){
for (j in 1: n_item){
log_lik[i, j] = ordered_logistic_log(Y[i, j],theta[i]*alpha[j],kappa[j]);
}}
}