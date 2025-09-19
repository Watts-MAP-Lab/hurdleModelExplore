data {
  int<lower=1> N;                 // persons
  int<lower=1> p;                 // items
  int<lower=2> k;                 // total categories including 0 (positive: 1..k-1)
  int<lower=0, upper=k-1> y[N, p];

  vector[2] Zero;                 // latent means (usually zeros)
}

parameters {
  // Latent traits: unit variances, free correlation
  vector[2] theta[N];
  cholesky_factor_corr[2] Lcorr;

  // 2PL item parameters
  vector[p] beta_pl;              // difficulty
  vector[p] log_alpha_pl;         // log discrimination (alpha_pl = exp(...))

  // GRM item parameters
  vector[p] log_alpha_grm;        // log discrimination (alpha_grm = exp(...))
  ordered[k - 2] kappa[p];        // thresholds for positive categories

  // Hyperparameters
  real mu_beta_pl;
  real<lower=0> sigma_beta_pl;
  real mu_kappa;
  real<lower=0> sigma_kappa;
  real<lower=0> sigma_log_alpha_pl;
  real<lower=0> sigma_log_alpha_grm;
}

transformed parameters {
  // Discriminations in original (positive) scale
  vector[p] alpha_pl  = exp(log_alpha_pl);
  vector[p] alpha_grm = exp(log_alpha_grm);

  // Precompute scaled GRM cutpoints: preserves order since alpha_grm[j] > 0
  vector[k - 2] cut[p];
  for (j in 1:p) {
    cut[j] = alpha_grm[j] * kappa[j];
  }
}

model {
  // Priors on latent traits
  Lcorr ~ lkj_corr_cholesky(2);                 // mild regularization
  theta ~ multi_normal_cholesky(Zero, Lcorr);   // unit variances by construction

  // Hyperpriors (more stable than half-Cauchy)
  mu_beta_pl            ~ normal(0, 1.5);
  sigma_beta_pl         ~ student_t(3, 0, 1);   // half-t
  mu_kappa              ~ normal(0, 1.5);
  sigma_kappa           ~ student_t(3, 0, 1);
  sigma_log_alpha_pl    ~ student_t(3, 0, 1);
  sigma_log_alpha_grm   ~ student_t(3, 0, 1);

  // Item priors
  beta_pl        ~ normal(mu_beta_pl, sigma_beta_pl);
  log_alpha_pl   ~ normal(0, sigma_log_alpha_pl);
  log_alpha_grm  ~ normal(0, sigma_log_alpha_grm);
  for (j in 1:p) kappa[j] ~ normal(mu_kappa, sigma_kappa);

  // Likelihood
  for (i in 1:N) {
    real th1 = theta[i, 1];
    real th2 = theta[i, 2];
    for (j in 1:p) {
      real eta_pl  = alpha_pl[j] * (th1 - beta_pl[j]);   // 2PL logit
      if (y[i, j] == 0) {
        target += bernoulli_logit_lpmf(0 | eta_pl);
      } else {
        int r = y[i, j];  // 1..k-1
        target += bernoulli_logit_lpmf(1 | eta_pl);                  // pass hurdle
        target += ordered_logistic_lpmf(r | alpha_grm[j] * th2, cut[j]); // GRM
      }
    }
  }
}