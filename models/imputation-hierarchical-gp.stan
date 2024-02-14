
data {
  int<lower=1> N; // number of observations
  int<lower=1> N_sites; // number of sites
  int<lower=1> N_countries; // number of countries
  int<lower=1> N_regions; // number of regions
  int<lower=1> N_times; // number of time points 
  int<lower=1> N_years; // number of years
  int<lower=1,upper=N_sites> site_ind[N]; // site index
  int<lower=1,upper=N_regions> region_ind[N]; // region index
  int<lower=1,upper=N_sites> site_region_ind[N_sites]; // site region index
  int<lower=1,upper=N_countries> country_ind[N]; // country index
  int<lower=1,upper=N_sites> site_country_ind[N_sites]; // site country index
  int<lower=1> time_ind[N]; // time index
  int<lower=1> year_ind[N]; // year index
  int<lower=1> quarter_ind[N]; // quarter index
  int<lower=1,upper=N_years> time_year_ind[N_times]; // time year index
  int<lower=1,upper=4> time_quart_ind[N_times]; // time quarter index
  vector<lower=-1>[N] y; // quarterly number of cases (-1: NA, >=0: quarterly cases or annual cases divided by 4)
  int<lower=-1,upper=1> missing[N]; // outcome missing (-1: neither quarterly nor annual outcome, 0: only annual outcome, 1: quarterly outcome)
  int<lower=0> missing_ind[N]; // missing index
  int<lower=0> K_missing; // number missing completely
  int<lower=0> J_missing_quart; // number missing only quarterly
  int<lower=0> J_missing_years; // total number of missing years
  int<lower=0> missing_quart_group_ind[N]; // indicator of the missing year group
  vector<lower=-1>[J_missing_years] y_annual; // annual number of cases for years with missing quarterly data
}

transformed data {
  real time[N_times + 4];
  vector[2 + N_regions + N_countries + 4] counts;
  int n_comps;
  int N_fcat;
  int N_times_fcat;
  
  N_fcat = N_sites * 4;
  N_times_fcat = N_times + 4;
  
  for (t in 1:N_times_fcat) 
    time[t] = t;
  
  n_comps = rows(counts);
  
  for (i in 1:n_comps)
    counts[i] = 2;
}

parameters {
  vector<lower=0>[K_missing] log_y_missing;
  simplex[4] alpha_missing_quart[J_missing_years];
  vector<lower=0>[J_missing_years] lambda_missing_quart;
  simplex[4] w_missing_quart[J_missing_years];
  matrix[N_times_fcat,N_sites] GP_site_short_std;
  matrix[N_times_fcat,N_sites] GP_site_long_std;
  matrix[N_times_fcat,N_regions] GP_region_short_std;
  matrix[N_times_fcat,N_regions] GP_region_long_std;
  real<lower=1,upper=N_times> length_GP_site_short;
  real<lower=1,upper=N_times> length_GP_site_long;
  real<lower=1,upper=N_times> length_GP_region_short;
  real<lower=1,upper=N_times> length_GP_region_long;
  real<lower=0> tot_var;
  simplex[n_comps] prop_var;
  vector[N_sites] site_std;
  vector[N_countries] country_std;
  vector[N_regions] region_std;
  real alpha;
}

transformed parameters {
  vector[4] y_missing_quart[J_missing_years];
  vector[N] log_y_imp;
  vector[N] mu;
  real sigma_region;
  real sigma_country;
  vector[N_regions] sigma_site_region;
  vector[N_countries] sigma_site_country;
  real sigma_GP_site_short;
  real sigma_GP_site_long;
  real sigma_GP_region_short;
  real sigma_GP_region_long;
  vector[n_comps] vars;
  vector[N_sites] site_re;
  vector[N_countries] country_re;
  vector[N_regions] region_re;
  matrix[N_times_fcat,N_sites] GP_site_short;
  matrix[N_times_fcat,N_sites] GP_site_long;
  matrix[N_times_fcat,N_regions] GP_region_short;
  matrix[N_times_fcat,N_regions] GP_region_long;

  for (j in 1:J_missing_years) {
    y_missing_quart[j] = w_missing_quart[j] * y_annual[j];
  }
  
  // transformed outcome with missing data
  {
    for (i in 1:N) {
      if (missing[i] < 0) {
        log_y_imp[i] = log(exp(log_y_missing[missing_ind[i]]) + .25); 
      } else if (missing[i] > 0) {
        log_y_imp[i] = log(y_missing_quart[missing_quart_group_ind[i],quarter_ind[i]] + .25);
      } else {
        log_y_imp[i] = log(y[i] + .25);
      }
    }
  }
  
  // variances and random effects
  vars = n_comps * prop_var * tot_var;
  sigma_region = sqrt(vars[1]);
  sigma_country = sqrt(vars[2]);
  for (i in 1:N_regions)
    sigma_site_region[i] = sqrt(vars[2 + 4 + i]);
  for (i in 1:N_countries)
    sigma_site_country[i] = sqrt(vars[2 + N_regions + i]);
  sigma_GP_site_short = sqrt(vars[2 + N_regions + N_countries + 1]);
  sigma_GP_site_long = sqrt(vars[2 + N_regions + N_countries + 2]);
  sigma_GP_region_short = sqrt(vars[2 + N_regions + N_countries + 3]);
  sigma_GP_region_long = sqrt(vars[2 + N_regions + N_countries + 4]);
  region_re = sigma_region * region_std;
  country_re = sigma_country * country_std;
  site_re = (sigma_site_region[site_region_ind] + sigma_site_country[site_country_ind]) .* site_std; 

  
  // site- and regional GPs
  {
    matrix[N_times_fcat, N_times_fcat] cov_site_short;
    matrix[N_times_fcat, N_times_fcat] cov_site_long;
    matrix[N_times_fcat, N_times_fcat] cov_region_short;
    matrix[N_times_fcat, N_times_fcat] cov_region_long;
    matrix[N_times_fcat, N_times_fcat] L_cov_site_short;
    matrix[N_times_fcat, N_times_fcat] L_cov_site_long;
    matrix[N_times_fcat, N_times_fcat] L_cov_region_short;
    matrix[N_times_fcat, N_times_fcat] L_cov_region_long;
    cov_site_short = cov_exp_quad(time, sigma_GP_site_short, length_GP_site_short);
    cov_site_long = cov_exp_quad(time, sigma_GP_site_long, length_GP_site_long);
    cov_region_short = cov_exp_quad(time, sigma_GP_region_short, length_GP_region_short);
    cov_region_long = cov_exp_quad(time, sigma_GP_region_long, length_GP_region_long);
    for (t in 1:N_times_fcat) {
      cov_site_short[t, t] = cov_site_short[t, t] + 1e-12;
      cov_site_long[t, t] = cov_site_long[t, t] + 1e-12;
      cov_region_short[t, t] = cov_region_short[t, t] + 1e-12;
      cov_region_long[t, t] = cov_region_long[t, t] + 1e-12;
    }
    L_cov_site_short = cholesky_decompose(cov_site_short);
    L_cov_site_long = cholesky_decompose(cov_site_long);
    L_cov_region_short = cholesky_decompose(cov_region_short);
    L_cov_region_long = cholesky_decompose(cov_region_long);
    GP_site_short = L_cov_site_short * GP_site_short_std;
    GP_site_long = L_cov_site_long * GP_site_long_std;
    GP_region_short = L_cov_region_short * GP_region_short_std;
    GP_region_long = L_cov_region_long * GP_region_long_std;
  }
  
  for (i in 1:N) {
    mu[i] = alpha 
      + site_re[site_ind[i]] 
      + region_re[region_ind[i]]
      + country_re[country_ind[i]]
      + GP_site_short[time_ind[i], site_ind[i]]
      + GP_site_long[time_ind[i], site_ind[i]]
      + GP_region_short[time_ind[i], region_ind[i]]
      + GP_region_long[time_ind[i], region_ind[i]];
  }
  
}



model {

  // priors
  lambda_missing_quart ~ exponential(0.1);
  for (j in 1:J_missing_years) {
    alpha_missing_quart[j] ~ dirichlet(rep_vector(1., 4));
    w_missing_quart[j] ~ dirichlet(alpha_missing_quart[j] * lambda_missing_quart[j]);
  }

  site_std ~ normal(0, 1);
  region_std ~ normal(0, 1);
  country_std ~ normal(0, 1);
  
  tot_var ~ gamma(5, 5);
  prop_var ~ dirichlet(counts);
  
  to_vector(GP_site_short_std) ~ normal(0, 1);
  to_vector(GP_site_long_std) ~ normal(0, 1);
  to_vector(GP_region_short_std) ~ normal(0, 1);
  to_vector(GP_region_long_std) ~ normal(0, 1);
  
  length_GP_site_short ~ weibull(32, 2);
  length_GP_site_long ~ weibull(32, 8);
  length_GP_region_short ~ weibull(32, 2);
  length_GP_region_long ~ weibull(32, 8);
  
  alpha ~ normal(0, 1);
  
  // likelihood
  log_y_imp ~ normal(mu, tot_var); 
}

generated quantities {
  vector[N+N_fcat] y_pred;
  
  // model fit and imputation
  for (i in 1:N) 
    y_pred[i] = poisson_rng(fmin(fmax(1e-5, exp(mu[i])-.25), 1e9));
    
  // model forecast 
  {
    int k = N;
    real mu_fcat;
    for (s in 1:N_sites) {
      for (t in (N_times+1):(N_times+4)) {
        k = k + 1;
        mu_fcat = alpha 
          + site_re[s] 
          + region_re[site_region_ind[s]]
          + country_re[site_country_ind[s]]
          + GP_site_short[t,s]
          + GP_site_long[t,s]
          + GP_region_short[t,site_region_ind[s]]
          + GP_region_long[t,site_region_ind[s]];
        y_pred[k] = poisson_rng(fmin(fmax(1e-5, exp(mu_fcat)-.25), 1e9));
      }
    }
  }
}
