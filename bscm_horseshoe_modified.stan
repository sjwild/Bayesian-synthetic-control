data{
  
  // Dimensions
  int N_pre;                                // Number of observations in the pre-treatment periods
  int N_post;                               // Number of observations in the post-treatment periods
  int p;                                    // Number of predictors (control units)
  
  // scales
  real<lower=0> scale_alpha;               // prior std for the intercept
  real<lower=0> scale_g;                   // scale for the half-t prior for tau
  
  // slabs
  real<lower=0> slab_scale;                // slab scale for the regularized horseshoe
  real<lower=0> slab_df;                   // slab degrees of freedom for the regularized

  // nu for degrees of freedom
  real<lower=1> nu_g;                      // degrees of freedom for the half-t prior
  real<lower=1> nu_l;                      // degrees of freedom for the half-t priors

  // inputs
  matrix[N_pre, p] X_pre;                  // design matrix pre-treatment (no intercept)
  matrix[N_post, p] X_post;                // design matrix in the post-treatment period (no intercept) 
  
  // y
  real y_pre[N_pre];                        // Treated unit in the pre-treatment periods

}

parameters{
  
  // alpha and beta
  real alpha;  
  vector[p] beta_raw;                       // Control unit weights (will be transformed)
  
  // shrinkage parameters
  real<lower=0> tau;                        // Global shrinkage
  vector<lower=0>[p] lambda;                // Local shrinkage
  
  // sigma
  real logsigma;
  
  //other parameters
  real<lower=0> caux;

}

transformed parameters{
  
  // containers
  vector[p] beta;                                // Control unit weights
  vector<lower=0>[p] lambda_tilde;               // Local shrinkage
  vector[N_pre] mu;                              // y_hat in the pre-treatment period
  real<lower=0> sigma;                           // error term for regression
  real<lower=0> c;                               // slab scale

  
  // transformed parameters
  sigma = exp(logsigma);
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(square(c) * square(lambda) ./ (square(c) + square(tau) * square(lambda)));
  beta = beta_raw .* lambda_tilde * tau;
  mu = alpha + X_pre * beta;
  
}

model{
  
  // priors
  beta_raw ~ std_normal();
  lambda ~ student_t(nu_l, 0, 1);
  tau ~ student_t(nu_g, 0, scale_g * sigma);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  alpha ~ normal(0, 1);

  // model
  y_pre ~ normal(mu, sigma);
  
}

generated quantities{
  //Post-treatment prediction & Log-likelihood
  vector[N_pre] y_fit;                        //Fitted synthetic control unit in the pre-treatment
  vector[N_post] y_post;                      //Predicted synthetic control unit in the post-treatment
  vector[N_pre] log_lik;                      //Log-likelihood
  for(f in 1:N_pre){
    y_fit[f] = normal_rng(alpha + X_pre[f,] * beta, sigma);
  }

  for(i in 1:N_post){
    y_post[i] = normal_rng(alpha + X_post[i,] * beta, sigma);
  }
  
  for (t in 1:N_pre){
    log_lik[t] = normal_lpdf(y_pre[t] | y_fit[t], sigma);
  }
  
}

