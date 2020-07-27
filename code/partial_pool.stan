data {
  int N; // number of observations
  int Pv; // number of strata
  int Pw; // #cols of confounder matrix (including constant)
  int V[N]; // strata indicators
  real W[N]; // matrix of confounders (including constant in 1st entry)
  vector[N] A; // treatment indicator
  int Y[N]; // outcome
}

parameters {
  real beta_w; // confounder effects and intercept
  vector[Pv] beta_v; // strata main effects
  vector[Pv] theta; // interaction effects
  real mu;
}

model {
  
  // specify priors
  
  beta_w ~ normal(0, 1); 
  beta_v ~ normal(0, 1); 
  
  mu ~ normal(0, 1);
  theta ~ normal( mu, .5);
  
  // specify likelihood 
  for(i in 1:N){
    Y[i] ~ bernoulli_logit( W[i]*beta_w + beta_v[ V[i] ] + theta[ V[i] ]*A[i] );
  }
}

// transform thetas to form Psi's (the dose-response curve)
generated quantities {
  
// take draw of bayesian bootstrap weights
vector[N] bb_weights = dirichlet_rng( rep_vector( 1, N) ) ;

vector[N] cond_mean_y1;
vector[N] cond_mean_y0;

real marg_mean_y1;
real marg_mean_y0;

real odds_1;
real odds_0;

vector[Pv] odds_ratio;

real overall_odds_ratio;

// cycle through strata of interest and compute causal Odds Ratio for each.
for( v in 1:Pv ){

  // compute mean under strata v, under both A=1 and A=0
  for(i in 1:N){
    cond_mean_y1[i] = inv_logit( W[i]*beta_w + beta_v[v] + theta[v] );
    cond_mean_y0[i] = inv_logit( W[i]*beta_w + beta_v[v] );
  }

  // taking average over bayesian bootstrap weights
  marg_mean_y1 = bb_weights' * cond_mean_y1 ;
  marg_mean_y0 = bb_weights' * cond_mean_y0 ;

  // compute odds ratio 
  odds_ratio[v] = (marg_mean_y1/(1 - marg_mean_y1)) / (marg_mean_y0/(1 - marg_mean_y0));
}


for(i in 1:N){
  cond_mean_y1[i] = inv_logit( W[i]*beta_w + beta_v[V[i]] + mu );
  cond_mean_y0[i] = inv_logit( W[i]*beta_w + beta_v[V[i]] );
}

marg_mean_y1 = bb_weights' * cond_mean_y1 ;
marg_mean_y0 = bb_weights' * cond_mean_y0 ;

// compute odds ratio 
overall_odds_ratio = (marg_mean_y1/(1 - marg_mean_y1)) / (marg_mean_y0/(1 - marg_mean_y0));

}
