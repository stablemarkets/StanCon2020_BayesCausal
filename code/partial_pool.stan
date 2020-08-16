data {
  int N; // number of observations
  int Pv; // number of strata
  int Pw; // #cols of confounder matrix (including constant)
  int V[N]; // strata indicators
  matrix[N, Pw] W; // matrix of confounders (including constant in 1st entry)
  vector[N] A; // treatment indicator
  int Y[N]; // outcome
  int n_v[Pv] ; // sample size in each stratum
  int ind[Pv+1];
}

parameters {
  vector[Pw] beta_w; // confounder effects and intercept
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
    Y[i] ~ bernoulli_logit( W[i,]*beta_w + beta_v[ V[i] ] + theta[ V[i] ]*A[i] );
  }
}

// transform thetas to form Psi's (the dose-response curve)
generated quantities {

  real marg_mean_y1;
  real marg_mean_y0;
  real odds_1;
  real odds_0;
  
  vector[N] overall_cond_mean_y1;
  vector[N] overall_cond_mean_y0;
  vector[N] overall_bb_weights;
  real overall_odds_ratio;
  
  // cycle through strata of interest and compute causal Odds Ratio for each.
  vector[Pv] odds_ratio;
  
  for( v in 1:Pv ){
    
    // compute mean under strata v, under both A=1 and A=0
    int nv = n_v[v];
    int v_start = ind[v]+1;
    int v_end = ind[v+1];
    
    vector[nv] cond_mean_y1;
    vector[nv] cond_mean_y0;
    vector[nv] bb_weights;
  
    
    matrix[nv, Pw] Wv = W[ v_start:v_end,  ];
    
    // compute conditional means.
    cond_mean_y1 = inv_logit( Wv*beta_w + beta_v[v] + theta[v] );
    cond_mean_y0 = inv_logit( Wv*beta_w + beta_v[v] );
    
    // Bayesian bootstrap weights
    bb_weights = dirichlet_rng( rep_vector(1, nv) ) ;
    
    // taking average over bayesian bootstrap weights
    marg_mean_y1 = bb_weights' * cond_mean_y1;
    marg_mean_y0 = bb_weights' * cond_mean_y0;
    
    // compute odds ratio 
    odds_1 = (marg_mean_y1/(1 - marg_mean_y1));
    odds_0 = (marg_mean_y0/(1 - marg_mean_y0));
    odds_ratio[v] = odds_1/odds_0;
  }
  
  for(i in 1:N){
    overall_cond_mean_y1[i] = inv_logit( W[i, ]*beta_w + beta_v[V[i]] + theta[ V[i] ] );
    overall_cond_mean_y0[i] = inv_logit( W[i, ]*beta_w + beta_v[V[i]] );
  }
  
  overall_bb_weights = dirichlet_rng( rep_vector(1, N ) ) ;
  marg_mean_y1 = overall_bb_weights' * overall_cond_mean_y1 ;
  marg_mean_y0 = overall_bb_weights' * overall_cond_mean_y0 ;
  
  // compute odds ratio
  overall_odds_ratio = (marg_mean_y1/(1 - marg_mean_y1)) / (marg_mean_y0/(1 - marg_mean_y0));

}
