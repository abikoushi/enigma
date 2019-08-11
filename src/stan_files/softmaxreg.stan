data{
  int N; //num of record
  int K; //num of category
  int L; //num of topics
  int D; //num of covariate
  int y[N,K];
  matrix[N,D] X;
  vector<lower=0>[L] alpha;//hyper parameter of pi
}
parameters{
  simplex[L] pi;
  matrix[D,K] beta;
  matrix[K,L] gamma;
  real<lower=0> sigma; //hyper parameter of beta
  real<lower=0> tau; //hyper parameter of gamma
}
transformed parameters{
  matrix[N,K] Xbeta;
  Xbeta = X*beta;
}
model{
  for(d in 1:D){
   target += normal_lpdf(beta[d,]|0,sigma);
  }
  for(l in 1:L){
    target += normal_lpdf(gamma[,l]|0,tau);
  }
  target += dirichlet_lpdf(pi|alpha);
  for(n in 1:N){
    vector[L] lp;
    for(l in 1:L){
      lp[l] = log(pi[l]) + multinomial_lpmf(y[n,]|softmax(gamma[,l]+Xbeta[n,]'));
    }
    target += log_sum_exp(lp);
  }
}
generated quantities{
  simplex[L] cate[N];
  for(n in 1:N){
    vector[L] lp;
    for(l in 1:L){
      lp[l] = categorical_lpmf(l|pi) + multinomial_lpmf(y[n,]|softmax(gamma[,l]+Xbeta[n,]'));
    }
    cate[n] = softmax(lp);
  }
}
