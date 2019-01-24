data {
  int<lower=1> maxm;
  int<lower=1> maxJ;
  int<lower=1> summ;
  int<lower=1> sumJ;
  int<lower=1> n;                    //number of individuals
  int<lower=1> D;                    //number of subtests (and dimensions)
  int<lower=1> J[D];                 //number of subtests' items
  int m[maxJ, D];                    //number of response categories
  int<lower=1> lenY;                 //number of observations (elements of Y)
  int<lower=1, upper=maxm> Y[lenY];  //matrix of responses
  int a_idx[maxJ, D];                //index of discriminations parameters
  int b_idx[maxJ, D, maxm-1];        //index of difficult parameters
  int theta0_idx[n];                 //index of global individuals parameters
  int theta_idx[D, n];               //index of specific individuals parameters
  int Y_idx[n, maxJ, D];             //index of response matrix and probs
}

parameters {
  real<lower=0> a0[sumJ];            //global discrimination parameters
  real<lower=0> a[sumJ];             //specific discrimination parameters
  real b[summ-sumJ];                 //difficulty parameters
  real theta[n*D];                   //specific latent traits
  real theta0[n];                    //global latent trait
}

transformed parameters {
  real eta[n, maxJ, D, maxm];
  vector[maxm] psum[n, maxJ, D];
  real b_aux;
  simplex[maxm] prob[lenY];
  for (i in 1:n) {
    for (d in 1:D) {
      for (j in 1:J[d]) {
        for (k in 1:m[j, d]) {
          if (k == 1) {
            b_aux <- 0;
          } else {
            b_aux <- b[b_idx[j, d, k-1]];
          }
		  eta[i, j, d, k] <- a0[a_idx[j, d]]*theta0[theta0_idx[i]] + a[a_idx[j, d]]*theta[theta_idx[d, i]] - b_aux;
          psum[i, j, d, k] <- sum(eta[i, j, d, 1:k]);
		}
		if (Y_idx[i, j, d] != 0){
		  prob[Y_idx[i, j, d], 1:m[j, d]] <- softmax(psum[i, j, d, 1:m[j, d]]);
		}
      }
    }
  }
}

model {
  for (i in 1:n) {
    theta0[theta0_idx[i]] ~ normal(0, 1);     //prior distribution for global latent traits
	for (d in 1:D)
	  theta[theta_idx[d, i]] ~ normal(0, 1);  //prior distribution for specific latent traits
  }
  for (d in 1:D) {
    for (j in 1:J[d]) {
      a0[a_idx[j, d]] ~ normal(0,1) T[0,];    //prior distribution for global discrimination parameters
	  a[a_idx[j, d]] ~ normal(0,1) T[0,];     //prior distribution for specific discrimination parameters
	  for (k in 1:(m[j, d]-1))
	    b[b_idx[j, d, k]] ~ normal(0, 1);     //prior distribution for difficulty parameters
	}
  }
  for (l in 1:lenY) {
    Y[l] ~ categorical(prob[l, 1:maxm]);
  }
}

generated quantities {
  int<lower=1, upper=maxm> y_rep[lenY];
  for (l in 1:lenY) {	
    y_rep[l] <- categorical_rng(prob[l, 1:maxm]);
  }
}
