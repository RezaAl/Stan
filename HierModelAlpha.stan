 data { 
int<lower=1> P; // # of individuals  
int<lower=1> D; //dimension of the meals' covariates 
int<lower=1> lX; //dimension of the meals' covariates 
int<lower=0> N[P]; 
int<lower=0> M[P];  
int<lower=1> aux1; 
int<lower=1> aux2; 
int<lower=1> max_tau; 
int indx[P, aux2]; // covariates' value 
matrix[P, aux1] y;  //observation 
matrix[P, aux1] tau;  //timestamp of observations 
vector[lX] tx; //reported time of meals   
//matrix[lX,D] X; // covariates' value 
vector[lX] X; // covariates' value 
}   

parameters { 
vector<lower=1,upper= max_tau>[lX] tx_star; // true time of meals  
real<lower=0> sig_time; // std of the time uncertainty  
real<lower=0> sig_y[P];   // std of the obs. likelihood  
real<lower=0>  beta_h; // 
real<lower=0> beta_h_p[P];   
real<lower=0> alpha; 
real<lower=0> alpha_p[P];  
//real<lower=0> trend; 
real<lower=0> trend_p[P]; 

} 

model {  
real fx[aux2]; 
fx=rep_array(0, aux2); 
beta_h ~ normal(0.15,0.1); //0,0.5
sig_time ~ normal(0,10);  
tx ~ normal(tx_star, sig_time);  
tx_star ~ uniform(1, max_tau);  
alpha ~ normal(20,10);  //30 
//trend ~ normal(3.5,0.25);  //30  

for (p in 1:P){ 

sig_y[p] ~ normal(0,0.5);  
alpha_p[p] ~ normal(alpha,10); 
beta_h_p[p] ~ normal(beta_h,0.05); 
//trend_p[p] ~ normal(trend,0.5);  
trend_p[p] ~ normal(3.5,1); 

  for (n in 1:N[p]){  

    for (m in 1:M[p]) { 

      if (fabs(tx[indx[p,m]]-tau[p,n])<150)  
      { 	
        fx[m] = (X[indx[p,m]]*beta_h_p[p]')*exp(-0.5*(tau[p,n]-tx_star[indx[p,m]]-3*(alpha_p[p]))^2/(alpha_p[p])^2); }  
        } 
    y[p,n] ~ normal(sum(fx)+trend_p[p], sig_y[p]); 
    } 
} 
} 

generated quantities{ 
real R[P, aux1, aux2]; 
matrix[P, aux1] Rout; 
matrix[P, aux1] log_lik1; 
row_vector[aux1] log_lik;  
log_lik1= rep_matrix(0,P, aux1); 
log_lik= rep_row_vector(0, aux1); 
R=rep_array(0,P, aux1, aux2); 
Rout=rep_matrix(0,P, aux1);  
  
for (p in 1:P){
  for (n in 1:N[p]){  
    for (m in 1:M[p]) { 
      if (fabs(tx[indx[p,m]]-tau[p,n]) < 150) {
        R[p,n,m] = (X[indx[p,m]]*beta_h_p[p]')*exp(-0.5*(tau[p,n]-tx_star[indx[p,m]]-3*(alpha_p[p]))^2/(alpha_p[p])^2); }  
        } 
    Rout[p,n] =sum(R[p,n,]); 
    log_lik1[p,n]= normal_lpdf(y[p,n] | Rout[p,n]+trend_p[p], sig_y[p]);  //Rout[p,n]+trend_p[p]
    } 
log_lik+=log_lik1[p]; 
} 
} 