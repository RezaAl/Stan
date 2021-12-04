In_data <- list(P,D,lX,N,M,aux1,aux2,max_tau,indx,y,tau,tx,X); #previously prepared based on the input data

fit0 <- stan(file = "<Path_to_Stan_file>",
              data = In_data, iter = itr, warmup=floor(0.5*itr), chains = 4, cores=4,
init = list(list(trend=5, trend_p = rep(6,P), beta_h = rep(0.05,D),beta_h_p = rep(0.05,P),beta_l = rep(0,D),tx_star = tx,alpha = c(22),alpha_p = rep(22,P),sig_y=rep(0.50,P),sig_time=c(30)),
            list(trend=5.25, trend_p = rep(5.25,P), beta_h = rep(0.04,D),beta_h_p = rep(0.04,P),beta_l = rep(0.01,D),tx_star = tx+3,alpha = c(20),alpha_p = rep(20,P),sig_y=rep(0.55,P),sig_time=c(34)),
            list(trend=4.5, trend_p = rep(5,P), beta_h = rep(0.03,D),beta_h_p = rep(0.03,P),beta_l = rep(-0.01,D),tx_star = tx-5,alpha = c(18),alpha_p = rep(18,P),sig_y=rep(0.63,P),sig_time=c(37)),
            list(trend=5.5, trend_p = rep(5.5,P), beta_h = rep(0.02,D),beta_h_p = rep(0.02,P),beta_l = rep(0.02,D),tx_star = tx+5,alpha = c(16),alpha_p = rep(16,P),sig_y=rep(0.65,P),sig_time=c(40))));

samples_alpha_before_RYGB <- extract(fit0, permuted = TRUE) 

saveRDS(samples_alpha_RYGB, "samples_alpha_RYGB.rds")
