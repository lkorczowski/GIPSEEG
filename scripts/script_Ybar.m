clear all
load('Paolo_echantillons.mat')
Ybar.wasserstein=mean_covariances(Y,'riemann')
Ybar.logdet=mean_covariances(Y,'logdet')
Ybar.logdet=mean_covariances(Y,'logdet')
mean_param=0.1
i_cfg = struct('tol',10^(-12),'N_itermax',100,'t',0.1,'phi',0.8);

[Ybar.power,crit.power,niter.power] = mean_power(Y,i_cfg)
[Ybar.power_dual,crit.power_dual,niter.power_dual] = mean_power_dual(Y,i_cfg)
Ybar.logeuc = logeuclid_mean(Y);
Ybar.ws = mean_wasserstein(Y);
[Ybar.logdet , crit.logdet, niter.logdet] = logdet_mean(Y,i_cfg.tol,i_cfg.N_itermax);
[Ybar.fisher, crit.fisher, niter.fisher] = riemann_mean(Y,i_cfg.tol,i_cfg.N_itermax,[],100);
[Ybar.fisher_c, crit.fisher_c, niter.fisher_c] = riemann_mean(Y,i_cfg.tol,i_cfg.N_itermax,[],1);


plot(crit.fisher)
save('Ybar.mat','Ybar','crit','niter')