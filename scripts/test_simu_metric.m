
% default parameters generation
cfg = struct('nChan',100,'nMat',500,'p',2,'nCommonNoiseMat',10,'cond_sig',30,'cond_noise',[0,inf],...
    'SNR_white',0.001,'SNR_struct',2);

[C,A,info] = simRandDat(cfg);
Geuc=mean(C,3)

disp('logdet')
[Gld, critere, niter]=logdet_mean(C);critere

disp('wasserstein')
tic ;[Gws, critere1, niter1]=wasserstein_mean(C);critere;toc
tic;[Gws2, critere2, niter2]=opttransp_mean(C);critere;toc
tic;[Gws3, critere3, niter3]=ws_congedo_mean(C);critere3;toc
tic;[Gws2, critere2, niter4]=opttransp_mean(C);critere;toc


distance_opttransp(Gws,Gws2)
distance_ws(Gws,Gws2)
distance_ws(Gws2,Gws3)

% distance_ha(Gld,Gws)
distance_ld(Gld,Gws)
distance_riemann(Gld,Gws)


% inv((inv(Gld)+inv(Gws))/2)

distance_ws(C(:,:,1),C(:,:,1))

N=200;A=randn(N)*randn(N,1)*100;sqrt(norm((A*A')/sqrt(2)/N^(3/2),'fro'))
N=2000;A=randn(N)*20;B=randn(N)*30;sqrt(norm((A*A'-B*B')/sqrt(2)/N^(3/2),'fro'))

N=1000;A=randn(N)*[100]';(norm(A*A','fro'))^(2/3)/N

N*(N)/2