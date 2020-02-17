%% Identification of strictly causal MVAR model and estimation of frequency domain causality
% REFERENCE: Luca Faes and Giandomenico Nollo (2011). Multivariate Frequency Domain Analysis of Causal Interactions in Physiological Time Series,
% Biomedical Engineering, Trends in Electronics, Communications and Software, Anthony N. Laskovski (Ed.), ISBN: 978-953-307-475-7, InTech,  
% Available from: http://www.intechopen.com/articles/show/title/multivariate-frequency-domain-analysis-of-causal-interactions-in-physiological-time-series 

% Here we generate a realization of Eq. 18, then we identify model coefficients
% and estimate frequency domain causality, comparing to theoretical functions

clear; close all; clc;
nfft=512; % number of frequency bins
fc=1; % sample frequency
T=1000; % length of simulated series
pcrit='aic'; % 'aic', 'mdl', or number for fixed model order
set(0,'RecursionLimit',10000)
do_surrogates='n'; %'y' performs surrogates, otherwise do not perform
numsurr=10; %number of surrogate time series
upalpha=95; %upper percentile of threshold

simunumber=1; % which simulation (either 1 or 4 are strictly causal)

idMode=0; %identification algorithm (0:Least Squares covariance; 1:Yule-Walker correlation; 2:Vieira-Morf partial correlation; 3:Nutall-Strandpartial correlation; other: see mvar.m)

%% alternative init
N=3;
pr=10;
A=zeros(N,N*pr);
    x0=zeros(pr,N);

% indP=;in=;ou=;A(ou,in+N*(indP-1))=;%at the order indP add the link from in to out
indP=1;in=1;ou=1;A(ou,in+N*(indP-1))=0.5;%at the order indP add the link from in to out
indP=1;in=2;ou=2;A(ou,in+N*(indP-1))=0.5;%at the order indP add the link from in to out
indP=1;in=3;ou=3;A(ou,in+N*(indP-1))=0.5;%at the order indP add the link from in to out

indP=3;in=1;ou=3;A(ou,in+N*(indP-1))=0.5;%at the order indP add the link from in to out
indP=10;in=3;ou=2;A(ou,in+N*(indP-1))=0.6;%at the order indP add the link from in to out
indP=5;in=3;ou=1;A(ou,in+N*(indP-1))=0.1;%at the order indP add the link from in to out

% indP=12;in=3;ou=1;A(ou,in+N*(indP-1))=0.2;%at the order indP add the link from in to out
% indP=13;in=3;ou=1;A(ou,in+N*(indP-1))=0.1;%at the order indP add the link from in to out

% indP=10;in=1;ou=2;A(ou,in+N*(indP-1))=0.6;%at the order indP add the link from in to out

% indP=1;in=1;ou=2;A(ou,in+N*(indP-1))=1.5;%at the order indP add the link from in to out
% indP=1;in=1;ou=3;A(ou,in+N*(indP-1))=1;%at the order indP add the link from in to out
B=eye(N);
C=eye(N);
D=eye(N)*0.1;
X=reshape(x0,pr,N);
Y=[];
w=randn(T,N);
u=zeros(T,N);
u(rand(T*N,1)<5e-2)=1;

%% Theoretical MVAR coefficients and spectral functions
% % coefficients of extended MVAR model
% [Bm,B0,Sw]=simuMVARcoeff(simunumber);
Bm=A;B0=zeros(N);Sw=C;
[Am, Su]=diag_coeff_rev(Bm,B0,Sw);
M=size(Bm,1);

%%% Theoretical spectral functions
[dc,dtf,pdc,gpdc,coh,pcoh,pcoh2,h,s,pp,f] = fdMVAR(Am,Su,nfft,fc);
S=abs(s); % spectral matrix
P=abs(pp); % inverse spectral matrix
DC=abs(dc).^2; % directed coherence
PDC=abs(gpdc).^2; % partial directed coherence
COH=abs(coh).^2; %coherence
PCOH=abs(pcoh).^2; % partial coherence
PCOH2=abs(pcoh2).^2; % partial coherence, other definition to check equality

%% generation of simulated time series
% U=randn(M,T); % uncorrelated gaussian innovations
U=InstModelfilter(T,Su,'StrictlyCausal'); % gaussian innovations with covariance Su
U=u';
[Y]=MVARfilter(Am,U); % realization of Eq. 18

%% Estimated MVAR coefficients and spectral functions 
%model order selection
if pcrit(1)=='a' | pcrit(1)=='m' 
    [pottaic,pottmdl,aic,mdl] = mos_idMVAR(Y,20,idMode);
    if pcrit(1)=='a', p=pottaic; else p=pottmdl; end
else
    p=pcrit;
end
% p=9

% model identification
[eAm,eSu,Yp,Up]=idMVAR(Y,p,idMode);

%%% Estimated spectral functions
[dc2,dtf2,pdc2,gpdc2,coh2,pcoh2,pcoh22,h2,s2,pp2] = fdMVAR(eAm,eSu,nfft,fc);
eS=abs(s2); % spectral matrix
eP=abs(pp2); % inverse spectral matrix
eDC=abs(dc2).^2; % directed coherence
ePDC=abs(gpdc2).^2; % partial directed coherence
eCOH=abs(coh2).^2; %coherence
ePCOH=abs(pcoh2).^2; % partial coherence

%% testing model assumptions
alpha=0.05; %significance for the tests

%%% WHITENESS: portmanteau Ljung-Box test
Up=Y-Yp; Up=Up(:,p+1:T); % residuals of strictly causal model
h=50; %number of lags
[pval,Qh,critlo,crithi,crit1tail,C0,stringWhite,flagWhite]=test_whiteness(Up,p,h,alpha);

%%% INDEPENDENCE among U residuals
testtype='Kendall';
[pM,rho,contarigett,stringInd]=test_independence(Up,testtype,alpha);

%%% NONGAUSSIANITY of U residuals (NOT REQUIRED HERE!)
[pGauss,ps,pk,lambdas,lambdak,crittresh,stringGauss,flagGauss]=test_gaussianity(Up,alpha);

%% surrogate time series - significance thresholds
if do_surrogates=='y'
    eCOHs=NaN*ones(M,M,nfft,numsurr); ePCOHs=eCOHs; eDCs=eCOHs; ePDCs=eCOHs;
    for is=1:numsurr 
        disp(['surrogate ' int2str(is) ' of ' int2str(numsurr)])
        Ys=surrVFT(Y); % FT surrogates
        [eAms,eSus,Yps,Ups]=idMVAR(Ys,p,idMode);
        [dcs,dtfs,pdcs,gpdcs,cohs,pcohs] = fdMVAR(eAms,eSus,nfft,fc);
        eCOHs(:,:,:,is)=abs(cohs).^2; %surro coherence
        ePCOHs(:,:,:,is)=abs(pcohs).^2; %surro partial coherence (maybe one should devise ad-hoc surrogates)
        
        % surrogates for DC and PDC
        for ii=1:M
            for jj=1:M
                if ii~=jj
                    Ys=surrVCFTf(Y,eAm,eSu,ii,jj); %CFTf surrogates
                    [eAms,eSus,Yps,Ups]=idMVAR(Ys,p,idMode);
                    [dcs] = fdMVAR(eAms,eSus,nfft,fc);
                    eDCs(ii,jj,:,is)=abs(dcs(ii,jj,:)).^2; %surro DC
                    
                    Ys=surrVCFTd(Y,eAm,eSu,ii,jj); %CFTd surrogates
                    [eAms,eSus,Yps,Ups]=idMVAR(Ys,p,idMode);
                    [dcs,dtfs,pdcs,gpdcs] = fdMVAR(eAms,eSus,nfft,fc);
                    ePDCs(ii,jj,:,is)=abs(gpdcs(ii,jj,:)).^2; %surro DC           
                end
            end
        end
        
        
    end
    eCOHsth=prctile(eCOHs,upalpha,4); ePCOHsth=prctile(ePCOHs,upalpha,4);; %thresholds
    eDCsth=prctile(eDCs,upalpha,4); ePDCsth=prctile(ePDCs,upalpha,4);

else
    eCOHsth=NaN*ones(M,M,nfft); ePCOHsth=eCOHsth;
    eDCsth=eCOHsth; ePDCsth=eCOHsth;
end



%% disps
clc;
disp(['Model order:']);
disp(['original: p=' int2str(size(Bm,2)/M) '; estimated: p=' int2str(p)]);
disp([' ']);
disp('Whiteness (Ljung-Box) test:');
disp(['p=' num2str(pval) ', ' stringWhite]);
disp([' ']);
disp(['Independence (' testtype ') test:']);
disp([int2str(contarigett) ' pairs of signals are dependent']);
disp(stringInd);
disp([' ']);
disp('Gaussianity (Jarque-Bera) test:');
disp(['p=' num2str(pGauss) ', ' stringGauss]);



%% graphs 
figure('numbertitle','off','name','Simulated time series'); % realizations of functions
for i=1:M
    subplot(M,1,i); plot(Y(i,:));
    xlim([0 T]);
    xlabel(['n']);ylabel(['y_' int2str(i)]);
end

h2a=figure('numbertitle','off','name','Spectra and Coherence (cyan:original; blue:estimated; blue-dashed:threshold)');
h2b=figure('numbertitle','off','name','Inverse spectra and Partial Coherence (magenta:original; red:estimated; red-dashed:threshold)');
h3=figure('numbertitle','off','name','Directed Coherence (cyan:original; blue:estimated; blue-dashed:threshold)');
h4=figure('numbertitle','off','name','Partial Directed Coherence (magenta:original; red:estimated; red-dashed:threshold)');
for i=1:M        
    for j=1:M
        figure(h2a); % spectra and coherence       
        subplot(M,M,(i-1)*M+j); 
        if i==j
            plot(f, squeeze(abs(S(i,j,:))),'k'); hold on; plot(f, squeeze(abs(eS(i,j,:))),'k:');
            xlim([0 fc/2]); title(['S' int2str(i)]);
        else
            plot(f, squeeze(COH(i,j,:)),'c'); hold on;
            plot(f, squeeze(eCOH(i,j,:)),'b'); plot(f, squeeze(eCOHsth(i,j,:)),'b:');
            axis([0 fc/2 -0.05 1.05]); title(['COH' int2str(i) int2str(j)]);
        end
        
        figure(h2b); % inverse spectra and partial coherence    
        subplot(M,M,(i-1)*M+j); 
        if i==j
            plot(f, squeeze(abs(P(i,j,:))),'k');hold on; plot(f, squeeze(abs(eP(i,j,:))),'k:');
            xlim([0 fc/2]); title(['P' int2str(i)]);
        else
            plot(f, squeeze(PCOH(i,j,:)),'m');hold on;
            plot(f, squeeze(ePCOH(i,j,:)),'r'); plot(f, squeeze(ePCOHsth(i,j,:)),'r:');
            axis([0 fc/2 -0.05 1.05]); title(['PCOH' int2str(i) int2str(j)]);
        end
                
        figure(h3); % DC
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(DC(i,j,:)),'c'); hold on;
        plot(f, squeeze(eDC(i,j,:)),'b'); plot(f, squeeze(eDCsth(i,j,:)),'b:');
        axis([0 fc/2 -0.05 1.05]); title(['DC' int2str(i) int2str(j)]);

        figure(h4); % PDC
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(PDC(i,j,:)),'m'); hold on;
        plot(f, squeeze(ePDC(i,j,:)),'r'); plot(f, squeeze(ePDCsth(i,j,:)),'r:');
        axis([0 fc/2 -0.05 1.05]); title(['PDC' int2str(i) int2str(j)]);        

    end
end
figure
eAmb=eAm;
eAmb(abs(eAm)<1e-1)=0;
subplot(221);imagesc(Am);caxis([0 1]);subplot(222);imagesc(eAmb)
caxis([0 1])
subplot(223);imagesc(sum(reshape(Am,[N,N,pr]),3));caxis([0 1]);subplot(224);imagesc(sum(reshape(eAmb,[N,N,p]),3))
