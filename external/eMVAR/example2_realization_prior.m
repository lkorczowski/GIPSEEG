%% Identification of extended MVAR model and estimation of frequency domain causality
% REFERENCE: Luca Faes and Giandomenico Nollo (2011). Multivariate Frequency Domain Analysis of Causal Interactions in Physiological Time Series,
% Biomedical Engineering, Trends in Electronics, Communications and Software, Anthony N. Laskovski (Ed.), ISBN: 978-953-307-475-7, InTech,  
% Available from: http://www.intechopen.com/articles/show/title/multivariate-frequency-domain-analysis-of-causal-interactions-in-physiological-time-series 

% Here we generate a realization of Eq. 30, then we identify model coefficients
% and estimate frequency domain causality, comparing to theoretical functions
% this script tests the identification method making use of prior information (approach 1)

clear; close all; clc;
nfft=512; % number of frequency bins
fc=1; % sample frequency
N=500; % length of simulated series
pcrit='aic'; % 'aic', 'mdl', or number for fixed model order
idMode=0; %identification algorithm (0:Least Squares covariance; 1:Yule-Walker correlation; 2:Vieira-Morf partial correlation; 3:Nutall-Strandpartial correlation; other: see mvar.m)

simunumber=2;
ki=[1 2 3 4]'; %works only if the order of inst effects is known (e.g., [1 2 3 4] for simu 2)!
% ki=[2 3 4 1]'; %works only if the order of inst effects is known (e.g., [1 2 3 4] for simu 2)!

%%% if simunumber=3 try also these, otherwise keep commented
% simunumber=3;
% % ki=[2 3 1]'; %correct
% ki=[1 2 3]';% wrong

do_surrogates='n'; %'y' performs surrogates, otherwise do not perform
numsurr=10; %number of surrogate time series
upalpha=95; %upper percentile of threshold

%% Theoretical MVAR coefficients and spectral functions
% coefficients of extended MVAR model
[Bm,B0,Sw]=simuMVARcoeff(simunumber);
% corresponding coefficients of strictly causal MVAR model
[Am,Su]=diag_coeff_rev(Bm,B0,Sw);

M=size(Bm,1);
% p=size(Bm,2)/M;

%%% Theoretical spectral functions
[edc,ddc,epdc,dpdc,coh,pcoh,g,s,pp,f] = fdMVAR0(Bm,B0,Sw,nfft,fc);
S=abs(s); % spectral matrix
P=abs(pp); % inverse spectral matrix
LDC=abs(ddc).^2; % lagged directed coherence
LPDC=abs(dpdc).^2; % lagged partial directed coherence
EDC=abs(edc).^2; % extended directed coherence
EPDC=abs(epdc).^2; % extended partial directed coherence
COH=abs(coh).^2; %coherence
PCOH=abs(pcoh).^2; % partial coherence

%% Estimated MVAR coefficients and spectral functions 
U=InstModelfilter(N,Sw,'ExtendedGauss',B0); % gaussian innovations W, U determined from B0
[Y]=MVARfilter(Am,U); % realization of Eq. 30

%model order selection
if pcrit(1)=='a' | pcrit(1)=='m' 
    [pottaic,pottmdl,aic,mdl] = mos_idMVAR(Y,20,idMode);
    if pcrit(1)=='a', p=pottaic; else p=pottmdl; end
else
    p=pcrit;
end

% model identification
[eBm,eB0,eSw,eAm,eSu,eU,eW]=idMVAR0prior(Y,p,ki,idMode);


%%% Estimated spectral functions
[edc2,ddc2,epdc2,dpdc2,coh2,pcoh2,g2,s2,pp2] = fdMVAR0(eBm,eB0,eSw,nfft,fc);
eS=abs(s2); % spectral matrix
eP=abs(pp2); % inverse spectral matrix
eLDC=abs(ddc2).^2; % lagged directed coherence
eLPDC=abs(dpdc2).^2; % lagged partial directed coherence
eEDC=abs(edc2).^2; % extended directed coherence
eEPDC=abs(epdc2).^2; % extended partial directed coherence
eCOH=abs(coh2).^2; %coherence
ePCOH=abs(pcoh2).^2; % partial coherence


%% testing model assumptions
alpha=0.05; %significance for the tests

%%% WHITENESS: portmanteau Ljung-Box test
h=50; %number of lags
[pval,Qh,critlo,crithi,crit1tail,C0,stringWhite,flagWhite]=test_whiteness(eU,p,h,alpha);

%%% INDEPENDENCE among U residuals
testtype='Kendall';
[pMU,rhoU,contarigettU,stringIndU]=test_independence(eU,testtype,alpha);
[pMW,rhoW,contarigettW,stringIndW]=test_independence(eW,testtype,alpha);

%%% NONGAUSSIANITY of W residuals
[pGauss,ps,pk,lambdas,lambdak,crittresh,stringGauss,flagGauss]=test_gaussianity(eW,alpha);



%% surrogate time series - significance thresholds
if do_surrogates=='y'
    eCOHs=NaN*ones(M,M,nfft,numsurr); ePCOHs=eCOHs; eEDCs=eCOHs; eEPDCs=eCOHs; eLDCs=eCOHs; eLPDCs=eCOHs;
    for is=1:numsurr 
        disp(['surrogate ' int2str(is) ' of ' int2str(numsurr)])
        Ys=surrVFT(Y); % FT surrogates
        [eBms,eB0s,eSws,eAms,eSus,eUs,eWs]=idMVAR0prior(Ys,p,ki,idMode);
        [edcs,ddcs,epdcs,dpdcs,cohs,pcohs] = fdMVAR0(eBms,eB0s,eSws,nfft,fc);
        eCOHs(:,:,:,is)=abs(cohs).^2; %surro coherence
        ePCOHs(:,:,:,is)=abs(pcohs).^2; %surro partial coherence (one should better devise ad-hoc surrogates)
        
        % surrogates for DC and PDC
        for ii=1:M
            for jj=1:M
                if ii~=jj
                    Ys=surrVCFTf0(Y,eBm,eB0,eSw,ii,jj); %CFTf0 surrogates
                    [eBms,eB0s,eSws,eAms,eSus,eUs,eWs]=idMVAR0prior(Ys,p,ki,idMode);
                    [edcs,ddcs] = fdMVAR0(eBms,eB0s,eSws,nfft,fc);
                    eEDCs(ii,jj,:,is)=abs(edcs(ii,jj,:)).^2; %surro EDC
                    eLDCs(ii,jj,:,is)=abs(ddcs(ii,jj,:)).^2; %surro LDC
                    
                    Ys=surrVCFTd0(Y,eBm,eB0,eSw,ii,jj); %CFTd0 surrogates
                    [eBms,eB0s,eSws,eAms,eSus,eUs,eWs]=idMVAR0prior(Ys,p,ki,idMode);
                    [edcs,ddcs,epdcs,dpdcs] = fdMVAR0(eBms,eB0s,eSws,nfft,fc);
                    eEPDCs(ii,jj,:,is)=abs(epdcs(ii,jj,:)).^2; %surro EPDC
                    eLPDCs(ii,jj,:,is)=abs(dpdcs(ii,jj,:)).^2; %surro LPDC      
                end
            end
        end
        
    end
    eCOHsth=prctile(eCOHs,upalpha,4); ePCOHsth=prctile(ePCOHs,upalpha,4);; %thresholds
    eEDCsth=prctile(eEDCs,upalpha,4); eLDCsth=prctile(eLDCs,upalpha,4);
    eEPDCsth=prctile(eEPDCs,upalpha,4); eLPDCsth=prctile(eLPDCs,upalpha,4);
 
else
    eCOHsth=NaN*ones(M,M,nfft); ePCOHsth=eCOHsth;
    eEDCsth=eCOHsth; eEPDCsth=eCOHsth;
    eLDCsth=eCOHsth; eLPDCsth=eCOHsth;
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
disp([int2str(contarigettU) ' pairs of U residuals are dependent']);
disp(stringIndU);
disp([int2str(contarigettW) ' pairs of W residuals are dependent']);
disp(stringIndW);
disp([' ']);
disp('Gaussianity (Jarque-Bera) test:');
disp(['p=' num2str(pGauss) ', ' stringGauss]);




%% graphs
figure('numbertitle','off','name','Simulated time series'); % realizations of functions
for i=1:M
    subplot(M,1,i); plot(Y(i,:));
    xlim([0 N]);
    xlabel(['n']);ylabel(['y_' int2str(i)]);
end

h2a=figure('numbertitle','off','name','Spectra and Coherence (cyan:original; blue:estimated; blue-dashed:threshold)');
h2b=figure('numbertitle','off','name','Inverse spectra and Partial Coherence magenta:original; red:estimated; red-dashed:threshold)');
h3=figure('numbertitle','off','name','Extended Directed Coherence (cyan:original; blue:estimated; blue-dashed:threshold)');
h4=figure('numbertitle','off','name','Extended Partial Directed Coherence magenta:original; red:estimated; red-dashed:threshold)');
h5=figure('numbertitle','off','name','Lagged Directed Coherence (cyan:original; blue:estimated; blue-dashed:threshold)');
h6=figure('numbertitle','off','name','Lagged Partial Directed Coherence magenta:original; red:estimated; red-dashed:threshold)');
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
                
        figure(h3); % eDC
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(EDC(i,j,:)),'c'); hold on;
        plot(f, squeeze(eEDC(i,j,:)),'b'); plot(f, squeeze(eEDCsth(i,j,:)),'b:');
        axis([0 fc/2 -0.05 1.05]); title(['eDC' int2str(i) int2str(j)]);

        figure(h4); % ePDC
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(EPDC(i,j,:)),'m'); hold on;
        plot(f, squeeze(eEPDC(i,j,:)),'r'); plot(f, squeeze(eEPDCsth(i,j,:)),'r:');
        axis([0 fc/2 -0.05 1.05]); title(['ePDC' int2str(i) int2str(j)]);
        
        figure(h5); % lDC
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(LDC(i,j,:)),'c'); hold on;
        plot(f, squeeze(eLDC(i,j,:)),'b'); plot(f, squeeze(eLDCsth(i,j,:)),'b:');
        axis([0 fc/2 -0.05 1.05]); title(['lDC' int2str(i) int2str(j)]);
        
        figure(h6); % lPDC
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(LPDC(i,j,:)),'m'); hold on;
        plot(f, squeeze(eLPDC(i,j,:)),'r'); plot(f, squeeze(eLPDCsth(i,j,:)),'r:');
        axis([0 fc/2 -0.05 1.05]); title(['lPDC' int2str(i) int2str(j)]);

    end
end
