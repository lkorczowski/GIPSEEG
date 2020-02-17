%% THEORETICAL EXAMPLE 1 (Sect. 3.3 Book Chapter)
% REFERENCE: Luca Faes and Giandomenico Nollo (2011). Multivariate Frequency Domain Analysis of Causal Interactions in Physiological Time Series,
% Biomedical Engineering, Trends in Electronics, Communications and Software, Anthony N. Laskovski (Ed.), ISBN: 978-953-307-475-7, InTech,  
% Available from: http://www.intechopen.com/articles/show/title/multivariate-frequency-domain-analysis-of-causal-interactions-in-physiological-time-series 

clear; close all; clc;
nfft=512; % number of frequency bins
fc=1; % sample frequency
simunumber=1;

%% MVAR coefficients
% coefficients of extended MVAR model
[Bm,B0,Sw]=simuMVARcoeff(simunumber);

% in this case B0=0, the model is strictly causal
Am=Bm;
Su=Sw;

M=size(Bm,1);
p=size(Bm,2)/M;

%%%%%%%%%%%%%%%%%%% tmp: show one realization of the process
% N=500;
% U=randn(M,N); % uncorrelated gaussian innovations
% [Y]=MVARfilter(Am,U); % realization of Eq. 18
% figure(1); % realizations of functions
% for i=1:M
%     subplot(M,1,i); plot(Y(i,:));
%     xlim([0 N]);
%     xlabel(['n']);ylabel(['y_' int2str(i)]);
% end
%%%%%%%%%%%%%%%%%%%

%% Theoretical spectral functions
[dc,dtf,pdc,gpdc,coh,pcoh,pcoh2,h,s,pp,f] = fdMVAR(Am,Su,nfft,fc);
S=abs(s); % spectral matrix
P=abs(pp); % inverse spectral matrix
DC=abs(dc).^2; % directed coherence
PDC=abs(gpdc).^2; % partial directed coherence
COH=abs(coh).^2; %coherence
PCOH=abs(pcoh).^2; % partial coherence
PCOH2=abs(pcoh2).^2; % partial coherence, other definition to check equality

%% Partial spectra
Sd=NaN*ones(M,M,nfft); % Decomposition of spectru
Pd=Sd; % Decomposition of inverse spectrum
for i=1:M
    for m=1:M
        Sd(m,i,:)=DC(i,m,:) .* S(i,i,:); % Eq. 13
        Pd(m,i,:)=PDC(m,i,:) .* P(i,i,:); % Eq. 17
    end
end


%% graphs & disp
h2a=figure('numbertitle','off','name','Fig. 2a'); % spectra and coherence
h2b=figure('numbertitle','off','name','Fig. 2b'); % inverse spectra and partial coherence
h3=figure('numbertitle','off','name','Fig. 3'); % Spectrum decomposition and directed coherence
h4=figure('numbertitle','off','name','Fig. 4'); % Inverse spectrum decomposition and partial directed coherence
% mappa=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 1 1 1];
mappa=[1 0 0; 1 1 0; 0 1 0; 0 1 1];
tmp=[];
for i=1:M    
    for m=1:M, Sim(:,m)=squeeze(Sd(m,i,:)); end 
    for m=1:M, Pim(:,m)=squeeze(Pd(m,i,:)); end 
    
    figure(h3); % Decomposition of spectrum
    subplot(M,M+1,(i-1)*M+i); area(f',Sim); title(['S' int2str(i) int2str(i)]); colormap(mappa);
    
    figure(h4); % Decomposition of inverse spectrum
    subplot(M,M+1,(i-1)*M+i); area(f',Pim); title(['P' int2str(i) int2str(i)]); colormap(mappa);

    
    for j=1:M
        figure(h2a); % spectra and coherence       
        subplot(M,M,(i-1)*M+j); 
        if i==j
            plot(f, squeeze(abs(S(i,j,:)))/max(squeeze(abs(S(i,j,:)))),'k');
            axis([0 fc/2 0 1]); title(['S' int2str(i)]);
        else
            plot(f, squeeze(COH(i,j,:)),'b');
            axis([0 fc/2 -0.05 1.05]); title(['COH' int2str(i) int2str(j)]); %line([0.1 0.1], [0 1], 'Color', [.8 .8 .8]);
        end
        
        figure(h2b); % inverse spectra and partial coherence    
        subplot(M,M,(i-1)*M+j); 
        if i==j
            plot(f, squeeze(abs(P(i,j,:)))/max(squeeze(abs(P(i,j,:)))),'k');
            axis([0 fc/2 -0.05 1.05]); title(['P' int2str(i)]);
        else
            plot(f, squeeze(PCOH(i,j,:)),'r');hold on; plot(f, squeeze(PCOH2(i,j,:)),'r--');
            axis([0 fc/2 -0.05 1.05]); title(['PCOH' int2str(i) int2str(j)]); %line([0.1 0.1], [0 1], 'Color', [.8 .8 .8]);
        end
                
        figure(h3); % DC
        subplot(M,M+1,(i-1)*M+j+i);
        plot(f, squeeze(DC(i,j,:)),'b'); 
        axis([0 fc/2 -0.05 1.05]); title(['DC' int2str(i) int2str(j)]); % line([0.1 0.1], [0 1], 'Color', [.8 .8 .8]);

        figure(h4); % PDC
        subplot(M,M+1,(i-1)*M+j+i);
        plot(f, squeeze(PDC(i,j,:)),'r'); 
        axis([0 fc/2 -0.05 1.05]); title(['PDC' int2str(i) int2str(j)]);        

    end
end

