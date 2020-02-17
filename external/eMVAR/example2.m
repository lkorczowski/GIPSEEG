%% THEORETICAL EXAMPLE 2 (Sect. 4.3)
clear; close all; clc;
nfft=512; % number of frequency bins
fc=1; % sample frequency
simunumber=2;

%% expected MVAR models
% coefficients of extended MVAR model
[Bm,B0,Sw]=simuMVARcoeff(simunumber);

M=size(Bm,1);
p=size(Bm,2)/M;%lag massimo presente nel mod. MVAR = ordine del modello mvar

% coefficients of strictly causal MVAR model
[Am, Su]=diag_coeff_rev(Bm,B0,Sw);

%% Theoretical spectral functions
%%% Strcitly causal MVAR model
[dc, dtf, pdc,gpdc,coh,pcoh,pcoh2,h,s,pp,f] = fdMVAR(Am,Su,nfft,fc);
S=abs(s);
PP=abs(pp); %inverse spectral matrix
DC=abs(dc).^2;
PDC=abs(gpdc).^2;
COH=abs(coh).^2;
PCOH=abs(pcoh).^2;
PCOH2=abs(pcoh2).^2;

%%% Extended MVAR model
[edc,ddc,epdc,dpdc,coh0,pcoh0,g,s0,pp0,f] = fdMVAR0(Bm,B0,Sw,nfft,fc);
S0=abs(s0);
PP0=abs(pp0);
EDC=abs(edc).^2;
DC0=abs(ddc).^2;
EPDC=abs(epdc).^2;
PDC0=abs(dpdc).^2;
COH0=abs(coh0).^2;
PCOH0=abs(pcoh0).^2;


%% graphs & disp
h6a=figure('numbertitle','off','name','Fig. 6a'); % spectra and coherence
h6b=figure('numbertitle','off','name','Fig. 6b'); % inverse spectra and partial coherence
h7a=figure('numbertitle','off','name','Fig. 7a'); % PDC
h7b=figure('numbertitle','off','name','Fig. 7b'); % ePDC
h8a=figure('numbertitle','off','name','Fig. 8a'); % DC
h8b=figure('numbertitle','off','name','Fig. 8b'); % eDC
for i=1:M
    for j=1:M              
        figure(h6a); % spectra and coherence
        subplot(M,M,(i-1)*M+j); 
        if i==j
            plot(f, squeeze(abs(S(i,j,:))),'k');  hold on; plot(f, squeeze(abs(S0(i,j,:))),'k--');
            xlim([0 fc/2]); zoom yon;
        else
            plot(f, squeeze(COH(i,j,:)),'b');hold on;
            plot(f, squeeze(COH0(i,j,:)),'b--');
            axis([0 fc/2 -0.05 1.05]);
        end
        
        figure(h6b); % inverse spectra and partial coherence    
        subplot(M,M,(i-1)*M+j); 
        if i==j
            plot(f, squeeze(abs(PP(i,j,:))),'k');  hold on; plot(f, squeeze(abs(PP0(i,j,:))),'k--');
            xlim([0 fc/2]); zoom yon;
        else
            plot(f, squeeze(PCOH(i,j,:)),'r');hold on;plot(f, squeeze(PCOH2(i,j,:)),'r--');
            plot(f, squeeze(PCOH0(i,j,:)),'r:');
            axis([0 fc/2 -0.05 1.05]);
        end

        figure(h7a); % PDC
        subplot(M,M,(i-1)*M+j); plot(f, squeeze(PDC(i,j,:)),'r');
        hold on; plot(f, squeeze(PDC0(i,j,:)),'k--');
        hold off; axis([0 fc/2 -0.05 1.05]);
        
        figure(h7b); % ePDC
        subplot(M,M,(i-1)*M+j); plot(f, squeeze(EPDC(i,j,:)),'k');
        hold off; axis([0 fc/2 -0.05 1.05]);

        figure(h8a); % DC
        subplot(M,M,(i-1)*M+j); plot(f, squeeze(DC(i,j,:)),'b');
        hold on; plot(f, squeeze(DC0(i,j,:)),'k--');
        hold off; axis([0 fc/2 -0.05 1.05]);
        
        figure(h8b); % eDC
        subplot(M,M,(i-1)*M+j); plot(f, squeeze(EDC(i,j,:)),'k');
        hold off; axis([0 fc/2 -0.05 1.05]);
                 

    end
end
