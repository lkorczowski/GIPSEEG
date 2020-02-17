%% FREQUENCY DOMAIN MVAR ANALYSIS WITH INSTANTANEOUS EFFECTS
% REFERENCE: Luca Faes and Giandomenico Nollo (2011). Multivariate Frequency Domain Analysis of Causal Interactions in Physiological Time Series,
% Biomedical Engineering, Trends in Electronics, Communications and Software, Anthony N. Laskovski (Ed.), ISBN: 978-953-307-475-7, InTech,  
% Available from: http://www.intechopen.com/articles/show/title/multivariate-frequency-domain-analysis-of-causal-interactions-in-physiological-time-series 

%%% inputs: 
% Bm=[B(1)...B(p)]: M*pM matrix of the extended MVAR model (with instantaneous effects) coeffs
% B0: M*M matrix of instantaneous effects
% Sw: M*M diagonal covariance matrix of the input noises
% N= number of points for calculation of the spectral functions (nfft)
% Fs= sampling frequency

%%% outputs:
% EDC= Extended Directed Coherence (Eq. 26)
% DDC= Delayed Directed Coherence (Eq. 29)
% EPDC= Extended Partial Directed Coherence (Eq. 27)
% DPDC= Delayed Partial Directed Coherence (Eq. 29)
% COH= Coherence (Eq. 3)
% PCOH= Partial Coherence (Eq. 3)
% G= Transfer Function Matrix (after Eq. 22)
% S= Spectral Matrix (Eq. 23)
% P= Inverse Spectral Matrix (Eq. 23)
% f= frequency vector

% note: Sw MUST be given diagonal as input, since inst. eff. are in the B0 coeffs!

function [EDC,DDC,EPDC,DPDC,COH,PCOH,G,S,P,f,St,Pt] = fdMVAR0(Bm,B0,Sw,N,Fs)

M= size(Bm,1); % Am has dim M*pM
p = size(Bm,2)/M; % p is the order of the MVAR model

if nargin<3, Sw = eye(M,M); end; % if not specified, we assume uncorrelated noises with unit variance as inputs 
if nargin<4, N = 512; end;
if nargin<5, Fs= 1; end;     
if all(size(N)==1),	 %if N is scalar
    f = (0:N-1)*(Fs/(2*N)); % frequency axis
else            % if N is a vector, we assume that it is the vector of the frequencies
    f = N; N = length(N);
end;

s = exp(i*2*pi*f/Fs); % vector of complex exponentials
z = i*2*pi/Fs;


%% Initializations: spectral matrices have M rows, M columns and are calculated at each of the N frequencies
G=zeros(M,M,N); % Transfer Matrix
S=zeros(M,M,N); St=S; % Spectral Matrix
P=zeros(M,M,N); Pt=P; % Inverse spectral Matrix
COH=zeros(M,M,N); %Coherence
EDC=zeros(M,M,N); % directed coherence in presence of instantaneous effects
DDC=zeros(M,M,N); % Delayed DC Faes book chapter 2011
EPDC=zeros(M,M,N); % Extended PDC Faes  BiolCyb 2010
DPDC=zeros(M,M,N); % Delayed PDC Faes  BiolCyb 2010
tmp1=zeros(M,1); tmp8=tmp1; %denominators for DC (column!)
tmp3=tmp1'; tmp4=tmp1'; %denominators for PDC (row!)

B = [eye(M)-B0 -Bm]; % matrix from which M*M blocks are selected to calculate the spectral functions
invSw=inv(Sw); % note: Sw is always diagonal, and so also invSw

%% computation of spectral functions
for n=1:N, % at each frequency
    
        %%% Coefficient matrix in the frequency domain
        Bs = zeros(M,M); % matrix Bs(z)=I-B0-sum(B(k))
        for k = 1:p+1,
                Bs = Bs + B(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n));  %indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B (B(1) is in the second block, and so on)
        end; 
        Bt=Bs+B0; %Btilda in Eq. 28
        Gt=inv(Bt); %Gtilda in Eq. 28
        
        %%% Transfer matrix (after Eq. 22)
        G(:,:,n)  = inv(Bs);
        
        %%% Spectral matrix (Eq. 23)
        S(:,:,n)  = G(:,:,n)*Sw*G(:,:,n)';
        
        %%% Spectral matrix without instantaneous effects
        St(:,:,n)  = Gt*Sw*Gt';
        
        %%% Inverse Spectral matrix (Eq. 23)
        P(:,:,n) = inv(S(:,:,n)); % P(:,:,n) = Bs'*invSw*Bs;
        
        %%% Inverse Spectral matrix without instantaneous effects
        Pt(:,:,n) = inv(St(:,:,n)); % Pt(:,:,n) = Bt'*invSw*Bt;
        
        %%% denominators of EDC, DDC, EPDC, DPDC for each m=1,...,num. channels
        for m = 1:M,
            tmp1(m)=sqrt((abs(G(m,:,n)).^2) * diag(Sw)); % for the EDC: m-th row of G * variance of W (Sw is diag)
            tmp8(m)=sqrt((abs(Gt(m,:)).^2) * diag(Sw)); % for the DDC: m-th row of Gt * variance of W (Sw is diag)
            tmpp1 = squeeze(Bs(:,m)); % this takes the m-th column of Bs...
            tmp3(m) = sqrt(tmpp1'*invSw*tmpp1); % for EPDC
            tmpp2 = squeeze(Bt(:,m)); % this takes the m-th column of Bt...
            tmp4(m) = sqrt(tmpp2'*invSw*tmpp2); % for DPDC
        end;
        
        %%% Extended Directed Coherence (Eq. 26)
        EDC(:,:,n) = G(:,:,n)*sqrt(Sw) ./ tmp1(:,ones(M,1));
        %nota: tmp1(:,ones(M,1)) crea M colonne tutte uguali a tmp1 - la riga (ossia il den) è la stessa - trova in un colpo solo tutti i denominatori per DC
        
        %%% Delayed Directed Coherence (Eq. 29)
        DDC(:,:,n) = Gt*sqrt(Sw) ./ tmp8(:,ones(M,1));
        %nota: tmp1(:,ones(M,1)) crea M colonne tutte uguali a tmp1 - la riga (ossia il den) è la stessa - trova in un colpo solo tutti i denominatori per DC
        
        %%% Extended Partial Directed Coherence (Eq. 27)
        EPDC(:,:,n) = (sqrt(invSw)*Bs) ./ tmp3(ones(1,M),:);
        %nota: tmp3(ones(1,M),:) crea M righe tutte uguali a tmp3 - la colonna (ossia il den) è la stessa - trova in un colpo solo tutti i denominatori per PDC
        
        %%% Delayed Partial Directed Coherence (Eq. 29)
        DPDC(:,:,n) = (sqrt(invSw)*Bt) ./ tmp4(ones(1,M),:);
        
end;

%%% COHERENCE and PARTIAL COHERENCE (Eq. 3)
for m1=1:M
    for m2=1:M
        COH(m1,m2,:) = (S(m1,m2,:)) ./ sqrt(abs(S(m1,m1,:).*S(m2,m2,:)));
        PCOH(m1,m2,:) = (P(m1,m2,:)) ./ sqrt(abs(P(m1,m1,:).*P(m2,m2,:)));
    end
end

