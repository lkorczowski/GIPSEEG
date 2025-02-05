function [Xm A Sm]=simEEG(T,N,M,Fs,S,Cond,PhaseShift,SpDen,SigA,SigN)
%%
% [Xm A Sm]=simEEG(T,N,M,Fs,S,Cond,PhaseShift,SpDen,SigA,SigN)
% EEG mixture model generator:
% Xm = A*Sm + w (instantaneous linear mixture with additive gaussian noise w)
% the sources Sm are generate with pure sine signals. Each source within subjects
% has prime number frequency in order to be non-proportionnal. Each source is only
% correlated with (M-1) other sources (i.e. one source for every other subject).
% If PhaseShift=0, then this correlation is always 1.
% If PhaseShift=/=0, then each pair has an unique correlation between 1 and 0.
%
% Remark : spectral approachs will only see correlation = 1 (same freq).
% To change that, add additive gaussian noise to the sources ?
%
% INPUT(S)
%T:     length data in secondes (default: 1 s)
%N:     number of electrodes (default: 16)
%M:     number of subjects (default: 2)
%Fs:    sample rate (default: 128 Hz)
%S:     if scalar : number of sources // if matrices: sources of size
%       Nx(T*Fs) (default S=N). The sources should be sorted by subject
%       such as Sm=[S1;...;Si;...;SM] where the sources Si are
%       uncorrelated.
%Cond:  maximum conditionning of the mixture matrix (default : 1e3)
%PhaseShift:    Shift between subjects sources (0= no shift, default: pi/2)
%SpDen: density of the mixture matrices (0<SpDen<1, default: 0.3)
%SigA:  standart deviation of the coefficients of A (default: 1)
%SigN: standart deviation of the additive gaussian noise w (default: 1e-2)
%
% NOTE : if S<N, then we add N-S gaussian noise sources.
%
% OUTPUT(S)
% Xm:    simulated EEG output Nx(Fs*T)
% A:     block-diagonal sparse random mixture matrix NxN of density SpDen 
%       and std SigA
% Sm:    sources Nx(Fs*T)
%
% exemple:
% simEEG %generate 32x128 linear mixture of 32 sources correlated by pairs

% uncomment : variable for internal tests
%{
T=1
Fs=128
N=16
M=2
S=16
SpDen=0.3
Cond=1000
PhaseShift=pi/2;
SigA=1;
SigN=0.01
%}
%important args
if nargin<1 || isempty(T); T=1; end
if nargin<2 || isempty(N); N=16; end
if nargin<3 || isempty(M); M=2; end
if nargin<4 || isempty(Fs); Fs=128; end

%optinal args
if nargin<10 || isempty(SigN) ; SigN=1e-2; end
if nargin<9 || isempty(SigA) ; SigA=.3; end
if nargin<8 || isempty(SpDen) ; SpDen=.3; end
if nargin<7 || isempty(PhaseShift); PhaseShift=pi/2; end
if nargin<6 || isempty(Cond); Cond=1e3; end
if nargin<5 || isempty(S); S=N; end





check=1;
while check
    
    %generate blockdiagonal mixture matrix A
    rankA=0;
    while rankA~=N*M || Cond<C %until fullrank mixture and good conditioning
        
        
        A=[];
        for indM=1:M
            A=(blkdiag(A,(SigA*sprandn(N,N,SpDen))));% generate random mixture for each subject
        end
        
        %check for the conditioning and rank
        A=(full(A));
        rankA=rank(A);
        C=cond(A);
        
        % uncomment for plot
        %{
spy(A);text(5,-1,['cond=' num2str(C)])
EIG=abs(eig(A));
C2=max(EIG)/min(EIG)
C
        %}
    end
    
    %Generate prime frequency (for non-proportionnals sources) up to Fs/2
    Freqs=primes(Fs*10/2);
    Freqs=Freqs(randperm(length(Freqs)))/10;
    
    if size(S,1)>1 %if the sources are given, don't generate them
        Sm=S;
    else %generate the sources
        
        %generate correlated sources between the M subjects
        Sm=[];
        for indM=1:M
            phi=PhaseShift*(indM-1);%for each subject, vary the maximum shift ?
            
            %asign randomly prime frequencies to sources for all subjects
            sources=[];
            for indS=1:S
                
                time=0:1/Fs:T-1/Fs;
                % generate the sources for each subject with phase-shift
                sources=[sources ; sin(2*pi*Freqs(indS)*time +phi*(sqrt(indS/S)))];
                % REMARK : IN ORDER TO MODIFY THE SPECTRAL CORRELATION, YOU
                % COULD NEED TO CHANGE THE SHAPE OF THE SOURCE (i.e. take
                % triangle, rectangle or other function) OR ADD SOME NOISE.
            end
            
            %complete with N-S gaussian noise sources for each subject
            if S<N
                for indN=1:N-S
                    time=0:1/Fs:T-1/Fs;
                    sources=[sources ; rand(1,T*Fs)];
                end
                
            end
            
            %add the sources of each subject
            Sm=[Sm;sources];
        end
    end
    
    Xm=A*Sm; %generate the artifial mixture
    Xm=Xm+SigN*randn(size(Xm)); %add gaussian noise
    
    check=any(eig(Xm*Xm')<0); %check for full spatial rank
    
    %uncomment for visual plot
    figure(2)
subplot(211)
CovA=log(Sm*Sm'+1-min(min(Sm*Sm')));
surf(CovA) %check the correlation between the sources
hold on

figure(1)
subplot(121);plotEEG(Sm,1,128);subplot(122);plotEEG(Xm,1,128)

hold on


figure(2)
subplot(212)
surf(Xm*Xm')
    hold off
    
end


end