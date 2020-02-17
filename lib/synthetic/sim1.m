function [Xm A Sm]=sim1(T,N,M,Fs,S,Cond,PhaseShift,SpDen,SigA)
%%
% INPUT
%T:     length data in secondes (default: 1 s)
%N:     number of electrodes (default: 16)
%M:     number of subjects (default: 2)
%Fs:    sample rate (default: 128 Hz)
%S:     if scalar : number of sources // if matrices: sources of size
%       Nx(T*Fs) (default S=N)
%Cond:  maximum conditionning of the mixture matrix (default : 1e3)
%PhaseShift:    Shift between subjects sources (0= no shift, default: pi/2)
%SpDen: density of the mixture matrices (0<SpDen<1, default : 0.3)

% NOTE : if S<N, then we add N-S gaussian noise sources.

% OUTPUT
%Xm:    simulated EEG output Nx(Fs*T)
%A:     mixture matrix NxN
%Sm:    sources Nx(Fs*T)

%important args
if nargin<1 || isempty(T); T=1; end
if nargin<2 || isempty(N); N=16; end
if nargin<3 || isempty(M); M=2; end
if nargin<4 || isempty(Fs); Fs=128; end

%optinal args
if nargin<8 || isempty(SpDen) ; SpDen=.3; end
if nargin<7 || isempty(PhaseShift); PhaseShift=pi/2; end
if nargin<6 || isempty(Cond); Cond=1e3; end
if nargin<5 || isempty(S); S=N; end


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

check=1
while check
    
    %generate blockdiagonal mixture matrix A
    rankA=0;
    while rankA~=N*M || Cond<C %until fullrank mixture and good conditioning
        
        A=[];
        for indM=1:M
            A=(blkdiag(A,(SigA*sprandn(N,N,SpDen))));% generate random mixture for each subject
        end
        
        %check for the conditioning and rank
        A=full(A);
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
    
    %generate correlated sources between the M subjects
    Sm=[];
    
    for indM=1:M
        phi=PhaseShift*(indM-1);
        
        %asign randomly prime frequencies to sources for all subjects
        sources=[];
        for indS=1:S
            
            time=0:1/Fs:T-1/Fs;
            sources=[sources ; sin(2*pi*Freqs(indS)*time +phi*(sqrt(indS/S)))];
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
    check=any(eig(Xm*Xm')<0);
    %uncomment for visual plot
    %{
figure
CovA=log(Sm*Sm'+1-min(min(Sm*Sm')));
surf(CovA)

figure
plotEEG(Sm,1,128)

Xm=A*Sm;
Xm=Xm+SigN*randn(size(Xm));


figure
surf(Xm*Xm')
    %}
end


end