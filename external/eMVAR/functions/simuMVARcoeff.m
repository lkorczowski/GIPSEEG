%% THEORETICAL COEFFICIENTS FOR SIMULATED MVAR PROCESSES

function [Bm,B0,Sw]=simuMVARcoeff(caso);

error(nargchk(0,1,nargin));
if nargin < 1, caso=1; end % default simulation 1

switch caso
case {1} % theoretical example 1 book chapter
        M=4;
        pmax=2; % maximum lag
        r1=0.9; f1=0.1; % oscillation channel 1
        r2=0; f2=0.1; % oscillation channel 2
        r3=0.9; f3=0.3; % oscillation channel 3
        r4=0; f4=0.1; % oscillation channel 4
        
        %Residual covariance matrix (DIAGONAL)
        Sw(1,:)=[1 0 0 0];
        Sw(2,:)=[0 1 0 0];
        Sw(3,:)=[0 0 1 0];
        Sw(4,:)=[0 0 0 1];
        
        %Matrix of instantaneous effects
        B0(1,:)=[0 0 0 0];
        B0(2,:)=[0 0 0 0];
        B0(3,:)=[0 0 0 0];
        B0(4,:)=[0 0 0 0];
        
        Bk=NaN*ones(M,M,pmax);
        %effects at lag 1
        Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 0 0];
        Bk(2,:,1)=[1 2*r2*cos(2*pi*f2) 0.5 0];
        Bk(3,:,1)=[0 0.5 2*r3*cos(2*pi*f3) 0];
        Bk(4,:,1)=[0 0 0 2*r4*cos(2*pi*f4)];
        %effects at lag 2
        Bk(1,:,2)=[-r1^2 0 0 0];
        Bk(2,:,2)=[0 -r2^2 0 0];
        Bk(3,:,2)=[0 0.5 -r3^2 0];
        Bk(4,:,2)=[1 0 0 -r4^2];
        
case {2} % theoretical example 2 book chapter
        M=4;
        pmax=2; % massimo ritardo nella simulazione
        r1=0.95; f1=0.125;% raggio primo processo AR (f=0.25)
        r2=0.8;
        
        %cov residui (DIAGONALE)
        Sw(1,:)=[1 0 0 0];
        Sw(2,:)=[0 2 0 0];
        Sw(3,:)=[0 0 8 0];
        Sw(4,:)=[0 0 0 1];
        
        %dipendenze a lag zero (modellizzate a destra in eq (2)!)
        B0(1,:)=[0 0 0 0];
        B0(2,:)=[1 0 0 0];
        B0(3,:)=[0 0.8 0 0];
        B0(4,:)=[0 0.6 0 0];
        
        Bk=NaN*ones(M,M,pmax);
        %dipendenze da n-1
        Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 -0.4 0];
        Bk(2,:,1)=[0.2 0 0 0];
        Bk(3,:,1)=[0 0 0 0];
        Bk(4,:,1)=[0 0 0 0];
        %dipendenze da n-2
        Bk(1,:,2)=[-r1^2 0 0.3 0];
        Bk(2,:,2)=[0 -r2^2 0 0];
        Bk(3,:,2)=[0 0 0 0];
        Bk(4,:,2)=[0 0 0 0];
        %%% errata corrige: the term +0.3*y3(n-2) is missing in Eq. 30 %%%
        %%% also Fig. 5 needs to be amended, adding b13(2)=0.3, a13(2)=0.3, a33(2)=0.24
            
    case {3} % simulation 2 from Faes et al, EMBC 2010 
        M=3;
        pmax=2; % massimo ritardo nella simulazione
        r2=sqrt(0.85); f2=0.125; % processo AR per y2
        r3=sqrt(0.75); % (f3=0.25) % processo AR per y3
        
        %cov residui (DIAGONALE)
        Sw(1,:)=[1 0 0];
        Sw(2,:)=[0 2 0];
        Sw(3,:)=[0 0 3];
        
        %dipendenze a lag zero (modellizzate a destra in eq (2)!)
        B0(1,:)=[0 0 -0.6];
        B0(2,:)=[0 0 0];
        B0(3,:)=[0 0.6 0];
        
        Bk=NaN*ones(M,M,pmax);
        %dipendenze da n-1
        Bk(1,:,1)=[0 0 0];
        Bk(2,:,1)=[0 2*r2*cos(2*pi*f2) 0];
        Bk(3,:,1)=[0 0.5 0];
        %dipendenze da n-2
        Bk(1,:,2)=[0 0 0];
        Bk(2,:,2)=[0 -r2^2 0.5];
        Bk(3,:,2)=[0 0 -r3^2];
        
        
    case {4} % theoretical example tutorial paper Faes-Nollo 2011
        M=5;
        pmax=2; % maximum lag
        r1=0.9; f1=0.1; % oscillation channel 1
        r2=0; f2=0.1; % oscillation channel 2
        r3=0; f3=0.3; % oscillation channel 3
        r4=0.8; f4=0.3; % oscillation channel 4
        
        %Residual covariance matrix (DIAGONAL)
        Sw=eye(M);
        
        %Matrix of instantaneous effects
        B0=zeros(M,M);
        
        Bk=NaN*ones(M,M,pmax);
        %effects at lag 1
        Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 0 0 0];
        Bk(2,:,1)=[0.5 2*r2*cos(2*pi*f2) 0 0.5 0];
        Bk(3,:,1)=[0 0.5 2*r3*cos(2*pi*f3) 0 0];
        Bk(4,:,1)=[0 0 0.5 2*r4*cos(2*pi*f4) 0];
        Bk(5,:,1)=[0.5 0 0 0 0];
        %effects at lag 2
        Bk(1,:,2)=[-r1^2 0 0 0 0];
        Bk(2,:,2)=[0 -r2^2 0 0 0];
        Bk(3,:,2)=[0 0.5 -r3^2 0 0];
        Bk(4,:,2)=[0 0 0.5 -r4^2 0];
        Bk(5,:,2)=[0.5 0 0 0 0];
        

    case {5} % Example 2 Faes and Nollo Biological Cybernetics 2010
        M=3;
        pmax=2; % massimo ritardo nella simulazione      
        r1=0.9; f1=0.1;
        r2=0.95; 
        %Residual covariance matrix (DIAGONAL)
        Sw=eye(M);
        %effects at lag 0
        B0(1,:)=[0 0 0];
        B0(2,:)=[1 0 0];
        B0(3,:)=[0 0.8 0];
        Bk=NaN*ones(M,M,pmax);
        %effects at lag 1
        Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 0];
        Bk(2,:,1)=[0 0 0];
        Bk(3,:,1)=[0 0.2 0];
        %effects at lag 2
        Bk(1,:,2)=[-r1^2 0 0];
        Bk(2,:,2)=[0 -r2^2 0];
        Bk(3,:,2)=[0 0 0];

    case {6} % Example 3 Faes and Nollo Biological Cybernetics 2010 
        M=4;
        pmax=2; % massimo ritardo nella simulazione
        r1=0.95; f1=0.125;% raggio primo processo AR (f=0.25)
        r2=0.8; f2=0.25;
        %Residual covariance matrix (DIAGONAL)
        Sw(1,:)=[1 0 0 0];
        Sw(2,:)=[0 2 0 0];
        Sw(3,:)=[0 0 8 0];
        Sw(4,:)=[0 0 0 1];
        %effects at lag 0
        B0(1,:)=[0 0 0 0];
        B0(2,:)=[1 0 0 0];
        B0(3,:)=[0 0.8 0 0];
        B0(4,:)=[0 0.6 0 0];
        Bk=NaN*ones(M,M,pmax);
        %effects at lag 1
        Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 -0.4 0];
        Bk(2,:,1)=[0.2 0 0 0];
        Bk(3,:,1)=[0 0 0 0];
        Bk(4,:,1)=[0 0 0 0];
        %effects at lag 2
        Bk(1,:,2)=[-r1^2 0 0.3 0];
        Bk(2,:,2)=[0 -r2^2 0 0];
        Bk(3,:,2)=[0 0 0 0];
        Bk(4,:,2)=[0 0 0 0];
               
        
%%
    otherwise
        pmax=0;
end

%% concateno in matrice Bm
Bm=[];
for kk=1:pmax
    Bm=[Bm Bk(:,:,kk)];
end
