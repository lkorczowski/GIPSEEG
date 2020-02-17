function [Bs Bt As At Class eigV]=CSTPover(X,Y,P,W,winElec,winTime,Pzf,Xnoise,Ynoise,PLOT)
% OLD VERSION OF THE CSTP, TO REMOVE
% [Bs Bt As At Class eigV]=CSTPover(X,Y,P,W,winElec,winTime,Pzf,PLOT)
% Compute Common Spatio-Temporal Pattern for 2 classes to maximize the SNR
% with Xhat=As*Bs'*X*Bt*At'
% Xhat is the estimated signal
% Bs is the spatial filter (features)
% Bt is the temporal filter (features)
% As and At allows to go back in the sensors spaces
%
%  *** History: 2015-03-19
%  *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
%  *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
if nargin<10
    PLOT=0;
end
if nargin<8
    Xnoise=X;
    Ynoise=Y;
end
if nargin<5 || isempty(winElec) || isempty(winTime)
    winElec=1:size(X,1);
    winTime=1:size(X,2);
end


if nargin<4 || isempty(W); 
    W=ones(1,size(X,3));
end


if nargin<3 || isempty(P) %check and verify the average P300
    P=[];Bool=1; %non overlapping method
else
    Bool=0; %use already computed mean P300
    if ~iscell(P) %check if the ERP is with a cell structure
        tmp=P;
        nbClass=size(P,3);
        clear P
       for z=1:nbClass
           P{z}=tmp(:,:,z);%one cell for each class
       end
    end
end
% Algo Parameters
eigWhite1=1:size(X,1);
eigWhite2=1:size(X,2);
eigV=1:12;%reduction
if nargin<2 || isempty(Y)
    Y=ones(size(X,3),1);
end


Class=unique(Y);
nbClass=length(Class);
Nbe=size(X,1);

%%%%%%%%%%%%%%%% START COMPUTING CSPT %%%%%%%%%%%%%%%%%%%%%%%%


%compute Sample covariance matrices
for k=1:size(Xnoise,3)
    Cs_k(:,:,k)=Xnoise(:,:,k)*Xnoise(:,:,k)'/(size(Xnoise,2)-1);
    Ct_k(:,:,k)=Xnoise(:,:,k)'*Xnoise(:,:,k)/(size(Xnoise,1)-1);
end

%mean covariance matrices
Cs=[];
Ct=[];


for z=1:nbClass%find(Class==1) %%%%%%%%%%%%%%%%%%%% WARNING
    if Bool==1
            P{z}=EnsembleAverage(X(:,:,Y==Class(z)));
    end
Cs=cat(3,Cs,mean(Cs_k(:,:,Ynoise==Class(z)),3)); %spatial signal + noise covariance matrix for each class
Ct=cat(3,Ct,mean(Ct_k(:,:,Ynoise==Class(z)),3)); %temporal signal + noise covariance matrix for each class
end

 % mean covariance matrix
% if nbClass>1
Cs=mean(Cs,3);
%Cs=mean(Cs_k,3);% Noise covariance matrix
%Cs=eye(size(Cs));

Ct=mean(Ct,3);
%Ct=mean(Ct_k,3);% Noise covariance matrix
%Ct=eye(size(Ct));
% surf(mean(Cs_k,3))
%  surf(Cs(6:end,6:end,1))

 %compute expected SNR ratio
%{ 
SNRden=zeros(1,nbClass);
for z=1:find(1==Class)%nbClass
    SNRnum=norm(P{z},'fro').^2;
    nbk=length(find(Y==Class(z)));
    Xkz=X(:,:,Y==Class(z));
    for k=1:nbk
        SNRden(z)=SNRden(z)+norm(squeeze(Xkz(:,:,k)),'fro').^2;
    end
    SNRden(z)=SNRden(z)/nbk;
end
SNR=SNRnum/SNRden(Class==1);
%}
[Ux Ws]=eig(Cs);
[Vx Wt]=eig(Ct);
[Ws,ind] = sort(diag(Ws),'descend');

sumWs=cumsum(Ws);
indWhite1tmp=find(sumWs>max(sumWs)*(1-1e-12),1);
if max(eigWhite1)>indWhite1tmp
    eigWhite1=1:indWhite1tmp;
end
if any(Ws(eigWhite1)<=0); disp('Warning negative SPATIAL eigenvalue(s)'); end

Ux = Ux(:,ind);
Ux = Ux(:,eigWhite1);
[Wt,ind] = sort(diag(Wt),'descend');
%remove too small eigenvalues (10e-9)

sumWt=cumsum(Wt);
indWhite2tmp=find(sumWt>max(sumWt)*(1-1e-12),1);
if max(eigWhite2)>indWhite2tmp
    eigWhite2=1:indWhite2tmp;
end
if any(Wt(eigWhite2)<=0); disp('Warning negative TEMPORAL eigenvalue(s)'); end

Vx = Vx(:,ind);
Vx = Vx(:,eigWhite2);
Fs = (sqrt(pinv(diag(Ws(eigWhite1)))) * Ux')'; 
Ft =(sqrt(pinv(diag(Wt(eigWhite2)))) * Vx')';
Gs =Ux*sqrt((diag(Ws(eigWhite1)))) ;
Gt =Vx*sqrt((diag(Wt(eigWhite2)))) ;

%compute CSPT
for z=1:length(P)
    Z{z}=Fs'*P{z}*Ft;
    [Uz{z},Wz{z},Vz{z}]=svd((Z{z}));
    %Uz{z}*Wz{z}*Vz{z}'
    %%############# FIND OPTIMAL Pz ##################
    %eigV=find(diag(Wz{z}.^2)>SNR); %%old version
    if exist('Pzf')
        eigV=1:Pzf;
    else
        Fdenum=norm(Z{z},'fro').^2;
        for Pz=Nbe*0.75:Nbe-1

            Uz_r{z}=Uz{z}(:,1:Pz);
            Vz_r{z}=Vz{z}(:,1:Pz);
            Zz=Uz_r{z}*Wz{z}(Pz,Pz)*Vz_r{z}';
            Fnum=norm(Zz(winElec,winTime),'fro').^2;
            indPz(Pz)=Fnum/Fdenum;

        end
        [tmp Pze]=max(indPz);
        eigV=1:Pze;
    end
    %%############# FIND OPTIMAL Pz END ##################

 
    if isempty(eigV)
        eigV=1;
    end
    
    Uz_r{z}=Uz{z}(:,eigV);
    Vz_r{z}=Vz{z}(:,eigV);
    

    
    Bs{z}=(Uz_r{z}'*Fs')';
    Bt{z}=(Vz_r{z}'*Ft')';
    As{z}=Gs*Uz{z}(:,eigV);
    At{z}=Gt*Vz{z}(:,eigV);

    
end
