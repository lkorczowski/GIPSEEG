function [Bs Bt As At Class eigV fullBs fullBt fullAs fullAt]=CSTP(X,Y,P,W,winElec,winTime,Pzf,Xnoise,Ynoise)
% [Bs Bt As At Class eigV]=ACSTP(X,Y,P,W,winElec,winTime,Pzf,Xnoise,Ynoise)
%
% ************************************************************************
% Compute the Common Spatio-Temporal Pattern for Z classes of ERP
% to maximize the SNR of each class such as :
%
% Xhat(:,:,k)=As*Bs'*W(k)*X(:,:,k)*Bt*At'; (0.5) & (0.8)
%
% With X(:,:,k) the raw single sweep ERP and
% Xhat(:,:,k) the estimation of a single sweep ERP that maximize the SNR.
%
% INPUTS :
% ------
% X: is an EEG signal [ nb electrodes x nb samples x nb epochs] with
%   corrected latency (or not), see LatencyEstimation.
% Y: is the class' numerical tags of the epochs of X [nb epochs]. By
%   default, all the epochs of X are in the same class with 1 as tag.
% W: are the Weights estimation [nb epochs]. Default is 1 for each epoch.
%   see WeightsEstimation.
% P: is a matrix containing the ensemble average for each class,
%   with size [nb electrodes x nb samples x nb classes]. The classes should
%   be sorted with respect to the numeral tag of Y (ascending order). see
%   meanOverlap.
% winElec, winTime: are the spatial and temporals masks for the selection of
%   the optimal subset (0.26). By default P will be the arithmetic ensemble
%   average of X for each class in Y. Useless if Pzf is given.
% Pzf: it's the desire Pz for the subspace reduction. By default, Pz is
%   estimated (poorly) such as (0.27) with masks above. For maximum
%   accuracy, we advice to compute ACSTP with several Pzf and choose the
%   best such as (0.26) is maximized. See Best_Pz.
% Xnoise, Ynoise: are the EEG epochs and tags in case that the user want to
%   use different epochs for estimation of the noise covariances matrixes in
%   (0.15) such as resting state or else. By default we use X and Y. Ynoise
%   is currently useless.
%
% OUTPUTS :
% ------
% Bs, Bt: are respectively the spatial and temporal filters that project X(k)
%   in the feature space, see (0.5)
% As, At: are respectively the spatial and temporal filters that project back
%   the features into the sensor space, see (0.8)
% Class: this output give the order of each class in the filters As At As
%   At
% eigV: gives the eigenvectors numbers in case of automatic subspace
%   dimension selection (Pzf is NOT given in input). With Pzf, then eigV is
%   (1... Pzf).
% fullBs, fullBt, fullAs, fullAt: are respectively the spatial and temporal filters without
%   the dimension reduction (for display purpose such as components
%   topoplot)
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also: ACSTP, applyCSTP, CSTPinfo

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

if nargin<2 || isempty(Y)
    Y=ones(size(X,3),1);
end



Class=unique(Y);
nbClass=length(Class);
Nbe=size(X,1);

if Bool==1
    for z=1:nbClass
        %estimate the EAE if needed
        P{z}=EnsembleAverage(X(:,:,Y==Class(z))); %(0.10)
    end
end

%%%%%%%%%%%%%%%% START COMPUTING CSPT %%%%%%%%%%%%%%%%%%%%%%%%
%compute Sample covariance matrices for each sweep
for k=1:size(Xnoise,3)
    Cs_k(:,:,k)=Xnoise(:,:,k)*Xnoise(:,:,k)'/(size(Xnoise,2)-1);
    Ct_k(:,:,k)=Xnoise(:,:,k)'*Xnoise(:,:,k)/(size(Xnoise,1)-1);
end

%mean covariance matrices
Cs=mean(Cs_k,3); %(0.15)
Ct=mean(Ct_k,3); %(0.15)

%%% noise covariance matrix, version before february 2015
%{
%compute Sample covariance matrices for each sweep
for k=1:size(Xnoise,3)
    Cs_k(:,:,k)=Xnoise(:,:,k)*Xnoise(:,:,k)'/(size(Xnoise,2)-1);
    Ct_k(:,:,k)=Xnoise(:,:,k)'*Xnoise(:,:,k)/(size(Xnoise,1)-1);
end

%mean covariance matrices
Cs=[];
Ct=[];


for z=1:nbClass
    if Bool==1
        %estimate the EAE if needed
        P{z}=EnsembleAverage(X(:,:,Y==Class(z))); %(0.10)
    end
    Cs=cat(3,Cs,mean(Cs_k(:,:,Ynoise==Class(z)),3));
    %mean spatial noise covariance matrix for each class
    
    Ct=cat(3,Ct,mean(Ct_k(:,:,Ynoise==Class(z)),3));
    %mean temporal noise covariance matrix for each class
end

% average the noise covariance matrix of the classes
Cs=mean(Cs,3); %(0.15)
Ct=mean(Ct,3); %(0.15)
%}

%%%%%%%%%%%%%%%%%%%%% BILINEAR WHITENING %%%%%%%%%%%%%%%%%%%%%%%%
[Ux Ws]=eig(Cs); %(0.18)
[Vx Wt]=eig(Ct); %(0.18)
[Ws,ind] = sort(diag(Ws),'descend');
sumWs=cumsum(Ws);
indWhite1tmp=find(sumWs>max(sumWs)*(1-1e-12),1);
% suppress only the smallest spatial eigenvectors to keep (1-1e-12)% of the
% power of the signal

if max(eigWhite1)>indWhite1tmp
    eigWhite1=1:indWhite1tmp;
end
if any(Ws(eigWhite1)<=0); disp('Warning negative SPATIAL eigenvalue(s)'); end

Ux = Ux(:,ind);%sort spatial eigenvector
Ux = Ux(:,eigWhite1);%whitening reduction
[Wt,ind] = sort(diag(Wt),'descend');
sumWt=cumsum(Wt);
indWhite2tmp=find(sumWt>max(sumWt)*(1-1e-12),1);
% suppress only the smallest temporal eigenvectors to keep (1-1e-12)% of the
% power of the signal

if max(eigWhite2)>indWhite2tmp
    eigWhite2=1:indWhite2tmp;
end
if any(Wt(eigWhite2)<=0); disp('Warning negative TEMPORAL eigenvalue(s)'); end

Vx = Vx(:,ind); %sort temporal eigenvector
Vx = Vx(:,eigWhite2);%whitening reduction

Fs = (sqrt(pinv(diag(Ws(eigWhite1)))) * Ux')'; %(0.19)
Ft =(sqrt(pinv(diag(Wt(eigWhite2)))) * Vx')'; %(0.19)
Gs =Ux*sqrt((diag(Ws(eigWhite1)))) ; %(0.20)
Gt =Vx*sqrt((diag(Wt(eigWhite2)))) ; %(0.20)
%%%%%%%%%%%%%% WHITENING MATRICES COMPUTED %%%%%%%%%%%%%%%


%% %%%%%%%%%%%% COMPUTE CSTP FILTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for z=1:length(P) %for each class
    Z{z}=Fs'*P{z}*Ft; %EAE whitening (0.22)
    [Uz{z},Wz{z},Vz{z}]=svd((Z{z})); % (0.23)

    %[Uz{z},Wz{z},Vz{z}]=svds((Z{z}),min(indWhite1tmp,indWhite2tmp)); % (0.23)
    
    %figure;plot(diag(Wz{z}))
    
    if exist('Pzf') %if Pzf has be given by the user (adviced)
        eigV=1:Pzf;
    else
        %%%%%%%%%%%%%% FIND OPTIMAL Pz (optinal) %%%%%%%%%%%%%%
        %eigV=find(diag(Wz{z}.^2)>SNR); (0.27)
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
        %%%%%%%%%%%%%% FIND OPTIMAL Pz END (optinal) %%%%%%%%%%%%%%
    end
    
    
    if isempty(eigV) %if no optimal Pz
        eigV=1;
    end
    
    Uz_r{z}=Uz{z}(:,eigV); %subspace reduction
    Vz_r{z}=Vz{z}(:,eigV); %subspace reduction
    
    
    % Compute Bs, Bt, As and At with both whitening and subspace reduction
    Bs{z}=Fs*Uz_r{z};%Bz (0.24)
    Bt{z}=Ft*Vz_r{z};%Dz (0.24)
    As{z}=Gs* Uz_r{z}; %Az (0.24)
    At{z}=Gt* Vz_r{z}; %Ez (0.24)
    
    fullBs{z}=Fs*Uz{z}; %only for display purpose
    fullBt{z}=Ft*Vz{z}; %only for display purpose
    fullAs{z}=Gs* Uz{z}; %only for display purpose
    fullAt{z}=Gt* Vz{z}; %only for display purpose
end
%topoplot_components(fullAt, fullAs,[0 1],250,9,5,'d:\tmp\') %needed to
%check the components of the filter.
end
