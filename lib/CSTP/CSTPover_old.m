function [Bs Bt As At Class eigV]=CSTP(X,Y,P,PLOT)
% OLD VERSION, TO REMOVE
% [Bs Bt As At]=CSPT(X,Y,PLOT)
% Compute Common Spatio-Temporal Pattern for 2 classes to maximize the SNR
% with Xhat=As*Bs'*X*Bt*At'
% Xhat is the estimated signal
% Bs is the spatial filter (features)
% Bt is the temporal filter (features)
% As and At allows to go back in the sensors spaces
%
if nargin<4
    PLOT=0;
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
eigWhite=1:size(X,1);
eigV=1:12;%reduction
if nargin<2 || isempty(Y)
    Y=ones(size(X,3),1);
end
Class=unique(Y);
nbClass=length(Class);
%PLOT PARAMETERS
if PLOT==1
Scale=8;
fs=128;
LABELS=[];
end

%%%%%%%%%%%%%%%% START COMPUTING CSPT %%%%%%%%%%%%%%%%%%%%%%%%




%compute Sample covariance matrices
for k=1:size(X,3)
    Cs_k(:,:,k)=X(:,:,k)*X(:,:,k)'/(size(X,1)-1);
    Ct_k(:,:,k)=X(:,:,k)'*X(:,:,k)/(size(X,2)-1);
end

%mean covariance matrices
Cs=[];
Ct=[];
for z=1:nbClass%find(Class==1) %%%%%%%%%%%%%%%%%%%% WARNING
    if Bool==1
    P{z}=mean(X(:,:,Y==Class(z)),3); % WITH NO OVERLAPPING
    end
Cs=cat(3,Cs,mean(Cs_k(:,:,Y==Class(z)),3)); %spatial signal + noise covariance matrix for each class

Ct=cat(3,Ct,mean(Ct_k(:,:,Y==Class(z)),3)); %temporal signal + noise covariance matrix for each class
end
 % mean covariance matrix
 if nbClass>1
Cs=mean(Cs,3);
Ct=mean(Ct,3);
 end
 
 %compute expected SNR ratio
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
SNR=SNRnum/SNRden(2);

[Ux Ws]=eig(Cs);
%norm(Cs-Ux*Ws*Ux','fro') %to check the estimate
[Vx Wt]=eig(Ct);
%norm(Ct-Vx*Wt*Vx','fro')  %to check the estimate
[Ws,ind] = sort(diag(Ws),'descend');
if any(Ws<=0); disp('Warning negative SPATIAL eigenvalue(s)'); end
Ux = Ux(:,ind);
Ux = Ux(:,eigWhite);
[Wt,ind] = sort(diag(Wt),'descend');
if any(Wt<=0); disp('Warning negative TEMPORAL eigenvalue(s)'); end


Vx = Vx(:,ind);
Vx = Vx(:,eigWhite);
Fs = (sqrt(pinv(diag(Ws(eigWhite)))) * Ux')'; 
Ft =(sqrt(inv(diag(Wt(eigWhite)))) * Vx')';
Gs =Ux*sqrt((diag(Ws(eigWhite)))) ;
Gt =Vx*sqrt((diag(Wt(eigWhite)))) ;
%{
disp('Fs_T*Gs');
diag(Fs'*Gs);
disp('Ft_T*Gt');
diag(Ft'*Gt);
%}
%plot normal P300 and power of covariance matrixes
%{
figure 
subplot(231);imagesc(log(abs(Cs)));title('Cs')
  subplot(232);  imagesc(log(abs(Ct)));title('Ct')
    subplot(233);  imagesc(log(abs(Fs)));title('Fs')

      subplot(234);  imagesc(log(abs(Ft)));title('Ft')
  subplot(235);  imagesc(log(abs(Gs)));title('Gs')

    subplot(236);  imagesc(log(abs(Gt)));title('Gt')
%}
if PLOT==1

figure;subplot(141)
plotEEG(P{1},Scale,fs,LABELS);title('X^{bar}_0')
subplot(142)
plotEEG(P{2},Scale,fs,LABELS);title('X^{bar}_1')
end
%compute CSPT
for z=1:length(P)
    Z{z}=Fs'*P{z}*Ft;
    [Uz{z},Wz{z},Vz{z}]=svd((Z{z}));
    %check
    %{
    disp('UzT*Z*Vz')
        Uz{z}'*Z{z}*Vz{z}
    %}
    eigV=find(diag(Wz{z})>SNR);
    Uz_r{z}=Uz{z}(:,eigV);
    Vz_r{z}=Vz{z}(:,eigV);
    
     %Uz_r{z}'*Uz{z}*Wz{z}*Vz{z}'*Vz_r{z}

    
    Bs{z}=(Uz_r{z}'*Fs')';
    Bt{z}=(Vz_r{z}'*Ft')';
    As{z}=Gs*Uz{z}(:,eigV);
        At{z}=Gt*Vz{z}(:,eigV);
    %check
    %{
        disp('BsT*Cs*Bs')
    Bs{z}'*Cs*Bs{z} %I
            disp('BtT*Ct*Bt')

    Bt{z}'*Ct*Bt{z} %I
    Bs{z}'*P{z}*Bt{z} %W
    %}
    %X_hat{z}=As{z}*Bs{z}'*P{z}*Bt{z}*At{z}';
    
end

if PLOT==1
subplot(143)
plotEEG(X_hat{1},Scale,fs,LABELS);title('X^{bar hat}_0 (CSPT)')
subplot(144)
plotEEG(X_hat{2},Scale,fs,LABELS);title('X^{bar hat}_1 (CSPT)')
eegplot([P{1} P{2} X_hat{1} X_hat{2}],'srate',128)
end