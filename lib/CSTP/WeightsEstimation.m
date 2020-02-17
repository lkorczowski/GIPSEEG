function [W Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At)
% [W Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At)
% Weights estimation for a set of epochs X based on the estimation of the
% SNR
%                         || Xhat(k) ||_f
% W(k)=SNR(Xhat(k))=  -----------------------
%                     || X(k) - Xhat(k) ||_f
%
% where || x || is the frobenius norm of x
% Xhat is the estimation of the signal component of X based on the bilinear
% filter CSTP such as :
% Xhat(k)=As*Bs'*X(k)*Bt*At'; (0.8)
%
% INPUTS :
% X is a 3D EEG signal [ nb electrodes x nb samples x nb epochs]
% Y is the class' tags of the epochs of X [nb epochs]
% Bs Bt are respectively the spatial and temporal filters that project X(k)
% in the feature space
% As At are respectively the spatial and temporal filters that project back
% the features into the sensor space.
%
% OUTPUTS :
% W are the Weights estimation [nb epochs]
% Xhatk are the filtered data used for the signal estimation
%
% Note : The weights are normalized such as the sum of the weights of each
% class is equal to the number of epoch of the class
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also: applyWeights, CSTP, ACSTP

    W=ones(1,size(X,3));
    
    if nargin<2
        Y=ones(size(X,3),1);
    end
        Class=unique(Y); % find all the classes

    
if nargin<3 %%no CSTP weights (initialization)
    for k=1:size(X,3) % for each epoch
        W(k)=1/norm(X(:,:,k),'fro'); %compute the initilized weights 
    end
else
    Xhatk=zeros(size(X));
    for k=1:size(X,3) % for each epoch
        Xk=X(:,:,k);
        c=(Y(k)==Class); % class of the current epoch
        Xhatk(:,:,k)=As{c}*Bs{c}'*Xk*Bt{c}*At{c}'; % signal estimation
        W(k)=norm(Xhatk(:,:,k) ,'fro')/norm( Xk-Xhatk(:,:,k),'fro'); % SNR estimation
    end
    

    
end
    % Weights normalization
    for z=1:length(Class) % for each class
        Wz=W(Y==Class(z));
        W(Y==Class(z))=length(Wz)*Wz/sum(Wz); % normalize each class weights
    end
end