function [COV P1] = covariances_p300_hyper(X,Y,NbUsers,P300_ref_orientation,method_cov,arg_cov)
% [COV P1] = covariances_p300_hyper(X,Y,users,P300_ref_orientation,method_cov,arg_cov)
%   compute the covariance matrix of EEG epoch X with the P1 reference added
% [COV] = covariances_p300_hyper(X,P1,users,P300_ref_orientation,method_cov,arg_cov)
%   compute the covariance matrix of EEG epoch X with the P1 reference added
%
% INPUTS
%       X :     matrix of epoch of size [Nbe*NbUsers x Nbs x Nbt] where Nbe is the
%               number of electrodes, NbUsers is the number of users, Nbs is the number of samples per
%               trials, Nbt is the number of trials
%               Thus it include recording of several users synchronized.
%
%       Y :     if [1xNbt] vector target class for each trial of X
%               if [Nbe x Nbs] matrix P1 of reference (usually the average
%               of P300 response of the (or all) subjet(s)
%
%       NbUsers : double, number of users
%
%       P300_ref_orientation :
%               'commonP1'     : COV is computed with one common
%                   reference computed with the mean of P300 in X of all
%                   players
%                   Thus
%                   P1 is [Nbe x Nbs]
%                   COV is [Nbe*(NbUsers+1) x Nbe*(NbUsers+1)] and computed
%                   from [P1;X]
%               'multiP1'   : COV is computed with one P1 for each user
%                   Thus
%                   P1 is [Nbe x Nbs x NbUsers]
%                   COV is [Nbe*NbUsers*2 x Nbe*NbUsers*2]
%               'noP1' : COV is computed only with data
%                   Thus
%                   COV is [Nbe*NbUsers x Nbe*NbUsers]
%
% COV

if nargin < 6
    arg_cov = {};
end
if nargin < 5 || isempty(method_cov)
    method_cov = 'scm';
end
if nargin < 4
    P300_ref_orientation = 'commonP1';
end
if nargin < 3
    NbUsers=1;
end

% BEGINS...
%if NbUsers==1 %if only 1 player, just use standard function
%    [COV P1] = covariances_p300(X,Y,method_cov,arg_cov);
    
%else %if more player
    
    %PHASE 1 : COMPUTE OR NOT THE AVERAGE P300 (TRAINING SET)
    Nbe=size(X,1)/NbUsers; % compute the number of electrods
    
    if all(size(Y(:,:,1))==[Nbe size(X(:,:,1),2)]) %check if Y is the P1 or the labels' vector 
                P1 = Y; % if case the input Y was the P300 (TEST SET)

            else

        switch P300_ref_orientation
            case 'commonP1' %for common P1, average between the P300 of all players
                P1 = mean(X(:,:,Y==1),3);
                
                P1=transpose(mean(reshape(P1',size(P1',1),Nbe,NbUsers),3));
                %mean P1 on all users
            case 'multiP1' %compute the average P300 for every players
                for Us=1:NbUsers
                    P1(:,:,Us) = mean(X((Nbe*(Us-1)+1):(Nbe*(Us)),:,Y==1),3);
                end
            case 'noP1'
                P1=[];
        end
        
    end
    
    % PHASE 2 : COMPUTE THE COVARIANCE MATRIX WITH SPECIFIC P1 ORGANIZATION
    if size(P1,3)>1 % CASE for P1_orentation is 'multiP1'
        %X2 = zeros(size(X,1)+size(P1,1)*size(P1,3),size(X,2),size(X,3));
        for i=1:size(X,3)
            tmp=[];
            for Us=1:NbUsers
               
                tmp=cat(1,tmp,P1(:,:,Us),X((Nbe*(Us-1)+1):(Nbe*(Us)),:,i));
            end
            X2(:,:,i) = tmp;
        end
    elseif isempty(P1)
        X2=X;
        
    else % CASE for P1_orientation is 'commonP1'
        X2 = zeros(size(X,1)+size(P1,1),size(X,2),size(X,3));
        for i=1:size(X2,3)
            X2(:,:,i) = cat(1,P1,X(:,:,i));
        end
    end
    
    % Covariance set
    COV = covariances(X2,method_cov,arg_cov);
%end