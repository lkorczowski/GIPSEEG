function Xnorm=normEEG(X,Method,Param)
% Xnorm=normEEG(X,Method,Param)
% X : 3D matrix of size [N,T,K] (space,time,repetition)
% Method :
%       'fro' : normalize each repetition for unitary trace
%       'baseline' : WARNING NOT FINISHED remove baselined based in Param
%                    Param being the time sampled used to compute the baseline
%                    example : Param=1:10 (ten first sample used for
%                    baseline)
%       '

if nargin<2
    Method='';
end

if strcmp(Method,'fro')
    
    for k=1:size(X,3)
        coef(k)=norm(X(:,:,k),'fro')/length(X(:,:,k));
        Xnorm(:,:,k)=X(:,:,k)/coef(k);
    end
    
elseif strcmp(Method,'baseline')
    if nargin<2 || isempty(Param)
        winBaseline=96:128; % HARDCODED NOT GOOD !!!!!!!!!!!!! ! ! !
    else
        winBaseline=Param;
    end
    for k=1:size(X,3)
        coef(:,:,k)=repmat(mean(X(:,winBaseline,k),2),1,size(X,2)); %compute the baseline
    end
            Xnorm=X-coef;
            
elseif strcmp(Method,'space')
  
      %Xnorm=(X)./repmat(max(abs(X),[],2), [1,size(X,2)]);
      Xnorm=(X)./repmat(max(abs(X),[],2), [1,size(X,2)]);

else
    % remove mean and var
    %mu=repmat(mean(X,2), [1,size(X,2)]);
      mu=repmat(mean(X,2), [1,size(X,2)]);
    Xnorm=(X-mu)./repmat(sqrt(var(X-mu,[],2)), [1,size(X-mu,2)]);
   % for k=1:size(Xnorm,3)
   %     Xnorm(:,:,k)=Xnorm(:,:,k)-repmat(mean(Xnorm(:,:,k),2),[1 size(Xnorm,2)]);
    %end
    
end