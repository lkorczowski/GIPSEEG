function [R Beta1 Pval]=Regression2D(Y,X)
% [R C Pval]=Regression2D(Y,X) 
% Compute regression between X and the rows of Y
% Return the pearson corrcoef R, the slope Beta1 and the p-value Pval

for n=1:size(Y,1)
[Rtmp,P]  =corrcoef(Y(n,:),X,'alpha',1e-10);
R(n)=Rtmp(2);
Pval(n)=P(2);
Beta1(n)=(X-mean(X))'\(Y(n,:)-mean(Y(n,:)))';
MaxMin(n,:)=[min(Y(n,:)) max(Y(n,:))  max(Y(n,:))/min(Y(n,:))];
end

