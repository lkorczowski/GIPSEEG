function Vred=mypca(X,B,PLOT)
% example
% load('D:\data\summary.mat')
% 
% P1all=mean(P1(:,:,:),3);
% P1all=P1all(randperm(size(P1all,1)),:);
% 
% V=mypca(P1all,8,1) %or V=mypca(P1all,0.999,1)
% 
% plot((V'*P1all)')


if nargin<3
    PLOT=0;
end
if nargin<2 || B==1
    B=0.999;
end


COV=cov(X');
[V, D]=eig(COV);
[Y I]=sort(diag(D),'descend');

if PLOT
    figure
    plot(cumsum(Y)./sum(Y))
    xlabel('nb parameters')
    ylabel('cumsum')
end

if B<1
IND=find((cumsum(Y)./sum(Y))>B,1);
else
IND=B;
end

Vred=V(I(1:IND),:);