%% test all deferent cohence calculus
clear all
close all
% define parameters
Len=128;
Elec=3;
NoiseVar=0.5;
SignalRatio=1;

% let's have two multimentionnals signals x and y
t=(1:Len)/(Len)*10;
x=NoiseVar*randn(Elec,Len);x(1,:)=x(1,:)+sin(2*pi*t);
y=NoiseVar*randn(Elec,Len);y(1,:)=y(1,:)+SignalRatio*sin(2*pi*t);
figure;subplot(211);plot(x');subplot(212);plot(y')

% we can compute the covariance matrix of x and y such as
Cx=cov(x');
Cy=cov(y');
% Thus a matrix containing both covariance matrix such as :
% Cref = [Cxx  0 ] 
%        [ 0 Cyy ]
Cref=blkdiag(Cx,Cy);
% We can also have the covariance matrix of the two concatenate signal
% C = [Cxx Cxy]
%     [Cyx Cyy]
C=cov([x;y]');

% that contains the cross-variance terms off-diag :
Covxy=1/(size(x,2)-1)*supressMean(x)*supressMean(y)';
Covyx=1/(size(x,2)-1)*supressMean(y)*supressMean(x)';

% the cross-variance terms could be normalized such as
R=(Cx)^(-1/2)*Covxy*(Cy)^(-1/2);
% the synchronization between the signal x and y is
SYNC=SYNC_Riem(x,y)
% that is equivalent to the riemannian distance between C and Cref:
% dr(C,Cref) :
D_C_Cref=distance(Cref,C,'riemann')

Rxy=Cref^(-1/2)*C*Cref^(-1/2);
D_R_I=distance(Cref^(-1/2)*C*Cref^(-1/2),Cref^(-1/2)*Cref*Cref^(-1/2),'riemann')
% that is equivalent to
% the normalized covariance matrix :
% Rxy = [ I  R ]
%       [ R* I ]
D_R_I=distance(Rxy,eye(size(Rxy)),'riemann')

% that is equivalent to
SumEig=sqrt(sum(log(eig(Rxy)).^2))
MeanMeanAbsCovxy=mean(mean(abs(Covxy)))
SumSumAbsCovxy=sum(sum(Covxy))
NormR=norm(R,'fro')
NormFroRxy=norm(Rxy-eye(size(Rxy)),'fro')
SumSumR=sqrt(sum(sum(R.^2)))
SumSumFro=sqrt(sum(sum((Rxy-eye(size(Rxy))).^2)))

for i=1:size(x,1)
    for j=1:size(y,1)
        CORRcoef(i,j)=xcorr(supressMean(x(i,:)/var(x(i,:))),supressMean(y(j,:)/var(y(i,:))),0);
    end
end
MeanMeanAbsCORRCoef=mean(mean(abs(CORRcoef)))

%% compute correlation matrix

eig(Covxy)



%%

% let's have two multimentionnals signals x and y
[x r m]=zscore(x,2);
[y r m]=zscore(y,2);


% we can compute the covariance matrix of x and y such as
Cx=cov(x');
Cy=cov(y');
% Thus a matrix containing both covariance matrix such as :
% Cref = [Cxx  0 ] 
%        [ 0 Cyy ]
Cref=blkdiag(Cx,Cy);
% We can also have the covariance matrix of the two concatenate signal
% C = [Cxx Cxy]
%     [Cyx Cyy]
C=cov([x;y]');

% that contains the cross-variance terms off-diag :
Covxy=1/(size(x,2)-1)*supressMean(x)*supressMean(y)';
Covyx=1/(size(x,2)-1)*supressMean(y)*supressMean(x)';

% the cross-variance terms could be normalized such as
R=(Cx)^(-1/2)*Covxy*(Cy)^(-1/2);

% the synchronization between the signal x and y is
SYNC_Riem(x,y)
% that is equivalent to the riemannian distance between C and Cref:
% dr(C,Cref) :
distance(C,Cref,'riemann')

Rxy=Cref^(-1/2)*C*Cref^(-1/2);
distance(Cref^(-1/2)*C*Cref^(-1/2),Cref^(-1/2)*Cref*Cref^(-1/2),'riemann')
% that is equivalent to
% the normalized covariance matrix :
% Rxy = [ I  R ]
%       [ R* I ]
distance(Rxy,eye(size(Rxy)),'riemann')

% that is equivalent to
sqrt(sum(log(eig(Rxy)).^2))

%sum sum
