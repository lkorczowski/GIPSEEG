% regularized LDA by shrinkage.
% inputs :  Ftest, test features
%           Ftrain, training features, Nvariables x Ntrials
%           Ytrain, training labels, Ntrials x 1
%
% outputs : class, test labels
%           y1, test scores
%           W, LDA coefficients
%           b, bias

function [class y1 W b s] = reglda(Ftest,Ftrain,Ytrain,METHOD_COV,ARG_COV)

if nargin <4
    METHOD_COV = 'shcov';
    ARG_COV = {};
end

if size(Ftrain,2)~=length(Ytrain)
    error('Features and labels must have the same length.');
end
if size(Ftrain,1)~=size(Ftest,1)
    error('Training and Test set must have the same number of variables.');
end

labels = unique(Ytrain);
Nclass = length(labels);

Nelec = size(Ftrain,1);

mu = zeros(Nelec,Nclass);
Covclass = zeros(Nelec,Nelec,Nclass);

for i=1:Nclass
    mu(:,i) = mean(Ftrain(:,Ytrain==labels(i)),2);
    Covclass(:,:,i) = covariances(Ftrain(:,Ytrain==labels(i)),METHOD_COV,ARG_COV);
end
  
mutot = mean(mu,2);
    
Sb = zeros(Nelec,Nelec);    
for i=1:Nclass
    Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
end
    
S = mean(Covclass,3);

[W Lambda] = eig(Sb,S);
[~, Index] = sort(diag(Lambda),'descend');

W = W(:,Index);

b = W(:,1)'*mutot;

s = sign(W(:,1)'*mu(:,2)-b);

y1 = s*(W(:,1)'*Ftest-b);
class = y1>0;