function [Xtr,Ytr,Xte,Yte,RND]=gp_gen_training_test_set(X,Y,ratio,RND)
% [Xtr,Ytr,Xte,Yte,RND]=gp_gen_training_test_set(X,Y,ratio,RND)
% stratified Monte Carlo Cross Validation
%
% Generate two sets for cross-validation
%
%% Syntaxe
% [Xtr,Ytr,Xte,Yte]=gp_gen_training_test_set(X,Y)
% [Xtr,Ytr,Xte,Yte]=gp_gen_training_test_set(X,Y,ratio)
% [Xtr,Ytr,Xte,Yte,RND]=gp_gen_training_test_set(X,Y,ratio,RND)

%% Details
% X : epochs EEG DATA
% Y : epochs' label
% ratio : if scalar : ratio training set data 0<P<1
% setBalance is 0 or 1)
% RND : is the permutation seed
[Xtr, Ytr,Xte, Yte]=deal([]);

if nargin<4 | RND==-1
    RND=1:length(Y);
    disp('generation of training test set is chronological')
    th=round(length(Y)*ratio);
    Xtr=X(:,:,1:th);Ytr=Y(1:th);
    Xte=X(:,:,th+1:end);Yte=Y(th+1:end);
else
    if size(Y,1)<size(Y,2)
        Y=Y';
    end
    if isempty(RND)
        RND=randperm(length(Y));
    end
    
    % randomization
    X=X(:,:,RND);
    Y=Y(RND);
    
    %for each class, find the epochs
    labels=unique(Y);
    Nclass=length(labels);
    indtr=[];
    indte=[];
    indicesC={};threshold={};
    for c=1:Nclass
        indicesC{c}=find(Y==labels(c));
        threshold{c}=round(length(indicesC{c})*ratio);
        indtr=[indtr; indicesC{c}(1:threshold{c})];
        indte=[indte; indicesC{c}(threshold{c}+1:end)];
    end
    
    Xtr=X(:,:,indtr);Ytr=Y(indtr);
    Xte=X(:,:,indte);Yte=Y(indte);
end