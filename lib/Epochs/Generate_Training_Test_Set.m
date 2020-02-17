function [Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(X,Y,P,RND,NbMax,setBalance)
% [Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(X,Y,P,RND,NbMax,setBalance)
% X : epochs EEG DATA
% Y : epochs' label
% P : if scalar : ratio training set data 0<P<1 OR P is [1 x 2] where P(1) is the number of training TA
% epochs and P(2) is the number of testing TA epochs (NT fit accordingly if
% setBalance is 0 or 1)
% RND : is the permutation seed, if it is empty, a new seed will be generated with
% the parameters NbMax and setBalance. 
% NbMax : is the maximum of samples used for the final Training+Test set
% (before balancing)
% setBalance : if 1 all the classes are balanced (same number of epoch) if
% 0, we keep the proportion in X. Not accuratly working if setBalance=1.
%
% *** History: 2014-06-01
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work:  L. KORCZOWSKI, M. CONGEDO, C. JUTTEN "Single-Trial Classification of Multi-User P300-Based
% Brain-Computer Interface Using Riemannian Geometry" (IEEE EMBC, 2015)
% 
%%
if nargin<6 || isempty(setBalance)
    setBalance=0;
end
if nargin<5 || isempty(NbMax)
    NbMax=length(Y);
end
if nargin <4
    disp('Training/Test set generated chronologically')
    RND='none';
elseif isempty(RND); %generate new training-test random seed
    %    Y=Ys(1:20)
    RND=randperm(length(Y));
    if NbMax<length(Y)
        RND=RND(1:NbMax);
    end
    X = X(:,:,RND);
    Y = Y(RND);
    
elseif length(RND)<=length(Y); % use existing random seed
    X = X(:,:,RND);
    Y = Y(RND);
else
    error('RANDOM NUMBER NOT GENERATE, RND seed not matching with data set dimension')
end

    C1=find(Y==1);%indice of class1
    C0=find(Y==0);%indice of class2
    n1=length(C1);
    if setBalance, n0=n1;else,n0=length(C0);end
    
if length(P)==1 %compute the amount for each classes
    nb1tr=round(n1*P);
    nb1te=round(n1*(1-P));
    nb0tr=round(n0*P);
    nb0te=round(n0*(1-P));
else
    nb1tr=P(1);
    nb1te=P(2);
    if setBalance
        nb0tr=P(1);
        nb0te=P(2);
    else %assuming there is 
        nb0tr=P(1)*round(n0/n1);
        nb0te=P(2)*round(n0/n1);
    end
end
    
    TrainInd=sort([C1(1:nb1tr) ;C0(1:nb0tr)]);
    TestInd=sort([C1(end-nb1te+1:end) ;C0(end-nb0te+1:end)]);

    Xtr=X(:,:,TrainInd);
    Xte=X(:,:,TestInd);
    Ytr=Y(TrainInd);
    Yte=Y(TestInd);

end

