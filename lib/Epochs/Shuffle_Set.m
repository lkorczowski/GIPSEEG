function [Xout SHUF]=Shuffle_Set(X,Y,SHUF)
% X epochs EEG DATA
% Y epochs' label
% SHUF the permutation vector for X that schuffle the signals from the
% second player with respect to the class Y
%%

if nargin <3
 disp('Training/Test set generated chronologically')
 return
elseif isempty(SHUF); 
    SHUF=zeros(size(Y));
    C1=find(Y==1);
    C1rnd=derangement(C1);
    SHUF(C1)=C1rnd;
    C0=find(Y==0);
    C0rnd=derangement(C0);
    SHUF(C0)=C0rnd;
elseif length(SHUF)==length(Y)
else
    disp('RANDOM NUMBER NOT GENERATE, RND seed not matching with data set dimension')
    return
end

    Player2=size(X,1)/2; %2players
    Xout= X;
    Xout(Player2+1:2*Player2,:,:) = X(Player2+1:2*Player2,:,SHUF);