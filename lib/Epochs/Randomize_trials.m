function [Xout y RND]=Randomize_trials(X,Y,RND)
    Nbclass=length(Y);
    y=[];
    Xout=[];
    for i=1:Nbclass
        y=[y;ones(size(X{i},3),1)*Y(i)];
        Xout=cat(3,Xout,X{i});
    end
    if nargin<3
    RND=randperm(size(Xout,3));
    end
    Xout=Xout(:,:,RND);
    y=y(RND);
end