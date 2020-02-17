function Xrnd=derangement(X)
    BOOL=1;
    while BOOL
    Xrnd=X(randperm(length(X)));
    BOOL=any(Xrnd==X);
    end