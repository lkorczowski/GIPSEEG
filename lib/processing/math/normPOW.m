function Xnorm=normPOW(X)
Xnorm=abs(X);
Xnorm=Xnorm-repmat((min(Xnorm,[],2)), [1,size(Xnorm,2)]);
Xnorm=Xnorm./repmat((max(Xnorm,[],2)), [1,size(Xnorm,2)]);