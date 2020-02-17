function Pout=supressMean(P1)

Pout=P1-repmat(mean(P1,2),1,size(P1,2));

end
