function [T labels]=gp_array2table(A)
% T=gp_array2table(A)
% Convert array of interger A into a binary table T in ascent order

labels=unique(A);
Nlabels=length(labels);
K=length(A);

T=zeros(K,Nlabels);
for indL=1:Nlabels
T(A==labels(indL),indL)=1;
end