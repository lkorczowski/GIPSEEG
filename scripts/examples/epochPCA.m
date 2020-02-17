
%%
load('PCA8.mat')
for i=1:length(X)
for j=1:size(X{i},3)
Xout{i}(:,:,j)=V'*X{i}(:,:,j);
end
end