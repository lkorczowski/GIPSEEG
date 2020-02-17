function [D V]=eigsort(A,B)    
if nargin <2
    [V,D] = eig(A);
else
[V,D] = eig(A,B);
end
    [D,I] = sort(diag(D),'descend');
    V = V(:, I);
    
end