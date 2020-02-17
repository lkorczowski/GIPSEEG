% function adapted from the LiNGAM package
% complete software may be downloaded from http://www.cs.helsinki.fi/group/neuroinf/lingam/


function [Bopt,optperm,bestval] = permslowertriagbrutal( B )

% questa è per strictly lower triangular su B0!

% Finds the best identical row and column permutation in terms of getting
% approximately strictly lower triangular matrix
    
%--------------------------------------------------------------------------
% Try all permutations, find best solution
%--------------------------------------------------------------------------

n = size(B,1);

bestval = Inf;
besti = 0;
allperms = perms(1:n);
nperms = size(allperms,1);

for i=1:nperms,
    Btilde = B(allperms(i,:),allperms(i,:));
    
    D = tril(Btilde,-1)-Btilde;
    c = sum(sum(D.^2));
%     c = sltscore(Btilde);

    if c<bestval,
	bestBtilde = Btilde;
	bestval = c;
	besti = i;
    end
end

Bopt = bestBtilde;
optperm = allperms(besti,:);

% % % % % function s = sltscore( B )
% % % % % % sltscore - strictly lower triangular score
% % % % % D = tril(B,-1)-B;
% % % % % s = sum(sum(D.^2));
