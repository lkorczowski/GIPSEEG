% function adapted from the LiNGAM package
% complete software may be downloaded from http://www.cs.helsinki.fi/group/neuroinf/lingam/

function [Wopt,rowp] = permnozeribrutal ( W )

%--------------------------------------------------------------------------
% Try all row permutations, find best solution
%--------------------------------------------------------------------------

n = size(W,1);

bestval = Inf;
besti = 0;

allperms = perms(1:n);
nperms = size(allperms,1);

for i=1:nperms,
    Pr = eye(n); Pr = Pr(:,allperms(i,:));
    Wtilde = Pr*W;
    c = sum(1./diag(abs(Wtilde)));
%     c = nzdiagscore(Wtilde);
    if c<bestval,
	bestWtilde = Wtilde;
	bestval = c;
	besti = i;
    end
end

Wopt = bestWtilde;
rowp = allperms(besti,:);
rowp = iperm(rowp);


% function s = nzdiagscore( W )
% s = sum(1./diag(abs(W)));