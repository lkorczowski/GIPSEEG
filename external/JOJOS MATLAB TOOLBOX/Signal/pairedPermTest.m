function [h_out,p_out] = pairedPermTest(x,y,alpha,nbPerms)
% function  [h,p] = pairedPermTest(x,y,alpha,nbPerms)
%
% ************************************************************************
% Perform paired-sample statistical testing using randomization procedure.
% Null hypothesis: differences x(i)-y(i) drawn from distribution with mean 
% O and unknown variance.
%
% Input:
% ------
% x,y           --> vectors of paired observations with same lengths (1,N)
%                   can also be matrices of size (M,N), in which case M 
%                   separate tests are performed along each column of x,y 
%                   and a vector of results is returned.
% alpha     	--> (optional) significance level, (must be >= 2^-N)
% nbPerms     	--> (optional) performs either systematic or random data permutation
%                   'all': compute all permutations (2^N), or
%                	scalar: compute a random subset only, must be: (1/alpha < scalar < 2^N) 
%
% Output:
% -------
% h_out, p_out --> scalars with null hyp. rejection, and associated p-value
%
% History:
% --------
% *** 2013-06-12
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% -- assess inputs
if(size(x,2) == 1)
    x = x';
    y = y';
end
N       = size(x,2);
nbTests = size(x,1);
if (nargin<3) || (isempty(alpha))
    alpha = .05;
end
if (nargin<4) || (isempty(nbPerms))
    nbPerms = 'all';
end
if length(x(:)) ~= length(y(:))
   error('input data x and y must have equal length.');
end
if alpha < 2^-N 
    error('not enough possible permutations for this significance level.');
end

if ~strcmp(nbPerms,'all') && ((nbPerms<(1/alpha))||(nbPerms>2^N))
    error('invalid requested number of permutations.');
end


% -- create shuffling matrix
if strcmp(nbPerms,'all')
    nbPerms = 2^N;
end
permMat = uniqueShuffle(nbPerms,N);  % permutation matrix used to shuffle observations


% -- Perform t-test for OBSERVATION SET
[h_obs,p_obs] = ttest(x,y);

% -- Perform t-test for PERMUTATION SET
perm_set1   = zeros(N,nbTests);
perm_set2   = zeros(N,nbTests);
p_perm      = zeros(nbPerms,nbTests);
for perm_ix = 1:nbPerms
    % draw specific permutation using 'permMat' logical indexes
    if(nbTests == 1)
        perm_set1(permMat(perm_ix,:))  = x(permMat(perm_ix,:));
        perm_set1(~permMat(perm_ix,:)) = y(~permMat(perm_ix,:));
        perm_set2(~permMat(perm_ix,:)) = x(~permMat(perm_ix,:));
        perm_set2(permMat(perm_ix,:))  = y(permMat(perm_ix,:));
    else
        perm_set1(permMat(perm_ix,:),:)  = x(permMat(perm_ix,:),:);
        perm_set1(~permMat(perm_ix,:),:) = y(~permMat(perm_ix,:),:);
        perm_set2(~permMat(perm_ix,:),:) = x(~permMat(perm_ix,:),:);
        perm_set2(permMat(perm_ix,:),:)  = y(permMat(perm_ix,:),:);
    end
    % compute t-test on this permutation
    [h_temp,p_perm(perm_ix,:)] = ttest(perm_set1,perm_set2);
end

% -- Process significance level from permutation values
p_out   = zeros(1,nbTests);
p_perm  = sort(p_perm,1,'ascend');
for test_ix = 1:nbTests
    p_out(test_ix) = (1/nbPerms)*length(find(p_perm(:,test_ix) < p_obs(test_ix)));
    if(p_out(test_ix)==0)
        p_out(test_ix) = 1/nbPerms;     % case when no permutation gives lower p value than observation
    end;
end    
h_out = (p_out <= alpha);


end