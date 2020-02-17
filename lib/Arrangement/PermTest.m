function p_val=PermTest(X,Y,nPerm,alpha)
% from Jonas Chatel-Goldman 2013

%--- set some parameters
if nargin<4
alpha   = 0.05;     % significance level
end
if nargin<3
nPerm   = 10000;    % number of permutations (i.e., size of surrogate data)
end
nObs    = length(X);     	% number of observations


vec_OBS_1   = X;
vec_OBS_2   = Y;

% --- perform permutation tests
if alpha < (1/nPerm)
    disp(['Not enough permutations for this significance level (' num2str(alpha) ')']);
    alpha = 1/nPerm;
    disp(['--> Significance level set to minimum possible (' num2str(alpha) ')']);
end
schuffle_index = uniqueShuffle2(nPerm, nObs, 1); % permutation matrix used to shuffle values between observation sets
diff_PERM = zeros(1,nPerm);
% compute OBSERVATION (reference/real value)


diff_OBS= mean(vec_OBS_2 - vec_OBS_1);  % here: average paired differences (CAN BE A SIMPLE DIFFERENCE, OR ANY CALCULATION!!!)
% compute PERMUTATION (surrogate values)
vec_PERM_1 = zeros(1,nObs);
vec_PERM_2 = zeros(1,nObs);
for perm_ix = 1:nPerm
    % draw specific permutation using 'schuffle_index' logical indexes
    vec_PERM_1(schuffle_index(perm_ix,:))   = vec_OBS_1(schuffle_index(perm_ix,:));
    vec_PERM_1(~schuffle_index(perm_ix,:))  = vec_OBS_2(~schuffle_index(perm_ix,:));
    vec_PERM_2(~schuffle_index(perm_ix,:))  = vec_OBS_1(~schuffle_index(perm_ix,:));
    vec_PERM_2(schuffle_index(perm_ix,:))   = vec_OBS_2(schuffle_index(perm_ix,:));
    % compute surrogate surrogate value for this specific permutationLrL
    diff_PERM(perm_ix) =  mean(vec_PERM_2 - vec_PERM_1);
end

% compute statistical significance and plot results
p_val = (1/nPerm)*(1+length(find(diff_PERM > diff_OBS)));
h_val = p_val < alpha;
%figure;dispBootstrap(diff_OBS, diff_PERM)