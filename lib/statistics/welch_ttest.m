function [tval, pval, df,sd,ci]=welch_ttest(x,y,varargin)
% [t, pval, df]=welch_ttest(x,y,varargin)
% perform Welch's ttest or unequal variances t-test, is a two-sample
% location test which is used to test the hypothesis that two populations have equal means.
% Better than ttest when the two samples have unequal variances and unequal sample sizes.
% Hypothesis1 : sig1~=sig2
% Hypothesis2 : K1+K2>40 or normal distribution
%
% OUTPUTS :
% t : t-value, pval : pvalue, df: degree of freedom
% sd, ci : std and confidence interval of x and y


% initialization
options = struct('alpha',0.05);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('plotEEGvariability needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = lower(pair{1}); % make case insensitive
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end
if ~exist('ttest2')
    mu1=mean(x);
    mu2=mean(y);
    sig1=std(x-mu1);
    sig2=std(y-mu2);
    K1=length(x);
    K2=length(y);
    
    % t-value
    tval=(mu1-mu2)/sqrt(sig1^2/K1 + sig2^2/K2);
    %Welch–Satterthwaite equation to estimate the degree of freedom
    df=(sig1^2/K1+sig2^2/K2)^2/...
        (sig1^4/(K1-1)/K1^2 + sig2^4/(K2-1)/K2^2);
    
    pval=NaN;
    error('welch_ttest : function not fully implemented, please use ttest2(x, y, ''Vartype'', ''unequal'') instead')
    
else
    
    [h,pval,ci,stats] = ttest2(x, y, 'Vartype', 'unequal','alpha',options.alpha);
    tval=stats.tstat;
    df=stats.df;
    sd=stats.sd;
end


