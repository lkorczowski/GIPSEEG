%% SPEARMAN or KENDALL test for independence
% null hypothesis: uncorrelation
% if p<alpha, reject -> presence of correlation

% input: X row column
% tipo: 'Spearman', 'Kendall'

function [pM,rho,contarigett,stringout,stringflag]=test_independence(X,tipo,alpha);

M=size(X,1); % n. of signals
N=size(X,2); % length of signals

[rho,pM] = corr(X','type',tipo);

contarigett=0;
for i=1:M
    for j=i+1:M
        if pM(i,j)<alpha
            contarigett=contarigett+1;
        end
    end
end

if contarigett == 0
    stringout='non-rejection: all signals are INDEPENDENT';
    stringflag=0;
else
    stringout='rejection: at least two signals are NOT INDEPENDENT';
    stringflag=1;
end


