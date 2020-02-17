function [PARTIAL_CORR, CORR, PVAL] = procPartialCorr(X,CORR_TYPE)
% [PARTIAL_CORR, CORR, P] = procPartialCorr(X)
%
% This function computes partial correlation between pairs of variables in
% X, through inversion of the concentration matrix.
% THIS FUNCTION USES MATLAB STATISTICS TOOLBOX.
%
% Inputs:
% - X         	--> N-by-P matrix with rows corresponding to observations, 
%                   and columns corresponding to variables.
% - CORR_TYPE 	--> (optionnal) type of correlation to compute:
%                       'Pearson' -> (linear) partial correlations (default)
%                       'Spearman' -> (rank) partial correlations
%
% Outputs:
% - PARTIAL_CORR	--> sample linear partial correlation coefficients
%                   (symetric P-by-P matrix)
% - CORR            --> sample linear correlation coefficients
%                   (symetric P-by-P matrix)
% - PVAL            --> matrix of p-values for testing the hypothesis of no
%                       correlation against the alternative that there is a 
%                       nonzero correlation
%
% Last version:  24/07/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% default correlation type is Pearson
if(nargin < 2)  
    CORR_TYPE = 'Pearson';  
end

% computing partial correlation through concentration matrix (inverse of
% correlation matrix)
[CORR PVAL] = corr(X,'type',CORR_TYPE);
PARTIAL_CORR = pinv(CORR);
% PARTIAL_CORR = -corrcov(inv(CORR)); 
% must be replaced with implementation next line, since most of the 
% time input matrix is not semi-definite 
PARTIAL_CORR = -(diag(diag(PARTIAL_CORR.^(-1/2)))*PARTIAL_CORR*diag(diag(PARTIAL_CORR.^(-1/2)))); 





end
