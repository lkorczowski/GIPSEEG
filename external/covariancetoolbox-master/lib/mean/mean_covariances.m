% [] = XXXX()
%
% -----------------------------Definition---------------------------------%
%
% description
%
% usage :
%
% -----------------------------Input--------------------------------------%
%
% XX : description + dimension
%
% -----------------------------Output-------------------------------------%
%
% XX : description + dimension
%
% -----------------------------References---------------------------------%
%
% [1] : XXX
%
%
%   Project : BCI-EEG
%
%   author : A. Barachant
%   date : 2011-XXXX
%   version : 1.0
%   status : a terminer
%   CEA/GRENOBLE-LETI/DTBS
%
%   See also MEAN, MEDIAN, RIEMANN_MEAN, REIMANN_MEDIAN, RIEMANN_TRIMMED_MEAN.

% [EOF: XXX.m]

function [C crit niter] = mean_covariances(COV,method_mean,arg_mean)

if (nargin<2)||(isempty(method_mean))
    method_mean = 'arithmetic';
end
if (nargin<3)
    arg_mean = {};
end
crit=0;niter=0;
switch method_mean
    case 'riemann'
        C = riemann_mean(COV,arg_mean);
    case 'riemann_step'
        C = riemann_mean_stepwise(COV,arg_mean);
    case 'riemanndiag'
        C = riemann_diag_mean(COV,arg_mean);
    case 'riemanntrim'
        C = riemann_trimmed_mean(COV,arg_mean);
    case 'logdet'
        C = logdet_mean(COV);
    case 'median'
        C = median(COV,3);
    case 'riemannmed'
        C = riemann_median(COV,arg_mean);
    case 'logeuclid'
        C = logeuclid_mean(COV);
    case 'opttransp'
        [C, crit, niter] = opttransp_mean(COV);
    case 'wasserstein'
        C = mean_wasserstein(COV);
    case 'ws'
        C = mean_wasserstein(COV);
    case 'power'
        [C, crit, niter] = mean_power(COV,arg_mean);
        
    case 'power_dual'
        [C, crit, niter] = mean_power_dual(COV,arg_mean);
    case 'congedo'
        [C, crit, niter] = ws_congedo_mean(COV);
    case 'ld'
        [C, crit, niter] = logdet_mean(COV);
    case 'geodesic'
        C = geodesic_mean(COV,arg_mean);
    case 'harmonic'
        iCOV = zeros(size(COV));
        for i=1:size(COV,3)
            iCOV(:,:,i) = inv(COV(:,:,i));
        end
        C = inv(mean(iCOV,3));
        
    case 'geometric'
        B = COV(:,:,1);
        for i=2:size(COV,3)
            B = B*COV(:,:,i);
        end
        C = B^(1/size(COV,3));
    otherwise
        C = mean(COV,3);
end
