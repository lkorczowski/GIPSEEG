function [C,A,info] = simRandDat(cfg)

% function simulating covariance data for AJD and barycenter Algorithms.
%
% The variance of the sources is assumed to decay with a power law.
% The covariance matrices also contain additive structured noise and small
% white noise.
%
% The model is:
%   C_k = Asig*Dsig_k*Asig' + alpha*I + Anoise_k*Dnoise_k*Anoise_k'
% where Asig is the mixing matrix of the signal, Dsig_k are diagonal
% matrices containing the variance of the sources, alpha is the white noise
% level, I is the identity matrix, Anoise_k are mixing matrices of the
% structured noise (which are common for some matrices and different for
% others) and Dnoise_k are diagonal matrices containing the variance of the
% structured noise.
%
% Input: 
%   - cfg: struct containing fields   
%           - nChan:           number of channels                       (4)
%           - nMat:            number of matrices                      (20)
%           - p:               power law order                        (1.5)
%           - nCommonNoiseMat: number of matrices where the structured 
%                              noise has the same mixing           (nMat/2)
%           - cond_sig:        range of acceptable condition number for the
%                              true mixing matrix                       (1)
%           - cond_noise:      range of acceptable condition number for the
%                              structured noise mixing matrices   ([0,Inf]) 
%           - SNR_white:       signal-noise ratio relative to the
%                              white noise                            (Inf)
%           - SNR_struct:      signal-noise ratio relative to the
%                              structured noise                       (Inf)
%
% Output:
%   - C:    nChan*nChan*nMat array containing the covariance matrices
%   - A:    nChan*nChan mixing matrix of the signal
%   - info: struct containing fields
%           - d_sig:   nChan*nMat array containing variances of the sources
%           - var_sig: expected variance of the sources
%           - alpha:   white noise level
%           - A_noise: nChan*nChan*(1+nMat-nCommonNoiseMat) array
%                      containing the mixing matrices of the structured
%                      noise
%           - d_noise: nChan*nMat array containing variances of the
%                      structured noise
%           - cfg:     the cfg structure used to simulate data
%
% authors: Marco Congedo, Florent Bouchard @Gipsa
%
% last modification : 12/07/2015

% deal whith input
% default parameters generation (copy/past to quick generation of your own parameters)
cfg_ = struct('nChan',4,'nMat',20,'p',1.5,'nCommonNoiseMat',10,...
       'cond_sig',1,'cond_noise',[0,inf],'SNR_white',inf,'SNR_struct',inf);

% read the acceptable names
optionNames = fieldnames(cfg_);

if  nargin<1 || isempty(cfg)
        warning('default cfg')
else
    nameArgs=fieldnames(cfg);
    nArgs = length(nameArgs);
    
    for indO = 1: nArgs
        inpName = nameArgs{indO};
        if any(strcmp(inpName,optionNames))
            cfg_.(inpName) = cfg.(inpName);
        else
            warning('%s is not a recognized parameter name',inpName)
        end
    end
end

% initialize variables
N          = cfg_.nChan;
K          = cfg_.nMat;
p          = cfg_.p;
L          = cfg_.nCommonNoiseMat;
cond_sig   = cfg_.cond_sig;
cond_noise = cfg_.cond_noise;
SNR_w      = cfg_.SNR_white;
SNR_s      = cfg_.SNR_struct;

% generate signal mixing matrix
if cond_sig==1
    A = orth(randn(N));
else
    A = randn(N);
    while cond(A)<cond_sig(1) && cond(A)>cond_sig(2)
        A = randn(N);
    end
end

% generate sources
d_sig = zeros(N,K);
for k=1:K
    for n=1:N
        d_sig(n,k) = randn^2/(n^p);
    end
end

% get expected variance of sources
var_sig = sum(1./((1:N).^p));

% white noise parameter (common to all matrices)
alpha = var_sig/N/SNR_w;

% structured noise
% generate noise mixing matrices
A_noise = zeros(N,N,K-L+1);
for k=1:K-L+1
    if cond_noise==1
        A_noise(:,:,k) = orth(randn(N));
    else
        while cond(A_noise(:,:,k))<cond_noise(1) && cond(A_noise(:,:,k))>cond_noise(2)
            A_noise(:,:,k) = randn(N);
        end
    end
end
% generate structured noise components
d_noise = zeros(N,K);
for k=1:K
    for n=1:N
        d_noise(n,k) = randn^2/(n^p)/SNR_s;
    end
end

% generate final matrices
for k=1:K
    if k<=L
        C(:,:,k) = A*diag(d_sig(:,k))*A' + alpha*eye(N) + A_noise(:,:,1)*diag(d_noise(:,k))*A_noise(:,:,1)';
    else
        C(:,:,k) = A*diag(d_sig(:,k))*A' + alpha*eye(N) + A_noise(:,:,k-L+1)*diag(d_noise(:,k))*A_noise(:,:,k-L+1)';
    end
end

% fill info struct
info.d_sig   = d_sig;
info.var_sig = var_sig;
info.alpha   = alpha;
info.A_noise = A_noise;
info.d_noise = d_noise;
info.cfg     = cfg_;

