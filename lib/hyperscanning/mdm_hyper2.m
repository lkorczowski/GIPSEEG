function [Ytest d C COVte] = mdm_hyper(COVte,COVtr,Ytraining,NbUsers,Stats,P300_ref_orientation,varargin)
%[Ytest d C] = mdm_hyper(COVte,COVtr,Ytraining,NbUsers,Stats,P300_ref_orientation,method_mean,method_dist)
% include different possibility for classification with several users
%
%
if isempty(varargin)
    method_mean = 'riemann';
    method_dist = 'riemann';
else
    method_mean = varargin{1};
    method_dist = varargin{2};
end
if nargin<6
    P300_ref_orientation='commonP1';
end
if nargin<5
    Stats='all';
end
if nargin<4
    NbUsers=1;
end

labels = unique(Ytraining);
Nclass = length(labels);
C = cell(Nclass,1);

if strcmp(Stats,'intra')
    NbE=size(COVte,1)/NbUsers; %number of electrodes
    %estimation of the center of each player
    for pl=1:NbUsers
        Mapping=1+(pl-1)*NbE:NbE*pl;
    for i=1:Nclass
        fprintf(['Compute mean cov for class' int2str(i) ' and player' int2str(pl) ' ...'])
        C{i}(Mapping,Mapping) = mean_covariances(COVtr(Mapping,Mapping,Ytraining==labels(i)),method_mean);
        fprintf('... DONE\n')
    end
    end
        % classification
    NTesttrial = size(COVte,3);
    COVte=covariances_p300_clean_inter(COVte,NbUsers,P300_ref_orientation);
    d = zeros(NTesttrial,Nclass);
    disp(['Searching minimum distance for classification'])
    for j=1:NTesttrial
        for i=1:Nclass
            d(j,i) = distance_riemann_hyper2(COVte(:,:,j),C{i},NbE);
        end
    end
else
    % estimation of center
    for i=1:Nclass
        fprintf(['Compute mean cov for class' int2str(i) ' ...'])
        C{i} = mean_covariances(COVtr(:,:,Ytraining==labels(i)),method_mean);
        fprintf('... DONE\n')
    end
    
    % classification
    NTesttrial = size(COVte,3);
    
    d = zeros(NTesttrial,Nclass);
    disp(['Searching minimum distance for classification'])
    for j=1:NTesttrial
        for i=1:Nclass
            d(j,i) = distance_riemann2(COVte(:,:,j),C{i});
        end
    end
end
disp('Classification DONE')

[~,ix] = min(d,[],2);
Ytest = labels(ix);