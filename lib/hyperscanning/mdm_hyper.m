function [Ytest d C COVte EIG] = mdm_hyper(COVte,COVtr,Ytraining,NbUsers,Stats,P300_ref_orientation,varargin)
%[Ytest d C] = mdm_hyper(COVte,COVtr,Ytraining,NbUsers,Stats,P300_ref_orientation,method_mean,method_dist)
% include different possibility for classification with several users
%
%P300_ref_orientation : 'commonP1', 'multiP1', 'noP1'
if isempty(varargin)
    method_mean = 'riemann';
    method_dist = 'riemann';
    method_dist_eig = 0;
else
    method_mean = varargin{1};
    method_dist = varargin{2};
    method_dist_eig = varargin{3};
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

% SELECT TYPE OF STATISTICS USED FOR CLASSIFICATION
if strcmp(Stats,'MDM-multi') %ONLY INTRA STATS
    NbE=size(COVte,1)/NbUsers; %number of electrodes
    
    %estimation of the center of class for each player
    for pl=1:NbUsers
        Mapping=1+(pl-1)*NbE:NbE*pl;
        for i=1:Nclass
            %fprintf(['Compute mean cov for class' int2str(i) ' and player' int2str(pl) ' ...'])
            C{i}(Mapping,Mapping) = mean_covariances(COVtr(Mapping,Mapping,Ytraining==labels(i)),method_mean);
            %fprintf('... DONE\n')
        end
    end
    
    % classification
    NTesttrial = size(COVte,3);
    COVte=covariances_p300_clean_inter(COVte,NbUsers,P300_ref_orientation);
    d = zeros(NTesttrial,Nclass);
    %disp(['Searching minimum distance for classification'])
    
    for j=1:NTesttrial
        for i=1:Nclass
            [d(j,i) EIG{j,i}] = distance_riemann_hyper(COVte(:,:,j),C{i},NbE,method_dist_eig);
            %d(j,i) = distance_riemann(COVte(:,:,j),C{i});
        end
    end
    
elseif strcmp(Stats,'MDM')
    NbE=size(COVte,1)/NbUsers; %number of electrodes
    
    %estimation of the center of class for each player
    for pl=1:NbUsers
        Mapping=1+(pl-1)*NbE:NbE*pl;
        for i=1:Nclass
            %fprintf(['Compute mean cov for class' int2str(i) ' and player' int2str(pl) ' ...'])
            C{i,pl} = mean_covariances(COVtr(Mapping,Mapping,Ytraining==labels(i)),method_mean);
            %fprintf('... DONE\n')
        end
    end
    
    % classification
    NTesttrial = size(COVte,3);
    COVte=covariances_p300_clean_inter(COVte,NbUsers,P300_ref_orientation);
    %d = zeros(NTesttrial,Nclass);
    %disp(['Searching minimum distance for classification'])
    for pl=1:NbUsers
        Mapping=1+(pl-1)*NbE:NbE*pl;
        for j=1:NTesttrial
            for i=1:Nclass
                [d{pl}(j,i) EIG{pl}{j,i}] = distance_riemann(COVte(Mapping,Mapping,j),C{i,pl},NbE);
                %d(j,i) = distance_riemann(COVte(:,:,j),C{i});
            end
        end
    end
    
else %MDM-hyper - INTRA+INTER STATS (interpersonnal)
    
    % estimation of center of each class
    for i=1:Nclass
        %fprintf(['Compute mean cov for class' int2str(i) ' ...'])
        C{i} = mean_covariances(COVtr(:,:,Ytraining==labels(i)),method_mean);
        %fprintf('... DONE\n')
    end
    
    % classification
    NTesttrial = size(COVte,3);
    
    d = zeros(NTesttrial,Nclass);
    %disp(['Searching minimum distance for classification'])
    for j=1:NTesttrial
        for i=1:Nclass
            %d(j,i) = distance_riemann(COVte(:,:,j),C{i},method_dist_eig);
            [d(j,i) EIG{j,i}] = distance_riemann_hyper(COVte(:,:,j),C{i},method_dist_eig);
        end
    end
end
%disp('Classification DONE')
if strcmp(Stats,'MDM')
    [~,ix]=cellfun(@(x) min(x,[],2), d,'UniformOutput', false);
    for i=1:length(ix)
        Ytest{i}=labels(ix{i});
    end
else
    [~,ix] = min(d,[],2);
    Ytest = labels(ix);
end