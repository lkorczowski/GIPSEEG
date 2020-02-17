function [Ytest d C COVte] = mdm(COVte,COVtr,Ytraining,varargin)
    
    if isempty(varargin)
        method_mean = 'ld';
        method_dist = 'riemann';
    else
        method_mean = varargin{1};
        method_dist = varargin{2};
    end
    
    labels = unique(Ytraining);
    Nclass = length(labels);
    C = cell(Nclass,1);
    
    % estimation of center
    for i=1:Nclass
        %disp(['Compute mean cov for class' int2str(i) ' ...'])
        C{i} = mean_covariances(COVtr(:,:,Ytraining==labels(i)),method_mean);
        %disp('... DONE')
    end

    % classification
    NTesttrial = size(COVte,3);
    
    d = zeros(NTesttrial,Nclass);
    %disp(['Computing distance from mean Cov'])
    for j=1:NTesttrial
        for i=1:Nclass
            [d(j,i)] = distance(COVte(:,:,j),C{i}(:,:),method_dist);
        end
    end
    %disp('Classification DONE')
    
    [~,ix] = min(d,[],2);
    Ytest = labels(ix);