function [Xout labels]=MultiSet_generator(varargin)
% Generator of multisubjets sets of EEG
% Xout=Multisubject_generator(X1,Y1,X2,Y2,.. ,Xn, Yn)
% INPUT :
%         X1, X2, .. Xn : [nb_electrodes, nb_samples, nb_trials] matrix
%         Y1, Y2, .. Yn : [nb_trials] array giving the class of each trial for the given Xi
%     
% OR
% Xout=Multisubject_generator(X,Y) where X and Y are cells of [1xNbsubjects] containing Xi and Yi         
% 
% OUTPUT :
%         XOUT [1xNbClass] cell
%         labels [1xNbClass] cell
%     
% see also epoch_p300

if mod(nargin,2)==1
    disp('Multisubject_generator needs both EEG set and class list')
    return
end
Xte=varargin((1:2:nargin));% extract Multiset
Yte=varargin((2:2:nargin));
%% generate multiplayers set
[trash IndF]=min(cellfun('length',Yte));

Yout=Yte{IndF};
 labels = unique(Yout);
 Nclass = length(labels);
    IndC = cell(Nclass,1);%index of the class
    
    % estimation of center
Xout={};
    for c=1:Nclass
        Xout{c}=[];
        Max(c)=length(find(Yout==labels(c)));
        for i = 1:length(Users)
        
        IndC{c}=find(Yte{i}==labels(c))
        RND{c}=randperm(length(IndC{c}))
        IndC{c} = IndC{c}(RND{c})
        disp('... DONE')
        Xout{c}=cat(1,Xout{c},Xte{i}(:,:,IndC{c}(1:Max(c))))
    end
    
end