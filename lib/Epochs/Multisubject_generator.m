function [Xout labels RND]=Multisubject_generator(varargin)
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

% if mod(nargin,2)==1
%     disp('Multisubject_generator needs both EEG set and class list')
%     return
% end
if nargin>3
    Xte=varargin((1:2:nargin)); %extract Multiset
    Yte=varargin((2:2:nargin));
end
if nargin==3
    RND=varargin{3};
else
    RND={};
end
if nargin<4
    Xte=varargin{1}; %extract Multiset
    Yte=varargin{2};
end
lUsers=length(Xte);
%% generate multiplayers set
for i=1:length(Xte)
   MAX=size(Xte{i},3);
   Yte{i}=Yte{i}(1:MAX); % check the number of epochs
end
[trash IndF]=min(cellfun('length',Yte));

Yout=Yte{IndF};
labels = unique(Yout);
Nclass = length(labels);
IndC = cell(Nclass,1);%index of the class



Xout={};
if isempty(RND);
    RND=cell(1,Nclass);
    New=1;
else
    New=0;
end
for c=1:Nclass
    Xout{c}=[];
    Max(c)=length(find(Yout==labels(c)));
    for i = 1:lUsers
        
        IndC{c}=find(Yte{i}==labels(c));
        if New
            RND{c}=randperm(length(IndC{c}));
            
        end
        IndC{c} = IndC{c}(RND{c});
        
        Xout{c}=cat(1,Xout{c},Xte{i}(:,:,IndC{c}(1:Max(c))));
    end
    
end