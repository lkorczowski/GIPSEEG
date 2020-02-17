function Yestimated=estimate_binaryclass(scores,pTA)
% Yestimated=estimate_class(scores,pTA)
% works only for boolean classifier (for know)
%estimate class with a prior knowledge (such as probability), here the
%prior knowledge is the probability of occurence of target trials
Yestimated=zeros(size(scores));
nbTA=pTA*length(scores);
[~,ix]=sort(scores,'descend');
Yestimated(ix(1:nbTA))=1;