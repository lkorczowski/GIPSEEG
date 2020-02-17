function [Yestimated, scores, errclass,irep,C,P1]=mdm_adaptive_hyper(X,Y,P1init,Cinit,counter,P)
%this function can be called recusively as long as X has at least one full repetition
%(ie 2 targets and 10 non-targets). Each repetition should be consecutive
%but the ordre of the target/non-target trials within are irrelevent.
labels=unique(Y);
Nclass=length(labels);
lrep=12;%trials in each repetition
pTA=1/6;%probability of target
if nargin<5 | isempty(counter)
    counter=0;
end

if nargin<6 | isempty(P)
    P = struct('nusers',1,'users',1);
end

if size(Y,2)>size(Y,1),Y=Y';end
%repetition
if (rem(length(Y)/lrep,1)~=0) | any(~sum(reshape(Y,lrep,length(Y)/lrep))==2)
    error('must be ordered by repetition (i.e. 2 targets, 10 non-targets');
end

nrep=length(Y)/lrep;
irep=0;
C=Cinit;
P1=P1init;
Yestimated=zeros(size(Y));
while irep<nrep
    irep=irep+1;
    indrep=(irep-1)*lrep+1:irep*lrep;%indices of the current repetition
    
    %generate features
    [COVte P1tmp]=covariances_p300_hyper(X(:,:,indrep),P1,length(P.users),P.P300_ref_orientation);
    
    % simulated online classification
    if 0 % only if P.classifieroptions.multiusers==multi
        COVte=covariances_p300_clean_inter(COVte,NbUsers,P.P300_ref_orientation);
    end
    d=[];
    for j=1:size(COVte,3)
        %classify the next trial
        for i=1:length(C)
            [d(j,i)] = distance(COVte(:,:,j),C{i},P.method_dist);
        end
    end
    scores(indrep)=-diff(d');
    Yestimated(indrep)=estimate_binaryclass(scores(indrep),pTA)';
    errclass(irep)=sum(abs(Yestimated(indrep)-Y(indrep)));
    % update the classifier with known tag (supervised) with the
    % forgetting factor alpha.
    %             alpha=1/(NTrainingtrial+j);
    alpha=1/(irep+counter);%alpha=0, no training, alpha=1, no memory
    for c=1:Nclass
        Clocal{c} = mean_covariances(COVte(:,:,Y(indrep)==labels(c)),P.method_mean);
        C{c} = geodesic(C{c},Clocal{c},alpha,P.method_mean);
    end
    P1=(1-alpha)*P1+alpha*P1tmp;

end

end