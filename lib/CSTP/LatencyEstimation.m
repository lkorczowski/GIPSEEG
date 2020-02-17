function [Latency Crit_TraceN]=LatencyEstimation(Xall,X,Y,W,Bs,Bt,As,At,winElec,winTime)
% [Latency Crit_TraceN]=LatencyEstimation(Xall,X,Y,W,Bs,Bt,As,At,winElec,winTime)
% estimate recursively the latency for every sweep by finding the maximum
% covariance between all possible delays and the arithmic ensemble average filtered by CSTP.
% more information : see TIME SHIFT ESTIMATION p19
%
% WARNING : this method is quite time consuming. There are some possible
% improvement in order to boost the computation of all the covariance.
%
% [Latency Crit_TraceN]=LatencyEstimation(Xall,X,Y,W,Bs,Bt,As,At,electrodes,Time)
%
% see also : CorrectLatency
%
% *** History: 19-Mar-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%


Window=size(X,2);
SizeMax=size(Xall,2);
finalInd=SizeMax-Window;
delay0=floor(finalInd/2)+1;
if nargin<3
    Y=ones(size(X,3),1);
end
if nargin<4
    W=ones(size(Y));
end
Xw=applyWeights(X,W);
Xallw=applyWeights(Xall,W);
Class=unique(Y);

if nargin<5 || isempty(Bs)
    for i=1:length(Class)
        Bs{i}=eye(size(X,1));
        Bt{i}=eye(size(X,2));
    end
end
if nargin<7 || isempty(As)
    for i=1:length(Class)
        
        As{i}=eye(size(X,1));
        At{i}=eye(size(X,1));
    end
end
if nargin<9
    winElec=1:size(X,1);
end
if nargin<10
    winTime=1:size(X,2);
end

for k=1:size(Xall,3)
    
    z=find(Class==Y(k));
    Xother=Xw(:,:,[1:k-1 k+1:end]);
    Tag=Y([1:k-1 k+1:end]);
    Xother=Xother(:,:,Tag==Class(z));
    Pother=mean(Xother,3); %P300 estimation
    Ybar=Bs{z}'*Pother*Bt{z};
    %Xbar=As{z}*Ybar*At{z}';
    
    %%%%%%% WARNING : this part can to be optimize %%%%%%%%%%%%%%%%%%%%%%%
    for lat=1:finalInd+1
        Xtmp=Xallw(:,lat:lat+Window-1,k);
        Xwk=zeros(size(Xtmp));
        Xwk(winElec,winTime)=Xtmp(winElec,winTime);
        Yk=Bs{z}'*Xwk*Bt{z};
        %Xk=As{z}*Yk*At{z}';
        
        Crit_Trace(k,lat)=trace(Ybar*Yk'); % change method ?
        %Crit_Trace(k,lat)=trace(Xbar(electrodes,Time)*Xk(electrodes,Time)'); % change method ?
        
    end
    %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Crit_TraceN(k,:)=Crit_Trace(k,:)-min(Crit_Trace(k,:));
    [Pik Peaks]=findpeaks(Crit_TraceN(k,:)); %criteria
    
    if isempty(Peaks)
        Peaks=delay0;
    else
        indP=find(Pik>max(Pik)*0.66); %accept only peak of 66% of max (heuristic criteria)
        Peaks=Peaks(indP);
    end
    
    
    [trash IndM]=min(abs(Peaks-delay0)); %choose closest peak from allowed ones
    Latency(k)=Peaks(IndM)-delay0;
    
end

Latency=Latency';