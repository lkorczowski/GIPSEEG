function Save=CSTP_chain(E,Flash,Y,Window,Delays,winElec,winTime)
% E : EEG signal
% Flash : indice of the begining of sweep (1 when start)
% Window : final length of the epoch
% Delays : delay allowed for recalibration of the flash (jitter correction)
% winElect : chosen spatial window of jitter calibration and Pz selection in the (for intance, you
% can select only Cz for P300 detection)


% OUTPUT Save{i} with i such as :
% (1) converge criteria : Zbarz(iter)-Zbarz(iter-1)
% (2) Pz (dimension reduction of the CSTP)
% (3) Zbarz{iter}
% (4) CSTP_all{iter} Bs Bt As At
% (5) Conv{iter} the latency correction convergence criteria
% (6) Weights{iter} of each sweep
% (7) LantencyCorrection{iter} of each sweep

Save={};% all the parameters

WindowB=Window+2*Delays;
offset=-Delays;
[Xall]=epoch_p300(E,Flash,WindowB,offset); %prepare epoch for full window+latency
maxIter=10;
CriteriaConvLatencies=1/length(find(Y))*Delays; %nb latencies correction allowed for convergence
if nargin<6 || isempty(winElec)
winElec=1:size(E,1); 
winTime=1:Window;
disp('WARNING : full electrods and time window for latency')
end


%[X]=epoch_p300(E',Flash,WindowB,offset,Artefacts);
LatencyCorrection=zeros(length(Y),2); %initialization latency at zero
%%(0) #################### Initialization Weights, CSTP ####################
iter=1;

X=epoch_p300(E,Flash,Window); %get epochs
%%(3) #################### Weights Initialization ####################
Weights=WeightsEstimation(X,Y);

%%(3) #################### Average Initialization ####################
[P Class]=meanOverlap(E,Flash,Window,Y,Weights);
%%(3) #################### CSTP Initialization ####################
[Bs Bt As At Class eigV]=CSTPover(X,Y,P,Weights,winElec,winTime);
%%(3) #################### Zbarz Initialization ####################
Zbarz{iter}=applyCSTP(P,Bs,Bt,As,At,Class);


%%#################### START LOOP ####################
Criteria=ConvergenceZbarz(Zbarz{iter}(:,:,2),P(:,:,2)); %while Zbarz converge

    Save{1}(iter)=Criteria;
    Save{2}(iter)=max(eigV);
    Save{3}{iter}=Zbarz{iter};
    Save{4}{iter}={Bs Bt As At};
    Save{5}=0;
    Save{6}{iter}=Weights;
    Save{7}{iter}=LatencyCorrection(:,iter);
cond=0;

while (Criteria>0.001 && iter<maxIter && (~(cond && Criteria<0.01))) || iter<2
    iter=iter+1
    tic
    
    %%(3) #################### Weights ####################
    [Weights Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At);
    
    
    %%(4) #################### Latency ####################
    classLat=Y==1;%ones(size(Y)); %chose class to compute latencies
    tempoLatency=[];
    STOPlat=1;
    iterLat=1;
    LatencyCorrection(classLat,iter)=zeros(size(LatencyCorrection(classLat,iter-1))); %initialize latencies
    
    while STOPlat %converge criteria
        [tempoLatency Crit_Trace]=LatencyEstimation(Xall(:,:,classLat),X(:,:,classLat),Y(classLat),Weights(classLat),Bs(2),Bt(2),As(2),Bt(2),winElec,winTime);
        
        %[LatencyCorrection(:,iter+1) Crit_Trace]=LatencyEstimation(Xall,X,Y,Weights,Bs,Bt,As,Bt);
        
        Conv(iter,iterLat)=ConvergenceLatencies(tempoLatency,LatencyCorrection(classLat,iter))
        
        if (iterLat>1 && Conv(iter,iterLat)<CriteriaConvLatencies) || iterLat>=maxIter 
            STOPlat=0;
        end
        %if iterLat>1 && STOPlat
        %    if Conv(iter,iterLat)<Conv(iter,iterLat-1)
        %          STOPlat=0;
        %    end
        %end
        LatencyCorrection(Y==1,iter)=tempoLatency;
        FlashCorrected=CorrectLatency(Flash,LatencyCorrection(:,iter));
        X=epoch_p300(E,FlashCorrected,Window);
        iterLat=iterLat+1;
        
    end

    
    %%(3) #################### Average Initialization ####################
    FlashCorrected=CorrectLatency(Flash,LatencyCorrection(:,iter));
    [P Class]=meanOverlap(E,FlashCorrected,Window,Y,Weights);
    X=epoch_p300(E,FlashCorrected,Window);
    %%(3) #################### CSTP Initialization ####################
    [Bs Bt As At Class eigV]=CSTPover(X,Y,P,Weights,winElec,winTime);
    %%(3) #################### Zbarz Initialization ####################
    Zbarz{iter}=applyCSTP(P,Bs,Bt,As,At,Class);
    
    
    Criteria=ConvergenceZbarz(Zbarz{iter}(:,:,2),Zbarz{iter-1}(:,:,2))
    toc
    Save{1}(iter)=Criteria;
    Save{2}(iter)=max(eigV);
    Save{3}{iter}=Zbarz{iter};
    Save{4}{iter}={Bs Bt As At};
    Save{5}=Conv;
    Save{6}{iter}=Weights;
    Save{7}{iter}=LatencyCorrection(:,iter);
    if iter>1
        cond=(Save{1}(iter)-Save{1}(iter-1))>0;
    end
end

%{
[Weights Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At);
Xhatkw=applyWeights(Xhatk,Weights);
P300{5}=EnsembleAverage(Xhatk,Y);

P300{6}=EnsembleAverage(Xhatkw,Y);
max(eigV)

[P Class]=EnsembleAverage(E,Y,FlashCorrected,Window,Weights);
P300{7}=applyCSTP(P,Bs,Bt,As,At,[0,1]);

norm(P300{6}(:,:,2)-P300{7}(:,:,2),'fro')
figure;subplot(411);plot(Save{1});title('Mean Zbar update');ylabel('�V');xlabel('iteration')
subplot(412);plot(Save{2});xlabel('iteration');ylabel('Pz');title('Nb Eigenvectors')
norm(P300{5}(:,:,2)-P300{6}(:,:,2),'fro')
%subplot(513); hist(LatencyCorrection(Y==1,2),[-Delays:Delays]); title('Latencies TA');xlabel('delay (samples)');ylabel('nb');axis([-Delays +Delays 0 80])
subplot(413); hist(LatencyCorrection(Y==1,iter),[-Delays:Delays]); title('Latencies TA');xlabel('delay (samples)');ylabel('nb');axis([-Delays +Delays 0 80])
subplot(414); bar(log(Weights(Y==1)));ylabel('log(W)');title('Weights TA');xlabel('Sweep (TA)')
text(40,-3,['Subject' num2str(Subjects)] ,'FontSize',20)

Weights(find(Y==1,2))
%}