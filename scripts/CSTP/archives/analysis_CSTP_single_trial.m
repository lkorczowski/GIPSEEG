clc
close all
clear all
Directory= 'D:\data\erwan\'; %change path if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'FullErwanData_RAW.mat'])
load('BImapping')
%%
Subjects=2
%%prepare data
        Session=1;
        Phase=find(cellfun(@(x) strcmp(x,'non-adaptive'), data(Subjects).session{Session}.phasetype));
        EEG=data(Subjects).session{Session}.phase{Phase}.training;
        
        Artefacts=[];
        %%%%%%%%%%%%%%%%%% function INPUTS %%%%%%%%%%%%%%%%%%%%%%
        decimationfactor=4; %put 1 to do nothing
        NOTCH=1; % put 1 to remove 50Hz
        BANDPASS=1; %put 1 to BANDPASS (see the filter bellow for the cutoff freq)
        
        E=EEG.s'; %(1) full EEG data
        Fs=EEG.Fs; %(2) sample frequency
        Flash=EEG.Flash;  %(3) Sweeps' indices
        Y=EEG.Y; %(4) Class of each sweep
        Window=1*Fs; % (5) the final sweep window will be 1s
        Delays=4; % (6) +/- nb shifted samples allowed for the jitter correction
        AddArte=0; % put 1 to simulate artefacts
        f1=1; %lowfreq cut (bandpass)
        f2=20; %highfreq cut (bandpass)
        N=4; %filter order (bandpass)
        %user parameters for the CSTP to improve converge :
        winElec=[7,10,11,12]; %electrodes used for latency calculation and Pz
        winTime=[6:96]; %time window used for latency calculation and Pz
        %winElec=[];
        %winTime=[];
        Ts=size(E,2);
        time1=40; %in second
        time2=180; %in second
        electrod=1;

        [Edec Fsd Flashd]=preprocessingEEG(E',Fs,[f1 f2 N decimationfactor],Flash);
        Fs=Fsd;
        Window=Window/decimationfactor;
        WindowB=Window+2*Delays;
        offset=-Delays;
        
        [Xall]=epoch_p300(Edec',Flashd,WindowB,offset); %prepare epoch for full window+latency
        [Flash Y]= removeEpochArtifact(Xall,Flashd,Y);
        clear E;
        E=Edec';
        P=EnsembleAverage(Xall,Y);
        %SNRmask(P(:,:,2),winElec,winTime)
        %%CSTP CHAIN
        %%%%%%%%%%%%%%%%%% prepare DATA accordingly %%%%%%%%%%%%%%%%%%%%%%
        
        
        %close all
        % Prepare the variable before the first iteration of the loop
        clear Zbarz Conv FlashCorrected P LatencyCorrection Weights Xhatk As At Bs Bt iterLat cond Bs As Bt At X
        %%
        %load(['D:\data\erwan\Results_variability_CSTP_s' num2str(Subjects) '.mat'])
        Trial=1;
        Noise=1;
        Epochs=1;
        % SAVE [ infos (see below) , nb trial , nb epochs, noise]
        % (1) converge criteria : Zbarz(iter)-Zbarz(iter-1)
% (2) Pz (dimension reduction of the CSTP)
% (3) Zbarz{iter}
% (4) CSTP_all{iter} Bs Bt As At
% (5) Conv{iter} the latency correction convergence criteria
% (6) Weights{iter} of each sweep
% (7) LantencyCorrection{iter} of each sweep
% (8) Y for the RND seed and nbTargets
% (9) RND seed
% (10) Ensemble Average
% (11) CSTP no weight


Pzfs=[size(Xall,1):-1:2];
iter=1;
maxIter=10;
CriteriaConvLatencies=1/length(find(Y))*Delays; %nb latencies correction allowed for convergence
for Pzf=Pzfs;
%[X]=epoch_p300(E',Flash,WindowB,offset,Artefacts);
LatencyCorrection(:,iter)=zeros(length(Y),1); %initialization latency at zero
%%(0) #################### Initialization Weights, CSTP ####################

X=epoch_p300(E,Flash,Window); %get epochs
Xnoise=X;
Ynoise=Y;
%%(3) #################### Weights Initialization ####################
Weights=WeightsEstimation(X,Y);

%%(3) #################### Average Initialization ####################
[P Class]=EnsembleAverage(X,Y);
Save{10}=P;
[Bs Bt As At Class eigV]=CSTPover(X,Y,P,[],winElec,winTime,Pzf,X,Y);
Save{11}{iter}=applyCSTP(P,Bs,Bt,As,At,Class);%
%%(3) #################### CSTP Initialization ####################
[Bs Bt As At Class eigV]=CSTPover(X,Y,P,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
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
    if iter>1
%while (Criteria>0.001 && iter<maxIter && (~(cond && Criteria<0.01))) || iter<2
    
    
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
        
        Conv(iter,iterLat)=ConvergenceLatencies(tempoLatency,LatencyCorrection(classLat,iter));
        
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
    X=epoch_p300(E,FlashCorrected,Window);
    [P Class]=EnsembleAverage(X,Y);
    %%(3) #################### CSTP Initialization ####################
    [Bs Bt As At Class eigV]=CSTPover(X,Y,P,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
    %%(3) #################### Zbarz Initialization ####################
    Zbarz{iter}=applyCSTP(P,Bs,Bt,As,At,Class);
    
    
    Criteria=ConvergenceZbarz(Zbarz{iter}(:,:,2),Zbarz{iter-1}(:,:,2));
    

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
            iter=iter+1;

%end
end
%% temproal info
close all
nbLINE=6
nbCOLUMN=4
Mapping={'Fp1';%1
    'Fp2';%2
    'F3';%3
    'AFz';%4
    'F4';%5
    'T7';%6
    'Cz';%7
    'T8';%8
    'P7';%9
    'P3';%10
    'Pz';%11
    'P4';%12
    'P8';%13
    'O1';%14
    'Oz';%15
    'O2'}%16
t = (0:(128-1))./128;
GFPmap=[1 128 0 80];
Pzind=4
[Bs]=Save{4}{Pzind}{1}
[Bt]=Save{4}{Pzind}{2}
[As]=Save{4}{Pzind}{3}
[At]=Save{4}{Pzind}{4}
Yw=Y
clims=[0 7]
ClassTag=1
Scale=4
figure(3)
subplot(nbLINE,4,[1 5 9]);plotEEG(mean(X(:,:,Y==1),3),Scale,128,Mapping) % Xbar TA

title('TA (A)')
subplot(nbLINE,4,[3 7 11]);plotEEG(mean(X(:,:,Y==0),3),Scale,128) %Xbar NT
title('NT (C)')
set(gca,'YtickLabel',[])

figure(1)
subplot(421)
GFP=global_field_power(X);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial')
subplot(423);plot(mean((GFP(:,Y==ClassTag)')))


subplot(422)
Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CSTP')
subplot(424);plot(mean((GFP(:,Y==ClassTag)')))

hold on
figure(3)

subplot(nbLINE,nbCOLUMN,[2 6 10]);plotEEG(mean(Xcstp(:,:,Y==1),3),Scale,128)
set(gca,'YtickLabel',[])
title('TA CSTP (B)')

subplot(nbLINE,nbCOLUMN,[4 8 12]);plotEEG(mean(Xcstp(:,:,Y==0),3),Scale,128)
set(gca,'YtickLabel',[])
title('NT CSTP (D)')

% global field power on average
%t=1:16
subplot(nbLINE,nbCOLUMN,13)
area(t,global_field_power(mean(X(:,:,Y==1),3)'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 6]);grid on;
ylabel('(1)','FontSize',16);set(gca,'YtickLabel',[])

subplot(nbLINE,nbCOLUMN,14)
area(t,global_field_power(mean(Xcstp(:,:,Y==1),3)'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) TA
ylim([0 6]);set(gca,'YtickLabel',[]);grid on

subplot(nbLINE,nbCOLUMN,15)
area(t,global_field_power(mean(X(:,:,Y==0),3)'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar NT
ylim([0 6]);set(gca,'YtickLabel',[]);grid on


subplot(nbLINE,nbCOLUMN,16)
area(t,global_field_power(mean(Xcstp(:,:,Y==0),3)'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) NT
ylim([0 6]);set(gca,'YtickLabel',[]);grid on


%single sweep global field power
subplot(nbLINE,nbCOLUMN,17)
imagesc(global_field_power(permute(X(:,:,Y==1),[2 1 3]))',clims); %mean GFP Xbar TA
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
ylabel('(2)','FontSize',16);
xticklabels = 0:.1:1;
xticks = linspace(1, 128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

subplot(nbLINE,nbCOLUMN,18)
imagesc(global_field_power(permute(Xcstp(:,:,Y==1),[2 1 3]))',clims); %mean GFP Xbar(cstp) TA
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
xticklabels = 0:.1:1;
xticks = linspace(1, 128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

subplot(nbLINE,nbCOLUMN,19)
imagesc(global_field_power(permute(X(:,:,Y==0),[2 1 3]))',clims); %mean GFP Xbar NT
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
xticklabels = 0:.1:1;
xticks = linspace(1, 128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

subplot(nbLINE,nbCOLUMN,20)
imagesc(global_field_power(permute(Xcstp(:,:,Y==0),[2 1 3]))',clims); %mean GFP Xbar(cstp) NT
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
xticklabels = 0:.1:1;
xticks = linspace(1, 128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

% average of single sweep global field power
%t=1:16
subplot(nbLINE,nbCOLUMN,21)
area(t,mean(global_field_power(permute(X(:,:,Y==1),[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 15]);grid on;
ylabel('(3)','FontSize',16);set(gca,'YtickLabel',[])
xlabel('Time (s)');

subplot(nbLINE,nbCOLUMN,22)
area(t,mean(global_field_power(permute(Xcstp(:,:,Y==1),[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) TA
ylim([0 15]);set(gca,'YtickLabel',[]);grid on
xlabel('Time (s)');

subplot(nbLINE,nbCOLUMN,23)
area(t,mean(global_field_power(permute(X(:,:,Y==0),[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar NT
ylim([0 15]);set(gca,'YtickLabel',[]);grid on
xlabel('Time (s)');


subplot(nbLINE,nbCOLUMN,24)
area(t,mean(global_field_power(permute(Xcstp(:,:,Y==0),[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) NT
ylim([0 15]);set(gca,'YtickLabel',[]);grid on
xlabel('Time (s)');

spaceplots

%% spatial info
close all
nbLINE=6
nbCOLUMN=4
Mapping={'Fp1';%1
    'Fp2';%2
    'F3';%3
    'AFz';%4
    'F4';%5
    'T7';%6
    'Cz';%7
    'T8';%8
    'P7';%9
    'P3';%10
    'Pz';%11
    'P4';%12
    'P8';%13
    'O1';%14
    'Oz';%15
    'O2'}%16
t = (0:(128-1))./128;
GFPmap=[1 80 1 16];
Pzind=4
[Bs]=Save{4}{Pzind}{1}
[Bt]=Save{4}{Pzind}{2}
[As]=Save{4}{Pzind}{3}
[At]=Save{4}{Pzind}{4}
Yw=Y
clims=[0 7]
ClassTag=1
Scale=4
figure(3)
subplot(nbLINE,4,[1 5 9]);plotEEG(mean(X(:,:,Y==1),3),Scale,128,Mapping) % Xbar TA
title('(A)')
%{
subplot(nbLINE,4,[3 7 11]);plotEEG(mean(X(:,:,Y==0),3),Scale,128) %Xbar NT
title('NT (A)')
set(gca,'YtickLabel',[])


figure(1)
subplot(421)
GFP=global_field_power(X);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial')
subplot(423);plot(mean((GFP(:,Y==ClassTag)')))


subplot(422)
Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CSTP')
subplot(424);plot(mean((GFP(:,Y==ClassTag)')))

hold on
figure(3)
%}
subplot(nbLINE,nbCOLUMN,[13 17 21]);plotEEG(mean(Xcstp(:,:,Y==1),3),Scale,128,Mapping)
%set(gca,'YtickLabel',[])
%title('TA (B)')
%{
subplot(nbLINE,nbCOLUMN,[4 8 12]);plotEEG(mean(Xcstp(:,:,Y==0),3),Scale,128)
set(gca,'YtickLabel',[])
title('NT (B)')
%}
% global field power on average
t=1:16
subplot(nbLINE,nbCOLUMN,[2 6 10])
area(t,global_field_power(mean(X(:,:,Y==1),3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 6]);grid on;
%ylabel('(1)','FontSize',16);set(gca,'YtickLabel',[])
%xticklabels = 1:1:16;
%xticks = linspace(1, 16, numel(xticklabels));
%set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
set(gca,'Fontsize',16)
title('(B)')
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);

subplot(nbLINE,nbCOLUMN,[14 18 22])
area(t,global_field_power(mean(Xcstp(:,:,Y==1),3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) TA
ylim([0 6]);set(gca,'YtickLabel',[]);grid on
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
set(gca,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);

%{
subplot(nbLINE,nbCOLUMN,15)
area(t,global_field_power(mean(X(:,:,Y==0),3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar NT
ylim([0 6]);set(gca,'YtickLabel',[]);grid on
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

subplot(nbLINE,nbCOLUMN,16)
area(t,global_field_power(mean(Xcstp(:,:,Y==0),3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) NT
ylim([0 6]);set(gca,'YtickLabel',[]);grid on
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
%}

%single sweep global field power
subplot(nbLINE,nbCOLUMN,[3 7 11])
imagesc(global_field_power((X(:,:,Y==1))),clims); %mean GFP Xbar TA
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
%ylabel('(2)','FontSize',16);
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);
set(gca,'Fontsize',16)

title('(C)')

subplot(nbLINE,nbCOLUMN,[15 19 23])
imagesc(global_field_power((Xcstp(:,:,Y==1))),clims); %mean GFP Xbar(cstp) TA
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);

%average single sweep global field power
subplot(nbLINE,nbCOLUMN,[4 8 12])
area(t,mean(global_field_power((X(:,:,Y==1))),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) TA
ylim([0 10]);set(gca,'YtickLabel',[]);grid on
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
set(gca,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);
title('(D)')

subplot(nbLINE,nbCOLUMN,[16 20 24])
area(t,mean(global_field_power((Xcstp(:,:,Y==1))),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) TA
ylim([0 10]);set(gca,'YtickLabel',[]);grid on
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
set(gca,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);
 %{
subplot(nbLINE,nbCOLUMN,19)
imagesc(global_field_power((X(:,:,Y==0)))',clims); %mean GFP Xbar NT
%set(gca,'XtickLabel',[]);
xlabel('# Electrod');axis(GFPmap);set(gca,'YtickLabel',[])
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

subplot(nbLINE,nbCOLUMN,20)
imagesc(global_field_power((Xcstp(:,:,Y==0)))',clims); %mean GFP Xbar(cstp) NT
set(gca,'XtickLabel',[]);xlabel('# Electrod');axis(GFPmap);set(gca,'YtickLabel',[])
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
%}
spaceplots

%% spatial info with topoplot
close all
nbLINE=6
nbCOLUMN=4
Mapping={'Fp1';%1
    'Fp2';%2
    'F3';%3
    'AFz';%4
    'F4';%5
    'T7';%6
    'Cz';%7
    'T8';%8
    'P7';%9
    'P3';%10
    'Pz';%11
    'P4';%12
    'P8';%13
    'O1';%14
    'Oz';%15
    'O2'}%16
t = (0:(128-1))./128;
GFPmap=[1 80 1 16];
Pzind=4
[Bs]=Save{4}{Pzind}{1}
[Bt]=Save{4}{Pzind}{2}
[As]=Save{4}{Pzind}{3}
[At]=Save{4}{Pzind}{4}
Yw=Y
clims=[0 7]
ClassTag=1
Scale=4
figure(3)
subplot(nbLINE,4,[1 5 9]);plotEEG(mean(X(:,:,Y==1),3),Scale,128,Mapping) % Xbar TA
title('(A)')
%{
subplot(nbLINE,4,[3 7 11]);plotEEG(mean(X(:,:,Y==0),3),Scale,128) %Xbar NT
title('NT (A)')
set(gca,'YtickLabel',[])


figure(1)
subplot(421)
GFP=global_field_power(X);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial')
subplot(423);plot(mean((GFP(:,Y==ClassTag)')))


subplot(422)
Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CSTP')
subplot(424);plot(mean((GFP(:,Y==ClassTag)')))

hold on
figure(3)
%}
subplot(nbLINE,nbCOLUMN,[13 17 21]);plotEEG(mean(Xcstp(:,:,Y==1),3),Scale,128,Mapping)
%set(gca,'YtickLabel',[])
%title('TA (B)')
%{
subplot(nbLINE,nbCOLUMN,[4 8 12]);plotEEG(mean(Xcstp(:,:,Y==0),3),Scale,128)
set(gca,'YtickLabel',[])
title('NT (B)')
%}
% global field power on average
t=1:16
subplot(nbLINE,nbCOLUMN,[2 6 10])
area(t,global_field_power(mean(X(:,:,Y==1),3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 6]);grid on;
%ylabel('(1)','FontSize',16);set(gca,'YtickLabel',[])
%xticklabels = 1:1:16;
%xticks = linspace(1, 16, numel(xticklabels));
%set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
set(gca,'Fontsize',16)
title('(B)')
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);

subplot(nbLINE,nbCOLUMN,[14 18 22])
area(t,global_field_power(mean(Xcstp(:,:,Y==1),3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar(cstp) TA
ylim([0 6]);set(gca,'YtickLabel',[]);grid on
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
set(gca,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
set(gca,'YtickLabel',[]);

subplot(nbLINE,nbCOLUMN,[3 4 7 8 11 12])
topoplot(global_field_power(mean(X(:,:,Y==1),3)),'Brain_Invaders_Erwan_16.locs')

subplot(nbLINE,nbCOLUMN,[15 16 19 20 23 24])
topoplot(global_field_power(mean(Xcstp(:,:,Y==1),3)),'Brain_Invaders_Erwan_16.locs')

%subplot(nbLINE,nbCOLUMN,[14 18 22])

spaceplots

%% mean + GFP temporal + GFP spatial + topoplot
for IND=1:4

close all
nbLINE=6
nbCOLUMN=4
Mapping={'Fp1';%1
    'Fp2';%2
    'F3';%3
    'AFz';%4
    'F4';%5
    'T7';%6
    'Cz';%7
    'T8';%8
    'P7';%9
    'P3';%10
    'Pz';%11
    'P4';%12
    'P8';%13
    'O1';%14
    'Oz';%15
    'O2'}%16
t = (0:(128-1))./128;
GFPmap=[1 128 0 80];
PZs=[16:-1:1];
Pzind=3
[Bs]=Save{4}{Pzind}{1}
[Bt]=Save{4}{Pzind}{2}
[As]=Save{4}{Pzind}{3}
[At]=Save{4}{Pzind}{4}
Yw=Y

Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
Xcsp=applyCSTP(X,Bs,{eye(128) eye(128)},As,{eye(128) eye(128)},Yw);
Xctp=applyCSTP(X,{eye(16) eye(16)},Bt,{eye(16) eye(16)},At,Yw);
GFP=global_field_power(Xcstp);
switch IND
    case 1
PLOTX=X(:,:,Y==1);
info='X';
    case 2
PLOTX=Xcstp(:,:,Y==1);
info='Xcstp';

    case 3
PLOTX=Xcsp(:,:,Y==1);
info='Xcsp';

    case 4
PLOTX=Xctp(:,:,Y==1);
info='Xctp';

end
xticklabels = 0:.2:1;

clims=[0 7]
ClassTag=1
Scale=5
hfig=figure(3)
set(hfig, 'PaperPosition', [0 0 15 30])
subplot(nbLINE,4,[1 2 5 6 9 10]);plotEEG(mean(PLOTX,3),Scale,128,Mapping) % Xbar TA
title('(A)')
%xlabel('Time (s)')
xticks = linspace(0, 1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

% global field power on average
%t=1:16
subplot(nbLINE,nbCOLUMN,[13 14])
area(t,global_field_power(mean(PLOTX,3)'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 6]);grid on;
ylabel('(1)','FontSize',16);set(gca,'YtickLabel',[])
xticks = linspace(0, 1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

%single sweep global field power
subplot(nbLINE,nbCOLUMN,[17 18])
[GFP1 GFP2]=global_field_power(permute(PLOTX,[2 1 3]));
imagesc(GFP1',clims); %mean GFP Xbar TA
colormap(gray)
%set(gca,'XtickLabel',[]);
axis(GFPmap);set(gca,'YtickLabel',[])
ylabel('(2)','FontSize',16);
xticks = linspace(1, 128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

% average of single sweep global field power
%t=1:16
subplot(nbLINE,nbCOLUMN,[21 22])
area(t,mean(global_field_power(permute(PLOTX,[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 15]);grid on;
ylabel('(3)','FontSize',16);set(gca,'YtickLabel',[])
xlabel('Time (s)');
xticks = linspace(0, 1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)    

% global field power on average
t=1:16
subplot(nbLINE,nbCOLUMN,[3 7 11])
area(t,global_field_power(mean(PLOTX,3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 3]);grid on;
set(gca,'Fontsize',16)
title('(B)')
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',16)
view(90,90)
set(gca,'XTickLabel',[])
text(0,3.5,['Subject ' num2str(Subjects)], 'FontSize',14)
text(1,3.5,[info], 'FontSize',14)
if IND~=1
text(2,3.5,['Pz ' num2str(PZs(Pzind))], 'FontSize',14)
end


%set(gca,'YtickLabel',[]);
%ylabel('µ')
subplot(nbLINE,nbCOLUMN,[15 16 19 20 23 24])
topoplot(global_field_power(mean(PLOTX,3))-min(global_field_power(mean(PLOTX,3))),'Brain_Invaders_Erwan_16.locs')
colormap(flipud(colormap))


spaceplots
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\CSTPvsCSPvsCTP\Subject_' num2str(Subjects) '_' info '.tiff'])
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\CSTPvsCSPvsCTP\Subject_' num2str(Subjects) '_' info '.fig'])

end
%%
figure(1)
hold off
subplot(425)
Xcstp=applyCSTP(X,Bs,{eye(128) eye(128)},As,{eye(128) eye(128)},Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CSP')
subplot(427);plot(mean((GFP(:,Y==ClassTag)')))



subplot(426)
Xcstp=applyCSTP(X,{eye(16) eye(16)},Bt,{eye(16) eye(16)},At,Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CTP')
subplot(428);plot(mean((GFP(:,Y==ClassTag)')))




figure
ClassTag=0
subplot(221)
subplot(421)
GFP=global_field_power(X);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial')
subplot(423);plot(mean((GFP(:,Y==ClassTag)')))


subplot(422)
Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CSTP')
subplot(424);plot(mean((GFP(:,Y==ClassTag)')))



subplot(425)
Xcstp=applyCSTP(X,Bs,{eye(128) eye(128)},As,{eye(128) eye(128)},Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==1)'./2),clims)
title('GFP single trial CSP')
subplot(427);plot(mean((GFP(:,Y==ClassTag)')))



subplot(426)
Xcstp=applyCSTP(X,{eye(16) eye(16)},Bt,{eye(16) eye(16)},At,Yw);
GFP=global_field_power(Xcstp);
imagesc((GFP(:,Y==ClassTag)'),clims)
title('GFP single trial CTP')
subplot(428);plot(mean((GFP(:,Y==ClassTag)')))
