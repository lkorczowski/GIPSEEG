clc
close all
clear all
Directory= 'D:\data\erwan\'; %change patvh if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'FullErwanData_RAW.mat'])
load('BImapping')
% Save
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
% (11) CSTP no weightSave
%%
for Subjects=1:24
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
        winElec=[7,10,11,12,13,14,15,16]; %electrodes used for latency calculation and Pz
        winTime=[floor((0.05)*128):ceil((0.550)*128)]; %time window used for latency calculation and Pz (50ms to 550ms)
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
        %P=EnsembleAverage(Xall,Y);
        %[P Class]=meanOverlap(E,Flash,Window,Y,Weights);
        %SNRmask(P(:,:,2),winElec,winTime)
        %%CSTP CHAIN
        %%%%%%%%%%%%%%%%%% prepare DATA accordingly %%%%%%%%%%%%%%%%%%%%%%
        
        
        %close all
        % Prepare the variable before the first iteration of the loop
        clear Zbarz Conv FlashCorrected P LatencyCorrection Weights Xhatk As At Bs Bt iterLat cond Bs As Bt At X

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
% (12) EA with toeplitz, weights, latencies


Pzfs=[size(Xall,1):-1:2];
iter=1;
maxIter=10;
CriteriaConvLatencies=1/length(find(Y))*Delays; %nb latencies correction allowed for convergence
[P Class]=meanOverlap(E,Flash,Window,Y); 
Save{10}=P; %EA
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

[P Class]=meanOverlap(E,Flash,Window,Y,Weights);
%subplot(121);plotEEG(P(:,:,1),3);subplot(122);plotEEG(P2(:,:,1),3);
[Bs Bt As At Class eigV]=ACSTP(X,Y,P,[],winElec,winTime,Pzf,X,Y);
Save{11}{iter}=applyCSTP(P,Bs,Bt,As,At,Class);%
%%(3) #################### CSTP Initialization ####################
[Bs Bt As At Class eigV]=ACSTP(X,Y,P,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
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
    %[P Class]=EnsembleAverage(X,Y);
    [P Class]=meanOverlap(E,FlashCorrected,Window,Y,Weights);

    %%(3) #################### CSTP Initialization ####################
    [Bs Bt As At Class eigV]=ACSTP(X,Y,P,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
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
    Save{12}=P;
%% all plot

close all
Scale=4 %init scale

for IND=1:2
nbLINE=9
nbCOLUMN=5
FontSize=20;
Mapping={'Fp1';%1
    'Fp2';%2
    'F5';%3
    'AFz';%4
    'F6';%5
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
PZs=[16:-1:1];
winElec2=[7 11 14 15 16];
winTime2=5:60;
[Pzind SNR]=best_Pz(Save{3},winElec,winTime);
%plot(SNR)
[Bs]=Save{4}{Pzind}{1}
[Bt]=Save{4}{Pzind}{2}
[As]=Save{4}{Pzind}{3}
[At]=Save{4}{Pzind}{4}
%Weights=Save{6}{Pzind};
 %   Latcorr=Save{7}{Pzind};
Yw=Y

Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
switch IND
    case 1
%PLOTX=mean(X(:,:,Y==1),3);
PLOTX=Save{10}(:,:,end);
info=['ss' num2str(Subjects) ' AEA'];
Scale=FindScaling(PLOTX,winElec,winTime);
    case 2
%PLOTX=Xcstp(:,:,Y==1);
PLOTX=mean(Save{3}{Pzind}(:,:,end),3);
info=['ss' num2str(Subjects) ' ACSTP(' num2str(PZs(Pzind)) ')'];

end
xticklabels = 0:.5:1;

clims=[0 7]
ClassTag=1
hfig=figure(2)
set(hfig, 'PaperPosition', [0 0 12.5 30])
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
subplot(nbLINE,nbCOLUMN,[1 2 6 7 11 12 16 17 21 22 26 27 31 32 36 37]+2*(IND-1));plotEEG(PLOTX,Scale,128,Mapping) % Xbar TA
set(gca,'yticklabel',[])
if IND>1
set(gca, 'color', [0.95 .95 .95])
end
set(gcf, 'color', [1 1 1])

title(info,'FontSize',FontSize)
%xlabel('Time (s)')
xticks = linspace(0, 1.01, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
set(gcf, 'InvertHardCopy', 'off');
% global field power on average
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
subplot(nbLINE,nbCOLUMN,[21 22]+20+2*(IND-1))
area(t,global_field_power(PLOTX'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 6]);grid on;
if IND==1
%ylabel('GFP','FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
end
set(gca,'YtickLabel',[])
xticks = linspace(0, 1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
set(gca, 'color', [0.95 .95 .95])


spaceplots
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\EAvsCSTP\7Subject_' num2str(Subjects) '.tiff'])
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\EAvsCSTP\7Subject_' num2str(Subjects) '.fig'])

end
%% mean + GFP temporal + GFP spatial + topoplot
for IND=1:4
FontSize=30
close all
nbLINE=12
nbCOLUMN=5
Mapping={'Fp1';%1
    'Fp2';%2
    'F5';%3
    'AFz';%4
    'F6';%5
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
GFPmap=[1 128 1 80];
PZs=[16:-1:1];
Pzind=9;
Yw=Y;

% find best Pz for each filter
clear Zbarz1 Zbarz2 Zbarz3
for tmp=1:length(Save{4})
[Bs]=Save{4}{tmp}{1};
[Bt]=Save{4}{tmp}{2};
[As]=Save{4}{tmp}{3};
[At]=Save{4}{tmp}{4};
Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
Zbarz1{tmp}=mean(Xcstp(:,:,Yw==1),3);
Xcsp=applyCSTP(X,Bs,{eye(128) eye(128)},As,{eye(128) eye(128)},Yw);
Zbarz2{tmp}=mean(Xcsp(:,:,Yw==1),3);
Xctp=applyCSTP(X,{eye(16) eye(16)},Bt,{eye(16) eye(16)},At,Yw);
Zbarz3{tmp}=mean(Xctp(:,:,Yw==1),3);
end
Pz1=best_Pz(Zbarz1,winElec,winTime); %CSTP
Pz2=best_Pz(Zbarz2,winElec,winTime); %CSP
Pz3=best_Pz(Zbarz3,winElec,winTime); %CTP

[Bs1]=Save{4}{Pz1}{1};
[Bt1]=Save{4}{Pz1}{2};
[As1]=Save{4}{Pz1}{3};
[At1]=Save{4}{Pz1}{4};
Xcstp=applyCSTP(X,Bs1,Bt1,As1,At1,Yw);
[Bs2]=Save{4}{Pz2}{1};
[As2]=Save{4}{Pz2}{3};
Xcsp=applyCSTP(X,Bs2,{eye(128) eye(128)},As2,{eye(128) eye(128)},Yw);
[Bt3]=Save{4}{Pz3}{2};
[At3]=Save{4}{Pz3}{4};
Xctp=applyCSTP(X,{eye(16) eye(16)},Bt3,{eye(16) eye(16)},At3,Yw);



GFP=global_field_power(Xcstp);
switch IND
    case 1
PLOTX=X(:,:,Y==1);
info='AEA';
Scale=FindScaling(mean(PLOTX,3),winElec,winTime);
    case 2
PLOTX=Xcstp(:,:,Y==1);
info=['ACSTP(' num2str(PZs(Pz1)) ')'];
Scale=FindScaling(mean(PLOTX,3),winElec,winTime);

    case 3
PLOTX=Xcsp(:,:,Y==1);
info=['ACSP(' num2str(PZs(Pz2)) ')'];
Scale=FindScaling(mean(PLOTX,3),winElec,winTime);

    case 4
PLOTX=Xctp(:,:,Y==1);
info=['ACTP(' num2str(PZs(Pz3)) ')'];
Scale=FindScaling(mean(PLOTX,3),winElec,winTime);

end
xticklabels = 0:.5:1;

clims=[0 7]
ClassTag=1
%Scale=5
hfig=figure(3)
set(hfig, 'PaperPosition', [0 0 15 nbLINE*5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
subplot(nbLINE,nbCOLUMN,[1 2 3 6 7 8 11 12 13 16 17 18 21 22 23 26 27 28 31 32 33 36 37 38]);plotEEG(mean(PLOTX,3),Scale,128,Mapping,FontSize) % Xbar TA
title('(A)')
%xlabel('Time (s)')
xticks = linspace(0, 1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

text(0.33,6*Scale,['ss ' num2str(Subjects)], 'FontSize',FontSize,'fontname','times new roman','FontAngle','italic')
text(0.33,4.5*Scale,[info], 'FontSize',FontSize,'fontname','times new roman','FontAngle','italic')

% global field power on average
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
subplot(nbLINE,nbCOLUMN,[21 22 23]+20)
area(t,global_field_power(mean(PLOTX,3)'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 6]);grid on;
ylabel('(B)','FontSize',FontSize,'fontname','times new roman','FontAngle','italic');set(gca,'YtickLabel',[])
xticks = linspace(0, 1.01, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'FontSize',FontSize,'fontname','times new roman','FontAngle','italic')

%single sweep global field power
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP single trial %%%%%%%%%%%%%%%
subplot(nbLINE,nbCOLUMN,[26 27 28 31 32 33]+20)
[GFP1]=global_field_power(permute(PLOTX,[2 1 3]));
imagesc(GFP1',clims); %mean GFP Xbar TA
colormap(gray)
%set(gca,'XtickLabel',[]);
%axis([1 1 1 1]);
set(gca,'YtickLabel',[])
ylabel('(C)','FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
xticks = linspace(1, 129, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'FontSize',FontSize,'fontname','times new roman','FontAngle','italic')

% average of single sweep global field power
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT mean GFP %%%%%%%%%%%%%%%
subplot(nbLINE,nbCOLUMN,[36 37 38]+20)
area(t,mean(global_field_power(permute(PLOTX,[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 15]);grid on;
ylabel('(D)','FontSize',FontSize,'fontname','times new roman','FontAngle','italic','fontname','times new roman','FontAngle','italic');set(gca,'YtickLabel',[],'FontSize',FontSize)
xlabel('Time (s)','FontAngle','italic','fontname','times new roman');
xticks = linspace(0, 1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'FontSize',FontSize,'fontname','times new roman','FontAngle','italic')    

% global field power on average
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP spatial %%%%%%%%%%%%%%%
%{
t=1:16
subplot(nbLINE,nbCOLUMN,[3 8 13 18 23 28 33 38])
area(t,global_field_power(mean(PLOTX,3)),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
ylim([0 4]);grid on;
set(gca,'Fontsize',FontSize)
title('(B)','fontname','times new roman','FontAngle','italic')
xticklabels = 1:1:16;
xticks = linspace(1, 16, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'Fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
view(90,90)
set(gca,'XTickLabel',[],'FontSize',FontSize)

%}

%set(gca,'YtickLabel',[]);
%ylabel('µ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT topoplot %%%%%%%%%%%%%%%
%{
subplot(nbLINE,nbCOLUMN,[61 62 63 64 66 67 68 69 71 72 73 74 76 77 78 79])
GFP=global_field_power(mean(PLOTX,3));
%GFP=(GFP-min(GFP))/max(GFP-min(GFP))*3-1.5;

topoplot((GFP),'Brain_Invaders_Erwan_16.locs','plotrad',0.55)
caxis([min(GFP) max(GFP)])
title('(E)','Fontsize',FontSize,'fontname','times new roman','FontAngle','italic','position',[-0.67 0],'rotation',+90)

colorbar('southoutside','FontSize',FontSize*1,'FontAngle','italic','fontname','times new roman')
%}
colormap(flipud(colormap))

%spaceplots
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\CSTPvsCSPvsCTP\7Subject_' num2str(Subjects) '_' info '.tiff'])
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\CSTPvsCSPvsCTP\7Subject_' num2str(Subjects) '_' info '.fig'])

end
%%
close all
CompLABELS={'c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12','c13','c14','c15','c16'}'
for Pzind=[1]
htopo=figure

[Bs]=(Save{4}{Pzind}{3}{2});
Bs=normPOW(Bs');
[Bt]=Save{4}{Pzind}{4}{2};
Bt=normEEG(Bt');
subplot(4,5,[1 6 11 16]);plotEEG(Bt,4,128,CompLABELS);xlabel('Time (s)')
for i=1:size(Bs,1)
    Count=[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20];
    %Count=[2 4 6 8 11 13 15 17 20 22 24 26 29 31 33 35];
    subplot(4,5,[Count(i)])
topoplot(abs(Bs(i,:)),'Brain_Invaders_Erwan_16.locs','plotrad',0.55);title(CompLABELS{i},'FontSize',FontSize)
set(htopo, 'PaperPosition', [0 0 30 25])
caxis([0 1])
end
hb=colorbar
set(hb,'Units','normalized', 'position', [0.92 0.3 0.02 0.6],'FontSize',FontSize);
colormap(gray)
colormap(flipud(colormap))
text(-5.8,5.5,['ss ' num2str(Subjects)], 'FontSize',FontSize*1.5,'fontname','times new roman','FontAngle','italic')

%spaceplots
saveas(gca,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\CSTP_components\7Subject_' num2str(Subjects) '_components_Pz' num2str(PZs(Pzind)) '.tiff'])
saveas(gca,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\final\CSTP_components\7Subject_' num2str(Subjects) '_components_Pz' num2str(PZs(Pzind)) '.fig'])
end
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
