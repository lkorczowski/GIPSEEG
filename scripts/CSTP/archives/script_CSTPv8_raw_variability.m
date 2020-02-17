%% LOAD ERWAN DATA

%%
clc
close all
clear all
Directory= 'D:\data\erwan\'; %change path if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'FullErwanData_RAW.mat'])
load('BImapping')
%%
%%(0) #################### prepare inputs ####################
% extract EEG signals and prepare the parameters
clc
close all
clear Flashc   FlashcInd FlashCorrected Edec
allSubjects=[1:24]
        Nbiter=100
        nbTargets=[40,20,10,5];
        nbNonTarget=[300 30];
for Subjects=allSubjects; %choose the subject to load
    try
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
        winElec=[7,9,10,11,12,13,14,15,16]; %electrodes used for latency calculation and Pz
        winTime=[floor((0.05)*128):ceil((0.550)*128)]; 
        %(14)* time window (in sample) used for latency calculation and Pz selection
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
        Save={};% all the parameters
        % (1) converge criteria
        % (2) Pz
        % (3) Zbarz{iter}
        % (4) CSTP_all{iter}
        % (5) Conv{iter} the latency correction convergence criteria
        

        tic
        for i=1:Nbiter
            RND=randperm(length(Y));
            maxRND=find(cumsum(Y(RND))==5,1);
            while maxRND<25
                RND=randperm(length(Y));
                maxRND=find(cumsum(Y(RND))==5,1);
            end
            for h=1:length(nbNonTarget);
                for j=1:length(nbTargets)
                    Save(:,i,j,h)=CSTP_chain_variability(E,Flash,Y,Window,Delays,winElec,winTime,RND,nbTargets(j),nbNonTarget(h))';
                    PercentShow(( find(Subjects==allSubjects)-1)*Nbiter*length(nbTargets)*length(nbNonTarget)+ (i-1)*length(nbNonTarget)*length(nbTargets)+(h-1)*length(nbTargets)+j, length(nbTargets)*length(nbNonTarget)*Nbiter*length(allSubjects), 'Compute bootstrap ... ', [], true);
                end
            end
        end
        save(['D:\data\erwan\2Results_variability_CSTP_s' num2str(Subjects) '.mat'],'Save')
        
        toc
        %{
        if any(Subjects==[1 5 8 12 23])
            %disp('saving...')
            save(['Results_variability_CSTP_s' num2str(Subjects) '.mat'],'Save')
            
        end
        %}
    catch e
        disp(['ERROR subject' num2str(Subjects)])
        continue
    end
    %%%
end
toc
%%
clc
Weights=Save{6};
LatencyCorrection=Save{7};
iter=length(Save{7});
%norm(P300{6}(:,:,2)-P300{7}(:,:,2),'fro')
figure;subplot(411);plot(Save{1});title('Mean Zbar update');ylabel('µV');xlabel('iteration')
subplot(412);plot(Save{2});xlabel('iteration');ylabel('Pz');title('Nb Eigenvectors')
%norm(P300{5}(:,:,2)-P300{6}(:,:,2),'fro')
%subplot(513); hist(LatencyCorrection(Y==1,2),[-Delays:Delays]); title('Latencies TA');xlabel('delay (samples)');ylabel('nb');axis([-Delays +Delays 0 80])
subplot(413); hist(LatencyCorrection{iter}(Save{8}==1),[-Delays:Delays]); title('Latencies TA');xlabel('delay (samples)');ylabel('nb');axis([-Delays +Delays 0 80])
subplot(414); bar(log(Weights{iter}(Save{8}==1)));ylabel('log(W)');title('Weights TA');xlabel('Sweep (TA)')
text(40,-3,['Subject' num2str(Subjects)] ,'FontSize',20)

%Weights(find(Y==1,2))

%%plot P300 results
%close all
plotclass=2;
Scale=6;
fs=128;
CLASStxt={'NT','TA'};
Zbarz=Save{3};
figure
X=epoch_p300(E,Flash,Window); %get epochs
P300{1}=EnsembleAverage(X(:,:,Save{9}),Y(Save{9}));
[Bs Bt As At Class eigV]=CSTPover(X(:,:,Save{9}),Y(Save{9}),P300{1});
[Weights Xhatk]=WeightsEstimation(X(:,:,Save{9}),Y(Save{9}),Bs,Bt,As,At);
P300{3}=EnsembleAverage(Xhatk,Y(Save{9}));
max(eigV)
subplot(151);plotEEG(P300{1}(:,:,plotclass),Scale,fs,LABELS);title('X^{bar}(0.8)')
subplot(153);plotEEG(P300{3}(:,:,plotclass),Scale,fs,LABELS);title('X^{barhat}_{noweight}');%set(gca,'YTickLabel',[])
subplot(154);plotEEG(Zbarz{1}(:,:,plotclass),Scale,fs,LABELS);title('X^{bar}_{\sigma}');%set(gca,'YTickLabel',[])
subplot(155);plotEEG(Zbarz{end}(:,:,plotclass),Scale,fs,LABELS);title(['X^{bar}_{\sigma,\epsilon} iter' num2str(iter)]);%set(gca,'YTickLabel',[])
text(-2.5,Scale*3,['Subject' num2str(Subjects) ' ' CLASStxt{plotclass}] ,'FontSize',20)
%spaceplots(gca,[0 0 0 0],[1 1])
%%
plotclass=2;
Scale=4;
nbPlot=length(Zbarz)
fs=128;
CLASStxt={'NT','TA'};
figure
for plotInd=1:nbPlot
    subplot(1,nbPlot,plotInd);plotEEG(Zbarz{plotInd}(:,:,plotclass),Scale,fs,LABELS);title(['X^{bar}_{\sigma,\epsilon} iter' num2str(plotInd)]);
end
text(-2.5,-130,['Subject' num2str(Subjects) ' ' CLASStxt{plotclass}] ,'FontSize',20)
%spaceplots(gca,[0 0 0 0],[1 1])
%% FFT plot
figure
maxE=9:12
nbPlot=1
for elect=maxE
    subplot(length(maxE),1,nbPlot)
    freq=0:fs/fs:fs-1/fs;
    plot(freq,(abs(fft(P300{1}(elect,:,plotclass)))));hold all;plot(freq,(abs(fft(P300{3}(elect,:,plotclass)))));plot(freq,(abs(fft(P300{6}(elect,:,plotclass)))))
    legend('EA','CSTP no W no L','CSTP W L')
    xlabel('Freq (Hz)')
    xlim([0 fs/2]);title(['Electrode' num2str(elect)])
    nbPlot=nbPlot+1;
end
%% modify figures
close all
 Subjects=[2]
openfig(['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\variability\Subject' num2str(Subjects) '_scale.fig'],'new','visible')

hold on
set(gca,'FontSize',20)

set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','XLabel'),[])
colormap(gray)
spaceplots