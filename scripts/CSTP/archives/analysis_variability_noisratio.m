%%
Directory= 'D:\data\erwan\'; %change path if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'FullErwanData_RAW.mat'])


%% ANALYSIS FINAL
clear P
for Subjects=1:24
    close
    load(['D:\data\erwan\2Results_variability_CSTP_s' num2str(Subjects) '.mat'])
    %% 40
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
        [X]=epoch_p300(E,Flash,Window); %prepare epoch for full window+latency
    %%
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
    ratios=[ 40 20 10 5]
    clear Zbarz_noweight Zbarz_weight Zbarz Zbarz2
    Fs=128
    
    time=0:1/Fs:Fs-1/Fs;
    
    close all
    nbconditions1=size(Save,3);
    Pzind=3
    
    electrode=1%[7 9 11 12 13 14 15 16]
    FontS=16
    class=2
    xticklabels = 0:.25:1;
    winElec=[7 9 10 11 12 13 14 15 16];
    winTime=5:80;
    Fig=figure('units','normalized','outerposition',[0.1 0.1 0.5 .5])
    set(Fig, 'PaperPosition', [0 0 20 15])
    xticks = linspace(0, 1.01, numel(xticklabels));
    
    for j=1:nbconditions1
        %h = text(-0.25, 0.5, 'row 2');
        %set(h, 'rotation', 90)
        Scale=25;%*(j/2);
        hold on
        
        for i=1:size(Save,2)
            RND=Save{9,i,j,1};
            Y=Save{8,i,j,1};
            Pcond=EnsembleAverage(X(:,:,RND),Y);
            P(:,:,i)=Pcond(:,:,class);
            %P(:,:,i)=Save{10,i,j,1}(:,:,class);
            %P_gfp=global_field_power(P);
            P_gfp=global_field_power(permute(P,[3 1 2]));
            P_gfp=permute(P_gfp,[3 2 1]);
            
            [Pzind]=best_Pz(Save{3,i,j,1},winElec,winTime);
            Zbarz(:,:,i)=Save{3,i,j,1}{Pzind}(:,:,class);
            %Zbarz_gfp=global_field_power(Zbarz);
            Zbarz_gfp=global_field_power(permute(Zbarz,[3 1 2]));
            Zbarz_gfp=permute(Zbarz_gfp,[3 2 1]);
            
            [Pzind]=best_Pz(Save{3,i,j,2},winElec,winTime);
            Zbarz2(:,:,i)=Save{3,i,j,2}{Pzind}(:,:,class);
            %Zbarz_gfp2=global_field_power(Zbarz2);
            Zbarz_gfp2=global_field_power(permute(Zbarz2,[3 1 2]));
            Zbarz_gfp2=permute(Zbarz_gfp2,[3 2 1]);
            
        end
        
        hold off
        meanZbarz=mean(Zbarz,3);
        varZbarz=var(Zbarz,[],3);
        hold on;
        subplot(nbconditions1,3,1+3*(j-1))
        plotEEGvariability(P_gfp,electrode,Fs);
        axis([0 1 -0 Scale])
        set(gca,'YTickLabel',[],'FontSize',FontS,'fontname','times new roman','FontAngle','italic', 'XTick', xticks, 'XTickLabel', xticklabels);
        %if j==1; title('X_{bar}'); end
        if j==1; title('AEA'); end
        if j==4; xlabel('time (s)');
        else
            set(gca,'XTickLabel',[]);
        end
        
        ylabel([num2str(ratios(j)) ' sweeps'])
        
        
        subplot(nbconditions1,3,2+3*(j-1))
        [h1 h2 h3]=plotEEGvariability(Zbarz_gfp,electrode,Fs);
        axis([0 1 -0 Scale])
        set(gca,'YTickLabel',[],'FontSize',FontS,'fontname','times new roman','FontAngle','italic', 'XTick', xticks, 'XTickLabel', xticklabels);
        
        %if j==1; title('X^{barhat}_{\sigma,\epsilon}, good noise estimation');
        if j==1; title('ACSTP(A)');
            
            text(-1,Scale*1.25,['ss' num2str(Subjects)] ,'FontSize',FontS,'fontname','times new roman','FontAngle','italic')
        end
        if j==4; xlabel('time (s)');
        else
            set(gca,'XTickLabel',[],'FontSize',FontS,'fontname','times new roman','FontAngle','italic');
        end
        
        subplot(nbconditions1,3,3+3*(j-1))
        [h1 h2 h3]=plotEEGvariability(Zbarz_gfp2,electrode,Fs);
        axis([0 1 -0 Scale])
        set(gca,'YTickLabel',[],'FontSize',FontS,'fontname','times new roman','FontAngle','italic', 'XTick', xticks, 'XTickLabel', xticklabels);
        
        %if j==1; title('X^{barhat}_{\sigma,\epsilon}, bad noise estimation'); end
        if j==1; title('ACSTP(B)'); end
        if j==4; xlabel('time (s)');
        else
            set(gca,'XTickLabel',[]);
        end
        if j==1; legend([h1 h2 h3],{'all','mean','90%'},'location',[0.62 0.87 0.1 0.1]); end
        
        colormap('gray')
    end
    set(gcf,'color','white')
    %set(gca,'fontsize',20);set(gca,'fontsize',20);
    %set(findobj(gcf, 'type','axes'), 'Visible','off')
    %set(gca,'XTickLabel',[]);
    spaceplots
    saveas(Fig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\variability\4Subject' num2str(Subjects)  '.fig'],'fig')
    saveas(Fig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\variability\4Subject_' num2str(Subjects) '.tiff'])
    
    
end