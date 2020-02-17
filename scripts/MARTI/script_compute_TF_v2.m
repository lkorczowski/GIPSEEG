%% (0) PREPARE Load the data player, filter, preprocessing DONT FUCKING DO IT (12h)
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\'
load([Directory 'MARTI_GroupsName.mat'])
LABELS=importdata([Directory 'MARTImapping.txt']);
DISP=0;startDISP=1700;endDISP=2000;

for indGroup=1:length(GroupsName)
    disp(['group' num2str(indGroup)])
    try
        load([Directory '\mat\dataPlayers\' GroupsName{indGroup} 'dataPlayers.mat'])
        for indPlayer=1:2
            EEG=dataPlayers(indPlayer);
            EEG.ElectrodesName=LABELS;
            offset=-EEG.Fs*1;
            
            
            if DISP, figure;
                
                Freqs=freqs_val(EEG.Fs,size(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:),1)/2);
                eegf=abs(fft(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:))).^2;
                subplot(211)
                plot(Freqs,log(mean(eegf(1:end/2,:),2)))
            end
            EEG=preprocessingEEG(EEG,[1 20 4 8])
            
            if DISP
                eegf=abs(fft(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:))).^2;
                subplot(212)
                plot(freqs_val(EEG.Fs,size(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:),1)/2),log(mean(eegf(1:end/2,:),2)))
            end
            
            %%
            prestim=1; % 3 seconds
            poststim=2;
            Overlap=1;
            EEGFT=EEG_LK2FT(EEG, prestim,poststim,Overlap);
            EEGFT.Conditions=EEG.Conditions;
            save([Directory '\mat\dataPlayers\' GroupsName{indGroup} 'AverageSeries_' num2str(indPlayer) '.mat'],'EEGFT')
            disp('file written')
        end
    catch e
        disp('error')
    end
end

%% (1) TF+COH COMPUTE the average and compute the TF (20min)
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\'
load([Directory 'MARTI_GroupsName.mat'])
LABELS=importdata([Directory 'MARTImapping.txt']); % or MARTI_ElectrodesName.mat
DISP=0;startDISP=1700;endDISP=2000;
close all
Conditions(1,:)=[107 109] ; %each conditions constrat (107,108,109,110)
% OR MULTI CONDITION MIXTURE
 Conditions=[107 109 ; 108 110 ; 107 108 ; 109 110]'; %each conditions constrats between SYNS/COCM
Methods=[1]% 2 4];
for indMethod=6
    %indMethod=1 (7 cycles) 4 (4 cycles) (5 (3 cycles)) (6(6 cycles,
    %collapsing conditions)
    %contrast SYNS collapsing COCM)
    DirectoryFig=['D:\data\Hyperscanning\MARTI\Groups\mat\COH' num2str(indMethod) '\']
    % groups to exclude 4 8 10 15
    for indGroup=1:length(GroupsName)
        clear COH MAmp Pcon UNI
        try
            for indPlayer=1:2
                load([Directory '\mat\dataPlayers\meanoverlap\' GroupsName{indGroup} 'AverageSeries_' num2str(indPlayer) '.mat'],'EEGFT')
                
                EEGFT.time=EEGFT.time(1:length(EEGFT.trial));
                EEGFT.trialinfo=EEGFT.trialinfo(find(EEGFT.sampleinfo),:);
                EEGFT.sampleinfo=EEGFT.Conditions(1:6:end);
                
                cfg              = [];
                cfg.output       = 'fourier';
                cfg.channel      = 'EEG';
                cfg.method       = 'wavelet';
                cfg.taper        = 'hanning';
                cfg.keeptrials = 'yes';
                cfg.foi          = 1:0.25:20;                         % analysis 2 to 30 Hz in steps of 2 Hz
                cfg.t_ftimwin    = 6./cfg.foi;   % length of time window = 0.5 sec (short time fourier)
                % OR
                cfg.width = 6; %number of the period of the wavelet
                cfg.toi          = -1:0.01:2;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                TFRhann = ft_freqanalysis(cfg, EEGFT);
                phi(:,:,:,:,indPlayer)=angle(TFRhann.fourierspctrm);
                
                UNI.MAmp(:,:,:,indPlayer)=squeeze(mean(abs(TFRhann.fourierspctrm),1));
                UNI.Pcon(:,:,:,:,indPlayer)=squeeze(abs(mean(TFRhann.fourierspctrm,1)));
            end
            UNI.label=TFRhann.label;
            UNI.dimord='chan_freq_time';
            UNI.time=TFRhann.time;
            UNI.freq=TFRhann.freq;
            
            save([DirectoryFig GroupsName{indGroup} '_univariate.mat'],'UNI')
            
            %% all conditions together
            close all
            COH.label=TFRhann.label;
            COH.dimord='chan_freq_time';
            COH.time=TFRhann.time;
            COH.freq=TFRhann.freq;
            COH.sampleinfo=EEGFT.sampleinfo;
            COH.powspctrm=abs(squeeze(mean(exp(i*(phi(:,:,:,:,1)-phi(:,:,:,:,2))),1)));
            save([DirectoryFig GroupsName{indGroup} '_PCoh_cond' num2str(0) '.mat'],'COH')
            
            cfg = [];
            cfg.baseline     = 'no';%[-0.5 -0.1];
            cfg.baselinetype = 'absolute';
            cfg.zlim         = [0 0.2];
            cfg.showlabels   = 'yes';
            cfg.layout       = 'elec1010.lay';
            cfg.channel = {'P4','PO8','P8'};%, 'Oz', 'O2'}
            cfg.figurename=['g' num2str(indGroup*10)];
            cfg.xlim          = [-0.5 1];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
            
            % ft_multiplotTFR(cfg, TFRhann
            ft_singleplotTFR(cfg,COH)
            print(gcf, [DirectoryFig 'TF\g' num2str(indGroup*10)],'-dtiff','-r450')
            
            cfg = [];
            cfg.baseline     = 'no';%[-0.5 -0.1];
            cfg.baselinetype = 'absolute';
            cfg.zlim         = [0 0.2];
            cfg.showlabels   = 'yes';
            cfg.layout       = 'elec1010.lay';
            cfg.channel = {'all'};%, 'Oz', 'O2'}
            cfg.figurename=['g' num2str(indGroup*10)];
            cfg.xlim          = [0 0.35];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
            cfg.colorbar='yes';
            cfg.marker = 'labels';
            cfg.ylim=[6 8];
            
            % ft_multiplotTFR(cfg, TFRhann
            ft_topoplotER(cfg, COH)
            print(gcf, [DirectoryFig 'Topo\g' num2str(indGroup*10)],'-dtiff','-r450')
            
            %% each condition
            close all
            indCond=1
            
            for Cond=Conditions %compute the Pcoh for each condition(s)
                %% save condition(s) Pcoh
                COH.powspctrm=abs(squeeze(mean(exp(i*(phi(find(ismember(COH.sampleinfo, Cond)),:,:,:,1)-phi(find(ismember(COH.sampleinfo, Cond)),:,:,:,2))),1)));%Phase coherence
                COH.condition=Cond;
                save([DirectoryFig GroupsName{indGroup} '_PCoh_cond' num2str(indCond) '.mat'],'COH')
                %% plot and save Pcoh time_freq for given electrodes
                % fig1
                cfg = [];
                cfg.baseline     = 'no';%[-0.5 -0.1];
                cfg.baselinetype = 'absolute';
                cfg.zlim         = [0 0.25];
                cfg.showlabels   = 'yes';
                cfg.layout       = 'elec1010.lay';
                cfg.channel = {'P4'};%, 'Oz', 'O2'}
                cfg.figurename=['g' num2str(indGroup*10) '_cond' num2str(indCond) '_' cfg.channel{:}];
                cfg.xlim          = [-0.5 1];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                ft_singleplotTFR(cfg,COH)
                print(gcf, [DirectoryFig 'TF\' cfg.figurename],'-dtiff','-r450')
                % fig2
                cfg.channel = {'PO8'};%, 'Oz', 'O2'}
                cfg.figurename=['g' num2str(indGroup*10) '_cond' num2str(indCond) '_' cfg.channel{:}];
                ft_singleplotTFR(cfg,COH)
                print(gcf, [DirectoryFig 'TF\' cfg.figurename],'-dtiff','-r450')
                
                %% plot and save Pcoh topo for given time and freq
                %fig1
                cfg = [];
                cfg.baseline     = 'no';%[-0.5 -0.1];
                cfg.baselinetype = 'absolute';
                cfg.zlim         = [0 0.25];
                cfg.showlabels   = 'yes';
                cfg.layout       = 'elec1010.lay';
                cfg.channel = {'all'};%, 'Oz', 'O2'}
                cfg.xlim          = [0 0.35]; % time window
                cfg.ylim=[6 8]; %frequency
                cfg.colorbar='yes';
                cfg.marker = 'labels';
                cfg.figurename=['g' num2str(indGroup*10) '_cond' num2str(indCond) 'time_' num2str(round(cfg.xlim*1000)) 'freq_' num2str(round(cfg.ylim))];
                ft_topoplotER(cfg,COH)
                print(gcf, [DirectoryFig 'Topo\' cfg.figurename ],'-dtiff','-r450')
                %fig2
                cfg.xlim          = [0.5 0.8]; % time window
                cfg.ylim=[2 3]; %frequency
                cfg.figurename=['g' num2str(indGroup*10) '_cond' num2str(indCond) 'time_' num2str(round(cfg.xlim*1000)) 'freq_' num2str(round(cfg.ylim))];
                ft_topoplotER(cfg,COH)
                print(gcf, [DirectoryFig 'Topo\' cfg.figurename ],'-dtiff','-r450')
                %fig3
                cfg.xlim          = [0.0 0.5]; % time window
                cfg.ylim=[2 3]; %frequency
                cfg.figurename=['g' num2str(indGroup*10) '_cond' num2str(indCond) 'time_' num2str(round(cfg.xlim*1000)) 'freq_' num2str(round(cfg.ylim))];
                ft_topoplotER(cfg,COH)
                print(gcf, [DirectoryFig 'Topo\' cfg.figurename ],'-dtiff','-r450')
                
                %% go to next condition(s)
                indCond=indCond+1;
            end
        catch e
            disp('error')
        end
    end
end

%% (2a) TMAP COMPUTE compute and save the Phase coherence (instantaneous)

clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\'
load([Directory 'MARTI_GroupsName.mat'])
LABELS=importdata([Directory 'MARTImapping.txt']);
DISP=0;startDISP=1700;endDISP=2000;
close all
Conditions=107:110;
DirectoryFig=['D:\data\Hyperscanning\MARTI\Groups\mat\COH6\']
errorgroups=[];
for indGroup=1:length(GroupsName)
    
    for indCond=1:4;
        try
            load([DirectoryFig GroupsName{indGroup} '_PCoh_cond' num2str(indCond) '.mat'],'COH')
            if indGroup==1 & indCond==1
                COHall.label=COH.label;
                COHall.time=COH.time;
                COHall.freq=COH.freq;
            end
            COHall.powspctrm(:,:,:,indGroup,indCond)=COH.powspctrm;
            
        catch e
            COHall.powspctrm(:,:,:,indGroup,indCond)=NaN;
            disp(['error in group' num2str(indGroup)])
            errorgroups=[errorgroups indGroup]
        end
    end
end
COHall.subj=GroupsName;
COHall.dimord='chan_freq_time_subj_rpt'
COHall.sampleinfo=[1 2 3 4]';
save([DirectoryFig 'COHall.mat'],'COHall')
%% (2b) TMAP PLOT t-map between conditions
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\'
DirectoryFig=['D:\data\Hyperscanning\MARTI\Groups\mat\COH6\']
load([Directory 'MARTI_GroupsName.mat'])
LABELS=importdata([Directory 'MARTImapping.txt']);
load([DirectoryFig 'COHall.mat'],'COHall')
close all
cfg = [];
cfg.subj=1:20
cfg.subj([4 8 9 10 15 17])=[];
cfg.baseline     = 'no';%[-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.zlim         = [-4 4];
cfg.showlabels   = 'yes';
cfg.layout       = 'elec1010.lay';
cfg.channel = {'all'};%, 'Oz', 'O2'}
cfg.xlim          = [0 0.30];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.colorbar='yes';
cfg.marker = 'labels';
cfg.ylim=[6 8];
% ft_multiplotTFR(cfg, TFRhann);
%print(gcf, [DirectoryFig 'Topo\g' num2str(indGroup*10) '_cond' num2str(indCond)],'-dtiff','-r450')
tic
SizeM=[size(COHall.powspctrm,1),size(COHall.powspctrm,2),size(COHall.powspctrm,3),size(COHall.powspctrm,4)*2];
coh1=COHall.powspctrm(:,:,:,cfg.subj,[1]);
coh1=permute(coh1(:,:,:,:),[4,1,2,3]);
coh2=COHall.powspctrm(:,:,:,cfg.subj,[2]);
coh2=permute(coh2(:,:,:,:),[4,1,2,3]);
[h,p,ci,stats] = ttest(coh1,coh2);
toc

COHstats=COHall;
COHstats.powspctrm=squeeze(stats.tstat)
COHstats.dimord='chan_freq_time';
% COHstats.powspctrm(find(abs(COHstats.powspctrm)<4))=0
% ft_multiplotTFR(cfg, TFRhann);
cfg.ylim=[2 3];
cfg.xlim          = [0 0.80];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
ft_topoplotER(cfg,COHstats)
print(gcf, [DirectoryFig 'tmap\tmap_topo' num2str(cfg.ylim)],'-dtiff','-r450')
cfg.xlim          = [0 0.80];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.ylim=[3 5];
ft_topoplotER(cfg,COHstats)
print(gcf, [DirectoryFig 'tmap\tmap_topo' num2str(cfg.ylim)],'-dtiff','-r450')
cfg.xlim          = [0 0.50];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.ylim=[5 7];
ft_topoplotER(cfg,COHstats)
print(gcf, [DirectoryFig 'tmap\tmap_topo' num2str(cfg.ylim)],'-dtiff','-r450')
cfg.xlim          = [0 0.50];                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.ylim=[7 9];
ft_topoplotER(cfg,COHstats)
print(gcf, [DirectoryFig 'tmap\tmap_topo' num2str(cfg.ylim)],'-dtiff','-r450')

cfg.ylim=[9 11];
ft_topoplotER(cfg,COHstats)
print(gcf, [DirectoryFig 'tmap\tmap_topo' num2str(cfg.ylim)],'-dtiff','-r450')

%% (3a) do clustering test on COH (needs 2a)
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\'
DirectoryFig=['D:\data\Hyperscanning\MARTI\Groups\mat\COH6\']
load([Directory 'MARTI_GroupsName.mat'])
LABELS=importdata([Directory 'MARTImapping.txt']);
load([DirectoryFig 'COHall.mat'],'COHall')
diff(COHall.time)
%remove the wrong subjects

cfg.subj=1:20;
cfg.subj([4 8 9 10 15 17])=[];
COHall.powspctrm=COHall.powspctrm(:,:,:,cfg.subj,:);
% COHall.subj=COHall.subj(cfg.subj);
elec=ft_read_sens('standard_1005.elc');
COHall.elec=elec;
COHall=removefields(COHall,'subj');
COHcond{1}=COHall;COHcond{1}.powspctrm=COHall.powspctrm(:,:,:,:,1);COHcond{1}.dimord='chan_freq_time_subj';
COHcond{2}=COHall;COHcond{2}.powspctrm=COHall.powspctrm(:,:,:,:,2);COHcond{2}.dimord='chan_freq_time_subj';
COHcond{3}=COHall;COHcond{3}.powspctrm=COHall.powspctrm(:,:,:,:,3);COHcond{3}.dimord='chan_freq_time_subj';
COHcond{4}=COHall;COHcond{4}.powspctrm=COHall.powspctrm(:,:,:,:,4);COHcond{4}.dimord='chan_freq_time_subj';


for indF=10:16 %4:length(COHall.freq)%:length(COHall.freq)
    
    
    %CO2.elec=elec;
    cfg = [];
    cfg.channel          = COHall.label;
    cfg.latency          = 'all';
    cfg.frequency        = COHall.freq(indF);
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT'; % changed to depsamples, for paired t-test
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;  % should be 0.05
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025; % should be 0.025 for two-sided
    cfg.numrandomization = 10000;
    cfg_neighb.method    = 'triangulation';
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, COHall);
    % cfg.parameter        = 'avg'; %contains avg cond values
    
    subj=size(COHall.powspctrm,4);
    design = zeros(2,2*subj);
    for i = 1:subj
        design(1,i) = i;
    end
    for i = 1:subj
        design(1,subj+i) = i;
    end
    design(2,1:subj)        = 1;
    design(2,subj+1:2*subj) = 2;
    
    
    cfg.design           = design;
    cfg.uvar             = 1;
    cfg.ivar             = 2;
    
    
    %[stats12{indF}] = ft_freqstatistics(cfg, COHcond{1}, COHcond{2});
    [stats34{indF}] = ft_freqstatistics(cfg, COHcond{3}, COHcond{4});
end
save([DirectoryFig 'stat_f4-6.mat'],'stats34') %'stats12',

%% (3b) Show Clustering results (needs 3a)
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\'
DirectoryFig=['D:\data\Hyperscanning\MARTI\Groups\mat\COH6\']
load([Directory 'MARTI_GroupsName.mat'])
LABELS=importdata([Directory 'MARTImapping.txt']);
load([DirectoryFig 'stat_f4-6.mat'],'stats34')
stat=stats34
for indF=15%:2:29%4:3:58
    try
        measure='Pcoh';
        
        %close all
        cfg = [];
        cfg.layout = 'EEG1010.lay' %'elec1020.lay';
        cfg.alpha  = 0.01;
        cfg.parameter = 'stat';
        cfg.zlim   = [-4 4];
        cfg.downsample= 10;
        
        
       % cfg.saveaspng=[DirectoryFig 'cluster_f4-6\' measure '_34_alpha_' num2str(cfg.alpha*100)  '_f' num2str(round(stat{indF}.freq))   ];
        
        ft_clusterplot(cfg, stat{indF});
    catch
        indF
    end
end