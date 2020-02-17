clear all
DirWave='D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\data\MARTI\TF\';
load('MARTI_ElectrodesName')
clear ToDoTest

elec=ft_read_sens('standard_1020.elc'); 
[IndicesElectrodes]=findElectrodes(elec.label, MARTImapping)%prepare elec structure

for Group=[(1:14), (16:21)]
for Player=1:2
for Condition=1

load([DirWave 'G' num2str(Group) '_p' num2str(Player) 'c' num2str(Condition) '.mat']);


NT.powspctrm((Group-1)*2+Player,:,:,:)=abs(double(TF.NT(:,3,:)));

TA.powspctrm((Group-1)*2+Player,:,:,:)=abs(double(TF.TA(:,3,:)));

end
end
end

NT.fsample=128;
NT.time=linspace(-1+1/128, 1, 256);
NT.dimord= 'rpt_chan_freq_time';
NT.label=MARTImapping';
NT.elec=elec;
NT.freq=3.5

TA.fsample=128;
TA.time=linspace(-1+1/128, 1, 256);
TA.dimord= 'rpt_chan_freq_time';
TA.label=MARTImapping';
TA.elec=elec;
TA.freq=3.5



% TF
%        label: {149x1 cell}
%        dimord: 'rpt_chan_freq_time'
%          freq: 20
%          time: [1x61 double]
%     powspctrm: [4-D double]
%     cumtapcnt: [77x1 double]
%           cfg: [1x1 struct]
%          grad: [1x1 struct] (or elec)
%
% allsubjFC{1} is
%       label: {152x1 cell}
%     fsample: 300
%         avg: [152x900 double]
%        time: [1x900 double]
%      dimord: 'chan_time'
%         cfg: [1x1 struct]
%        grad: [1x1 struct]
%type = ft_chantype(TA{1}.label)


load('ERF_orig');
allsubjFC{1}.cfg


cfg = [];
cfg.channel          = TA.label;
cfg.latency          = 'all';
cfg.frequency        = 3.5;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 3000;
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, TA);

design = zeros(1,size(TA.powspctrm,1) + size(NT.powspctrm,1));
design(1,1:size(TA.powspctrm,1)) = 1;
design(1,(size(TA.powspctrm,1)+1):(size(TA.powspctrm,1)+...
  size(NT.powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;
cfg.elec=TA.elec
[stat] = ft_freqstatistics(cfg, TA, NT);
save([DirWave 'stats.mat'],'stat')
%stat.raweffect = TA.powspctrm - NT.powspctrm;
%layout = ft_prepare_layout(cfg, TA)
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'raweffect';
cfg.zlim   = [-1e-27 1e-27];
cfg.layout = 'elec1020.lay';
ft_clusterplot(cfg, stat);

cfg=[]
cfg.xlim = [-0.4:0.2:1.4];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
TA.raweffect = TA.powspctrm - NT.powspctrm;
figure; ft_topoplotTFR(cfg,TA);
figure; ft_topoplotTFR(cfg,NT);
%%
%channel = ft_channelselection( 'all',{'Pz'})
cfg=[]
  cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 100;


cfg.design = design;
cfg.ivar = 1;

[stat] = ft_timelockstatistics(cfg, ToDoTest{1}, ToDoTest{2})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIELD TRIP TEST %%%%%%%%%%%%%%%%%%%
%% load individual subject data
load('ERF_orig');
 
% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_FC         = ft_timelockgrandaverage(cfg,allsubjFC{:});  
GA_FIC        = ft_timelockgrandaverage(cfg,allsubjFIC{:});

cfg = [];
cfg.showlabels  = 'yes';
cfg.layout    	= 'CTF151.lay';
figure; ft_multiplotER(cfg,GA_FC, GA_FIC)
 
cfg = [];
cfg.channel = 'MLT12';
figure; ft_singleplotER(cfg,GA_FC, GA_FIC)

time = [0.3 0.7];
% Scaling of the vertical axis for the plots below
ymax = 1.9e-13; 
figure; 
for isub = 1:10
    subplot(3,4,isub)
    % use the rectangle to indicate the time range used later 
    rectangle('Position',[time(1) 0 (time(2)-time(1)) ymax],'FaceColor',[0.7 0.7 0.7]);
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubjFC{isub}.time,allsubjFC{isub}.avg(52,:));
    plot(allsubjFIC{isub}.time,allsubjFIC{isub}.avg(52,:),'r');
    title(strcat('subject ',num2str(isub)))
    ylim([0 1.9e-13])
    xlim([-1 2])
end 
subplot(3,4,11); 
text(0.5,0.5,'FC','color','b') ;text(0.5,0.3,'FIC','color','r')
axis off

chan = 52;
time = [0.3 0.7];
 
% find the time points for the effect of interest in the grand average data
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));
 
% select the individual subject data from the time points and calculate the mean
for isub = 1:10
    values_FIC(isub)  = mean(allsubjFIC{isub}.avg(chan,timesel_FIC));   
    values_FC(isub)  = mean(allsubjFC{isub}.avg(chan,timesel_FC));
end
 
% plot to see the effect in each subject
M = [values_FC',values_FIC'];
figure; plot(M','o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        'subj7', 'subj8', 'subj9', 'subj10'}, 'location','EastOutside');
    
    FCminFIC = values_FC - values_FIC;
[h,p,ci,stats] = ttest(FCminFIC, 0, 0.05) % H0: mean = 0, alpha 0.05

% define the parameters for the statistical comparison
cfg = [];
cfg.channel     = 'MLT12';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,allsubjFC{:},allsubjFIC{:});   % don't forget the {:}!
% "{:}" means to use data from all elements of the variable
% Procedure
% In this tutorial we will consider a between-trials experiment, in which we analyze the data of a single subject. The statistical analysis for this experiment we perform both on axial and planar ERFs. The steps we perform are as follows:
% 
% Preprocessing and time-locked analysis with the ft_definetrial, ft_preprocessing and ft_timelockanalysis functions
% (Calculation of the planar gradient with the ft_megplanar and ft_combineplanar functions)
% Permutation test with the ft_timelockstatistics function
% Plotting the result with the ft_topoplotER function
%% test with fieldtrip data
% find the interesting segments of data
cfgFIC = [];                                           % empty configuration
cfgFIC.dataset                 = fullfile('D:', 'Mes Documents GIPSA', 'MATLAB', 'LK_TOOLBOX', 'data', 'Subject01.ds');       % name of CTF dataset  
cfgFIC.trialdef.eventtype      = 'backpanel trigger';
cfgFIC.trialdef.prestim        = 1;
cfgFIC.trialdef.poststim       = 2;
cfgFIC.trialdef.eventvalue     = 3;                    % trigger value for fully incongruent (FIC)
cfgFIC = ft_definetrial(cfgFIC);            

% remove the trials that have artifacts from the trl
cfgFIC.trl([15, 36, 39, 42, 43, 49, 50, 81, 82, 84],:) = []; 

% preprocess the data
cfgFIC.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfgFIC.demean     = 'yes';
cfgFIC.baselinewindow  = [-0.2 0];
cfgFIC.lpfilter   = 'yes';                              % apply lowpass filter
cfgFIC.lpfreq     = 35;                                 % lowpass at 35 Hz.

dataFIC_LP = ft_preprocessing(cfg);  

cfgFC = [];                                           % empty configuration
cfgFC.dataset                 = fullfile('D:', 'Mes Documents GIPSA', 'MATLAB', 'LK_TOOLBOX', 'data', 'Subject01.ds');       % name of CTF dataset  
cfgFC.trialdef.eventtype      = 'backpanel trigger';
cfgFC.trialdef.prestim        = 1;
cfgFC.trialdef.poststim       = 2;
cfgFC.trialdef.eventvalue     = 9;                    % trigger value for fully incongruent (FC)
cfgFC = ft_definetrial(cfg);            

% remove the trials that have artifacts from the trl
%cfgFC.trl([2, 3, 4, 30, 39, 40, 41, 45, 46, 47, 51, 53, 59, 77, 85],:) = []; 

% preprocess the data
cfgFC.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfgFC.demean     = 'yes';
cfgFC.baselinewindow  = [-0.2 0];
cfgFC.lpfilter   = 'yes';                              % apply lowpass filter
cfgFC.lpfreq     = 35;                                 % lowpass at 35 Hz.

dataFC_LP = ft_preprocessing(cfg); 
%%
cfg = [];
cfg.keeptrials = 'yes';
timelockFIC = ft_timelockanalysis(cfg, dataFIC_LP);
timelockFC  = ft_timelockanalysis(cfg, dataFC_LP);
%%

%load ERF_orig;
cfg = [];
cfg.planarmethod   = 'sincos';
  cfg.method='distance'
  cfg.neighbours = ft_prepare_neighbours(cfg, dataFIC_LP);
  timelockFIC_planar = ft_megplanar(cfg, timelockFIC);
timelockFC_planar  = ft_megplanar(cfg, timelockFC);
       
timelockFIC_planar_cmb = ft_combineplanar(cfg, timelockFIC_planar);
timelockFC_planar_cmb  = ft_combineplanar(cfg, timelockFC_planar);
  
timelockFIC_planar_cmb.grad = timelockFIC.grad;  % add the gradiometer structure
timelockFC_planar_cmb.grad  = timelockFC.grad;
     
%%
%cfg = [];
cfg.channel = {'MEG'};
cfg.latency = [0 1];
  cfg.feedback = 'yes';

    %ft_neighbourplot(cfg, allsubjFC{1});

  cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 100;

design = zeros(1,size(timelockFIC_planar_cmb.trial,1) + size(timelockFC_planar_cmb.trial,1));
design(1,1:size(timelockFIC_planar_cmb.trial,1)) = 1;
design(1,(size(timelockFIC_planar_cmb.trial,1)+1):(size(timelockFIC_planar_cmb.trial,1) + size(timelockFC_planar_cmb.trial,1)))= 2;

cfg.design = design;
cfg.ivar = 1;

[stat] = ft_timelockstatistics(cfg, timelockFIC_planar_cmb, timelockFC_planar_cmb)

save stat_ERF_planar_FICvsFC stat

clear all
clc
load stat_ERF_planar_FICvsFC


%% compute the average 
cfg = [];
cfg.keeptrials = 'no';   % now only the average, not the single trials
avgFIC_planar = ft_timelockanalysis(cfg, timelockFIC_planar);
avgFC_planar  = ft_timelockanalysis(cfg, timelockFC_planar);
cfg = [];
avgFIC_planar_cmb = ft_combineplanar(cfg, avgFIC_planar);
avgFC_planar_cmb  = ft_combineplanar(cfg, avgFC_planar);

% subtract avgFC from avgFIC
cfg = [];
cfg.operation = 'subtract'
cfg.parameter = 'avg';
raweffectFICvsFC     = ft_math(cfg,avgFIC_planar_cmb,avgFC_planar_cmb);


%% Highlight the channels (stat)
figure;  
timestep = 0.05;		%(in seconds)
sampling_rate = dataFC_LP.fsample;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples

pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% Remember to do the same for negative clusters if you want them!

for k = 1:20;
     subplot(4,5,k);   
     cfg = [];
     cfg.xlim =[j(k) j(k+1)];
     cfg.zlim = [-1.0e-13 1.0e-13];   
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     cfg.highlight = 'on';
     cfg.highlightchannel = find(pos_int);
     cfg.comment = 'xlim';
     cfg.commentpos = 'title';
     cfg.layout = 'CTF151.lay';
     ft_topoplotER(cfg, raweffectFICvsFC);
end