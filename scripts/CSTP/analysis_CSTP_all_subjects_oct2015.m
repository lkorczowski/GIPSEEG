clc
close all
clear all
Directory= 'D:\data\erwan\'; %change patvh if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'mat\FullErwanData_RAW.mat'])
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
%% compute and save the ACSTP
for Subjects=1:24
    %%prepare data
    Session=1;
    Phase=find(cellfun(@(x) strcmp(x,'non-adaptive'), data(Subjects).session{Session}.phasetype));
    EEG=convert_MAT_gdf2lk(data(Subjects).session{Session}.phase{Phase}.training);
    decimationfactor=4; %put 1 to do nothing
    f1=1; %lowfreq cut (bandpass)
    f2=20; %highfreq cut (bandpass)
    N=4; %filter order (bandpass)
    [EEG]=preprocessingEEG(EEG,[f1 f2 N decimationfactor]);
    
    % ACSTPoptions is a structure with
    
    %      Epoch_size: scalar, the length of the epoch window (in samples)
    % LatencyCorr_max: scalar, the maximum of samples for the latency
    %                   correction. Set 0 to disable the Latency correction.
    % Mask_Electrodes: vector of the selectionned electrodes. Usefull for the
    %                  automatic subspace selection (BestPz) and latency
    %                  correction (Latency).
    %       Mask_Time: vector of the selectionned sample. Usefull for the
    %                  automatic subspace selection (BestPz) and latency
    %                  correction (Latency).
    % MaxIterLatency*: scalar, the maximum iteration allowed to compute the
    %                   latency correction. Default: 10.
    %    SubspaceDim*: vector, containing all the subspace dimension (<nb electrodes)
    %                  to test in descent order.
    %                   By default, it is equal to (nb_channels:-1:(nb_channels/2))
    %computeClassLat*: vector, containing all the class tag in which you want
    %                   to compute the latency correction. By default, it
    %                   computes it for all classes but it could be long (for
    %                   instance you can skip the non-target).
    %        Weights*: Default: true.
    %                  option1(given) [nb epochs x1] vector, containing the weights for each
    %                   epoch if it is computed from an external function.
    %                  option2 (true/false) boolean. If true (default) the ACSTP compute the
    %                   weight for each epoch. If false, the weights are
    %                   desactivated (i.e. set to 1).
    %        DISPLAY*: Boolean, should the comparative result between the arithmetic ensemble
    %                   average and the ACSTP should be display at the end. Default: true.
    %
    %   *optional
    
    ACSTPoptions.Epoch_size=1*EEG.Fs;%(5) the final sweep window will be 1s
        ACSTPoptions.LatencyCorr_max=4; % (6) +/- nb shifted samples allowed for the jitter correction
        
        %user parameters for the CSTP to improve converge :
        ACSTPoptions.Mask_Electrodes=[7,10,11,12,13,14,15,16]; %electrodes used for latency calculation and Pz
        ACSTPoptions.Mask_Time=[floor((0.05)*EEG.Fs):ceil((0.550)*EEG.Fs)]; %time window used for latency calculation and Pz (50ms to 550ms)
        ACSTPoptions.computeClassLat=[1]; %just compute for TA
        Epochs(Subjects).X=epoch_p300(EEG.Channels',EEG.Trigger,EEG.Fs,0);
        Epochs(Subjects).EpochClass=EEG.EpochClass;
        [Epochs(Subjects).Xhat ACSTPstruct(Subjects)]=ACSTP(EEG,ACSTPoptions);
        


%     Save{1}(iter)=Criteria;
%     Save{2}(iter)=max(eigV);
%     Save{3}{iter}=Zbarz{iter};
%     Save{4}{iter}={Bs Bt As At};
%     Save{5}=Conv;
%     Save{6}{iter}=Weights;
%     Save{7}{iter}=LatencyCorrection(:,iter);
% 
%     Save{12}=P;
end

save([Directory 'results\ACSTPresults.mat'],'Epochs','ACSTPstruct','ACSTPoptions')
%% plot AEA versus ACSTP / topolot / components
clc
close all
clear all
Directory= 'D:\data\erwan\'; %change patvh if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'results\ACSTPresults.mat'])
load('BImapping')
%%plot AEA versus ACSTP
ClassesName={'NT','TA'}
for Subjects=1:24%length(ACSTPstruct)
    for indClass=1:2
    close all
    figure
    outputfig=['D:\Mes Documents GIPSA\MATLAB\figures\CSTP\ERWAN\EAvsCSTP\ss' num2str(Subjects) '_'];
    PrintOptions = struct('outputfolder',outputfig,'scale',5,'fontsize',20,'fs',128,'offset',0,'label',{LABELS},'paperfactor',1,'indclass',indClass,...
        'title',['ss' num2str(Subjects) ' ' ClassesName{indClass} ]);
    ACSTPprint(ACSTPstruct(Subjects),PrintOptions)
    end
end



%%mean + GFP temporal + GFP spatial + topoplot
for Subjects=1:length(ACSTPstruct)
  close all
    figure
    outputfig=['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ERWAN\CSTPvsCSPvsCTP\ss' num2str(Subjects) '_'];
    PrintOptions = struct('outputfolder',outputfig,'scale',3.5,'fontsize',16,'fs',128,'offset',0,'label',{LABELS},'paperfactor',1,'indclass',2,'title',['ss' num2str(Subjects)]);
ACSTPcspctp(Epochs(Subjects),ACSTPstruct(Subjects),ACSTPoptions,PrintOptions)

end
%%components
for Subjects=1:length(ACSTPstruct)
  close all
    figure
    outputfig=['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ERWAN\CSTP_components\ss' num2str(Subjects)];
    PrintOptions = struct('outputfolder',outputfig,'scale',3,'fontsize',16,'fs',128,'offset',0,'label',{LABELS},'paperfactor',1,'indclass',2,'title',['ss' num2str(Subjects)]);
ACSTPcomponents(ACSTPstruct(Subjects),PrintOptions)

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
