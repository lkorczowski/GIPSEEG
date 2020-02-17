EEG=data.session(3); %concatenate the conditions first

FILT=[1 20 4 4];
[EEG]=preprocessingEEG(EEG,FILT)
EEG.ElectrodesName={'Fp1';%1
                        'Fp2';%2
                        'AFz';%3
                        'F7';%4
                        'F3';%5
                        'F4';%6
                        'F8';%7
                        'FC5';%8
                        'FC1';%9
                        'FC2';%10
                        'FC6';%11
                        'T7';%12
                        'C3';%13
                        'Cz';%14
                        'C4';%15
                        'T8';%16
                        'Cp5';%17
                        'Cp1';%18
                        'Cp2';%19
                        'Cp6';%20
                        'P7';%21
                        'P3';%22
                        'Pz';%23
                        'P4';%24
                        'P8';%25
                        'PO7';%26
                        'O1';%27
                        'Oz';%28
                        'O2';%29
                        'PO8';%30
                        'PO9';%31
                        'PO10' };%32
%let's call X [nb_channels x Epoch_size x nb_Epochs] the epochs
for player=2
Artefacts=[];
[X isBad] = epoch_p300(EEG.Signal(:,1+(player-1)*32:(player)*32)',EEG.Trigger,EEG.Fs,0,Artefacts);
end

Conditions=ones(length(EEG.EpochClass),1);
TA= mean(X(:,:,EEG.EpochClass==1 & Conditions==1),3);
subplot(10,2,1:2:16)
plotEEG(TA,4,EEG.Fs,EEG.ElectrodesName)

nbpoints=10
% you choose
AvWindow=ones(nbpoints,1);%uniform average (moving average)
AvWindow=window(@hamming,nbpoints);%uniform average (moving average)

AvWindow=AvWindow/sum(AvWindow);

winTime=(0.250*EEG.Fs):(0.500*EEG.Fs);

Z = conv2(1,AvWindow,TA,'same')
subplot(10,2,2:2:16)
plotEEG(Z,4,EEG.Fs,EEG.ElectrodesName)

timeVector=(0:EEG.Fs-1)/EEG.Fs;
% here there different ways:

%% (1) find the maximum (maximum in time of the maximum in electrode)
[~,idx]=max(max(Z)); %if you want to extract WHICH electrode you need to get the indice of maximum amplitude along the electrodes first.
idx/EEG.Fs*1000 %in ms

%% (2) find the maximum on spatial average
Feat=mean(Z,1);
[~,idx]=max(Feat);
idx/EEG.Fs*1000 %in ms
subplot(10,2,[17 18])
plot(timeVector,Feat)

%% (3) find the maximum on the GFP
Feat=global_field_power(Z');
[~,idx]=max(Feat)
idx/EEG.Fs*1000 %in ms
subplot(10,2,[19 20])
plot(timeVector,Feat)

    [PKS,LOCS] = findpeaks(Feat,'THRESHOLD',0)%,'NPEAKS',2)
        max(PKS)
    LOCS(PKS>(max(PKS)*0.66))/EEG.Fs*1000 %in ms
%% let's call Xhat [nb_channels x Epoch_size x nb_Epochs] the ACSTP's epochs
[Xhat Struct]=ACSTP(fsdgfd)
% With the output ACSTPstruct, you can manually apply the filter on the EEG
% such as :
% TriggerCorrected=CorrectLatency(EEG.Trigger,ACSTPstruct.Latency); % apply latency correction
% X=epoch_p300(EEG.Signal,TriggerCorrected,ACSTPstruct.Epoch_size); % extract the epochs
% Xw=applyWeights(X,ACSTPstruct.Weights); %apply the weights
% Xhat=applyCSTP(Xw,ACSTPstruct.Bs,ACSTPstruct.Bt,ACSTPstruct.As,ACSTPstruct.At,EEG.EpochClass); %apply the ACSTP


