function SaveErwandata(AUsers,SAVEPATH)
%% Save all data

if(nargin<2 || isempty(SAVEPATH))
    SAVEPATH='ALLdata.mat';
    end
 nusers=length(AUsers);
 allcombi=nchoosek(AUsers,nusers);
se=1; %selection session for users data 
ph=2; %selection phase for users data

% P300 moyen et structure de bruit
    Xtr={};
    Xte={};
    Yte={};
    data={};
    %load data set
    
for i = 1:length(AUsers)
    % load data
    data.subjects{i} = load_erwan_data(AUsers{i});
    % Epoch signal
    training = data.subjects{i}.session{se}.phase{ph}.training;
    test = data.subjects{i}.session{se}.phase{ph}.online;
    Fs = test.Fs;
    window = Fs; % 1s window
    Xtr{i} = epoch_p300(training.s,training.Flash,window);
    Xte{i} = epoch_p300(test.s,test.Flash,window);
    Yte{i}=test.Y;
    Ytr{i}=training.Y;
    
end

ALLdata.Xtr=Xtr;
    ALLdata.Xte=Xte;
    ALLdata.Yte=Yte;
    ALLdata.Ytr=Ytr;
    save(SAVEPATH,'ALLdata')