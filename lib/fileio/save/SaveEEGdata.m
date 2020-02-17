function e=SaveEEGdata(Directory,AUsers,SAVEPATH, Artefacts)
%% Save all data
if(nargin<4 || isempty(Artefacts))
    Artefacts=[];
end

if(nargin<3 || isempty(SAVEPATH))
    SAVEPATH='ALLdata.mat';
    end
 nusers=length(AUsers);

% P300 moyen et structure de bruit
Xtr={};
    Xte={};
    Yte={};
    data={};
    %load data set
    se=1; %selection session for users data 

    
for i = 1:length(AUsers)
    % load data
    errors=[];
    try
    data.subjects{i} = load_EEG_data(Directory,AUsers{i});
    % Epoch signal
    se=1; %selection session for users data 

    test = data.subjects{i}.session{se}.online;
    Fs = test.Fs;
    window = Fs; % 1s window
    Xte{i} = epoch_p300(test.s,test.Flash,window,0,Artefacts);
    Yte{i}=test.Y;
    catch e
        %error=[error; AUsers{i}];
    end
end

    ALLdata.Xte=Xte;
    ALLdata.Yte=Yte;
    save(SAVEPATH,'ALLdata')