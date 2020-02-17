function e=SaveEEGepoch2p(Directory,AUsers,SAVEPATH,Artefacts)
%% Save all data
if(nargin<4 || isempty(Artefacts))
    Artefacts=[];
end

if(nargin<3 || isempty(SAVEPATH))
    disp('Please set save path')
    return
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
        for se=3; %selection session for users data
            %session 1 : player 1
            %session 2 : player 2
            %session 3 : all players
            
            test = data.subjects{i}.session{se}.online;
            Fs = test.Fs;
            window = Fs; % 1s window
            [Xte{i} isBad{i}] = epoch_p300(test.s,test.Flash,window,0,Artefacts);
            Yte{i}=test.Y;
            name{i}=[AUsers{i} '_' test.name];
            
        end
        catch e
                %error=[error; AUsers{i}];
    end
    end
    
    ALLdata.Xte=Xte;
    ALLdata.isBad=isBad;
    ALLdata.Yte=Yte;
    ALLdata.name=name;
    save(SAVEPATH,'ALLdata')