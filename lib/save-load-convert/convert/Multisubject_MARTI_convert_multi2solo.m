function dataPlayers=Multisubject_MARTI_convert_multi2solo(data,nb_subjects)
%% convert, concatenate the sessions
% INPUT data structure format for each session
%            Fs: scalar
%        StimPos: [nb_epochs x1 double]
%       StimCode: [nb_epochs x1 double]
%        Trigger: [nb_samples x1 double]
%         Target: [nb_samples x1 double]
%     EpochClass: [nb_epochs x1 double]
%       Channels: [nb_samples x nb_Channels*nb_Users+1 double]
%              h: [1x1 struct]
%           name: 'Session1'

% OUPUT dataPlayers structure format for each Player
%            Fs: scalar
%        StimPos: [nb_epochs*nb_sessions x1 double]
%       StimCode: [nb_epochs*nb_sessions x1 double]
%        Trigger: [nb_samples*nb_sessions x1 double]
%         Target: [nb_samples*nb_sessions x1 double]
%     EpochClass: [nb_epochs*nb_sessions x1 double]
%         Channels: [nb_samples*nb_sessions x nb_Channels]
%              h: [1x nb_sessions struct]
%           name: 'Session1'
%     Conditions: [nb_epochs*nb_sessions x1 double]
%   Stimulations: [nb_samples*nb_sessions x1 double]
%
% *** History: 2015-07-07
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015

if nargin<2
    nb_subjects=2;
end

dataPlayers=struct;


for indUser=1:nb_subjects
    dataPlayers(indUser).Fs=data.session(1).Fs;
    dataPlayers(indUser).Fs=data.session(1).Fs;
    nbchannels=(size(data.session(1).Channels,2)-1)/2;
    userChannels=(nbchannels*(indUser-1)+1):nbchannels*indUser;
    dataPlayers(indUser).Channels=[];
    dataPlayers(indUser).EpochClass=[];
    dataPlayers(indUser).StimPos=[];
    dataPlayers(indUser).Conditions=[];
    dataPlayers(indUser).Stimulations=[];
    dataPlayers(indUser).Trigger=[];
    
    for indSession=1:length(data.session)
        %concatenate all the sessions
        dataPlayers(indUser).h(indSession)=data.session(indSession).h;
        dataPlayers(indUser).Channels=[dataPlayers(indUser).Channels; data.session(indSession).Channels(:,userChannels)];
        dataPlayers(indUser).Trigger=[dataPlayers(indUser).Trigger;data.session(indSession).Trigger];
        [tmp1 tmp2 tmp3 tmp4]=Multisubject_MARTI_triggers(data.session(indSession).Channels(:,65),indUser);
        dataPlayers(indUser).EpochClass=[dataPlayers(indUser).EpochClass ; tmp1];
            dataPlayers(indUser).StimPos=[dataPlayers(indUser).StimPos ;tmp2];
            dataPlayers(indUser).Conditions=[dataPlayers(indUser).Conditions ; tmp3];
            dataPlayers(indUser).Stimulations=[dataPlayers(indUser).Stimulations; tmp4];
    end
    
    
    
end