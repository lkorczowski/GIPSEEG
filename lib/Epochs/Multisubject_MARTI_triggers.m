function [Y flashpos Conditions Stimulations]= Multisubject_MARTI_triggers(events, player)
%% Function which outputs .Y vector containing 0 (NT flash) and 1 (TA flash) for different conditions
% Input: vector with event-codes (length: nb samples) (for instance
% s(:,65))
% Input: 1 or 2 (for player)
% Output: .Y vector containing 0 (NT flash) and 1 (TA flash) for different
% occasions (length:nbflash)
% Ouptut: flashpos, determining position of each flash (regardless class(
% NT/TA1/TA2)) within event vecotr (length: nbsamples)


% works good.
ct=1;

ConditionsCodes=[107 108 109 110]; %COOP sync nosync COMP sync nosync
for indCond=1:length(ConditionsCodes)
    try
        ConditionPos(indCond)=find(events==ConditionsCodes(indCond),1);
    catch
        ConditionPos(indCond)=-1;
        disp(['ERROR : condition ' num2str(ConditionsCodes(indCond)) ' is missing'])
    end
end

Stimulations=zeros(length(events),1);
for ni=2:length(events)
    if events(ni)==events(ni-1)
    else
        Stimulations(ni)=events(ni);
    end
end

if player == 2
    for ni=2:length(events)
        if events(ni)==events(ni-1)
            %next...
        elseif events(ni)==0
            %next...
            
            % Store current condition
        elseif events(ni) ==107 % COOP sync
            cond= 107;
        elseif events(ni) ==108 % COOP nosync
            cond= 108;
        elseif events(ni) ==109 % COMP sync
            cond= 109;
        elseif events(ni) ==110 % COMP nosync
            cond= 110;
            
        elseif events(ni)>=160 && events(ni)<=175 && (cond == 108 || cond == 110) % Nosync,
            fp(ct,1)=ni;
            ct=ct+1;
        elseif events(ni)>=180 && events(ni)<=195 && (cond == 108 || cond == 110) % Nosync,
            fp(ct,1)=ni;
            ct=ct+1;
        elseif events(ni)>=60 && events(ni)<=75 && (cond == 107 || cond == 109)   % Sync, so TA1 = TA2
            fp(ct,1)=ni;
            ct=ct+1;
        elseif events(ni)>=80 && events(ni)<=95 && (cond == 107 || cond == 109)   % Sync, so TA1 = TA2
            fp(ct,1)=ni;
            ct=ct+1;
        end
    end
    
elseif player ==1
    for ni=2:length(events)
        if events(ni)==events(ni-1)
        elseif events(ni)>=60 && events(ni)<=75
            fp(ct,1)=ni;
            ct=ct+1;
        elseif events(ni)>=80 && events(ni)<=95
            fp(ct,1)=ni;
            ct=ct+1;
        end
    end
    
end




%%

% Flash
% get vector flashpos comparable to session.Flash  , length of nb flashes, 1 =
% flash at sample, 0 = nonflash, but now only for first sample of flash
flashpos=zeros(size(events,1),1);
for ni=2:length(flashpos)
    if (events(ni)==events(ni-1)) || (events(ni) ==38) % the 38 can be from a error in openvibe
    elseif (events(ni)>=20 & events(ni)<=95)
        flashpos(ni)=1;
    elseif events(ni)>=160 & events(ni)<=195
        flashpos(ni)=1;
    end
end


Ylong= zeros(length(events),1);
Ylong(fp)=1;
Y=Ylong(logical(flashpos)); % Finds where fp and flashpos coincide, so which flashes are targets
indflashpos=find(flashpos);
%% give the condition for each flash

[ConditionPos,IndCond]=sort(ConditionPos);
ConditionsCodes=ConditionsCodes(IndCond);
Conditions=zeros(length(indflashpos),1);
for indCond=1:length(ConditionsCodes)
    cond=ConditionsCodes(indCond);
    Conditions(indflashpos>ConditionPos(indCond))=cond;
end


