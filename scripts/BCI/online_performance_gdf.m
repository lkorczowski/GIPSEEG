
clear all
load('MARTI_GroupsName.mat');users=GroupsName;
% load('EKATE_GroupsName.mat');users=ALLgroups;
% users=Generate_Users_Numbers([1:71]);
Directory='D:\data\Hyperscanning\MARTI\Groups'
nsession=4
% dbstop in load_EEG_data at 30 if indU==4

for indU=1:length(users)
    try
        %bi2015a and earlier
        data = load_EEG_data(Directory,['' users{indU}]);
        for indS=1:nsession
            try
                StimCode=data.session(indS).h.EVENT.TYP;
                nrep2(indU,indS)=count(StimCode==33040 | StimCode==1409);
                nrep(indU,indS)=count(StimCode==33279);
                nlvl(indU,indS)=count(StimCode==786);
                nlvl2(indU,indS)=count(StimCode==781);
                nNT(indU,indS)=count(StimCode==33286);
                nTA(indU,indS)=count(StimCode==33285);
                
                nstart(indU,indS)=count(StimCode==100);
                nstop(indU,indS)=count(StimCode==101);
                StimCodeC=struct;
                % ONLY bi2015b BEGIN (65th channel is stimulation code)
                if size(data.session(indS).Channels,2)==65,
                    StimCode2=findpeaks(data.session(indS).Channels(:,65));
                   % StimCodeC(indS)=triggers_BI2015b_extractONLINE(StimCode);
                    %COSY : Cooperation Sync
                    %CONS : Cooperation NonSync
                    %CMSY : Competition Sync
                    %CMNS : Competition NonSync
                    Conditions=[107 108 109 110]; %COSY CONS CMSY CMNS
                    for indCond=1:length(Conditions)
                        [Edges tmp]=find(StimCode2==Conditions(indCond));
                        if length(Edges)==2
                            nP1c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==105);
                            nP2c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==106);
                            nP0c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==112);
                            nP12c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==111);
                        elseif length(Edges)==3
                            nP1c(indU,indS,indCond)=count(StimCode2(Edges(2):Edges(3))==105);
                            nP2c(indU,indS,indCond)=count(StimCode2(Edges(2):Edges(3))==106);
                            nP0c(indU,indS,indCond)=count(StimCode2(Edges(2):Edges(3))==112);
                            nP12c(indU,indS,indCond)=count(StimCode2(Edges(2):Edges(3))==111);
                            disp(['error user' num2str(indU) 'session' num2str(indS) ' condition'  num2str(Conditions(indCond))  ' wrong number of stim'])
                            
                        elseif length(Edges)==4
                            nP1c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==105)+count(StimCode2(Edges(3):Edges(4))==105);
                            nP2c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==106)+count(StimCode2(Edges(3):Edges(4))==106);
                            nP0c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==112)+count(StimCode2(Edges(3):Edges(4))==112);
                            nP12c(indU,indS,indCond)=count(StimCode2(Edges(1):Edges(2))==111)+count(StimCode2(Edges(3):Edges(4))==111);
                            disp(['error user' num2str(indU) 'session' num2str(indS) ' condition'  num2str(Conditions(indCond))  ' twice'])
                            
                        else
                            nP1c(indU,indS,indCond)=NaN;
                            nP2c(indU,indS,indCond)=NaN;
                            nP0c(indU,indS,indCond)=NaN;
                            nP12c(indU,indS,indCond)=NaN;
                            disp(['error user' num2str(indU) 'session' num2str(indS) ' condition' num2str(Conditions(indCond)) ' missing'])
                        end
                    end
                    nP1(indU,indS)=count(StimCode2==105); %tot should be 20
                    nP2(indU,indS)=count(StimCode2==106); %tot should be 20
                    nP0(indU,indS)=count(StimCode2==112); %tot should be 20
                    nP12(indU,indS)=count(StimCode2==111); %tot should be 20
                    % 105 : ONLY Player1 succeeded
                    % 106 : ONLY Player2 succeeded
                    % 112 : NONE succeeded
                    % 111 : BOTH succeeded
                    COSY(indU,indS)=count(StimCode2==107);
                    CONS(indU,indS)=count(StimCode2==108);
                    CMSY(indU,indS)=count(StimCode2==109);
                    CMNS(indU,indS)=count(StimCode2==110);
                    
                end
                % ONLY bi2015b END
                
            catch
                disp(['error user' num2str(indU) 'session' num2str(indS)])
                
            end
        end
    catch
        disp(['error user' num2str(indU)])
    end
end

%% multi martine
tot=nP1+nP2+nP0+nP12 %should be nlvl
totc=nP1c+nP2c+nP0c+nP12c
tot2=sum(totc,3)

perfP1=(nP1+nP12)./tot
perfP2=(nP2+nP12)./tot


perfP1c=nP1c+nP12c
perfP2c=nP2c+nP12c

% perf per session and per condition
perfP1cCO=(perfP1c(:,:,1)./totc(:,:,1)+perfP1c(:,:,2)./totc(:,:,2))/2
perfP2cCO=(perfP2c(:,:,1)./totc(:,:,1)+perfP2c(:,:,2)./totc(:,:,2))/2
perfP1cCM=(perfP1c(:,:,3)./totc(:,:,3)+perfP1c(:,:,4)./totc(:,:,4))/2
perfP2cCM=(perfP2c(:,:,3)./totc(:,:,3)+perfP2c(:,:,4)./totc(:,:,4))/2

%save results
if 1
    clear bi2015b
    bi2015b.nblevels=tot;
    bi2015b.targetPlayers(:,:,1)=(nP1+nP12);
    bi2015b.targetPlayers(:,:,2)=(nP2+nP12);
    bi2015b.accPlayers(:,:,1)=perfP1;bi2015b.accPlayers(:,:,2)=perfP2;
    for indS=1:4
        bi2015b.targetP1= perfP1c;
        bi2015b.targetP2= perfP2c;
        bi2015b.nblevels=totc;
        bi2015b.cond=[107 108 109 110];
        bi2015b.condName={'CoopSync' 'CoopNonSync' 'CompSync' 'CompNonSync'};
        bi2015b.groups=GroupsName;
        bi2015b.session=[1:4];
    end
    
end
%% other (EKATE)
nbrep=nlvl(nlvl(:,1)>8,1)
mean(1./(nbrep(2:end)/9))
std(1./(nbrep(2:end)/9))

diff(find(StimCode==786))
indrow= (find(StimCode==33040))

StimCode(indrow(1):indrow(2))
nrep(nlvl==27)./nlvl(nlvl==27)
perf=(nlvl(nlvl==27)./nrep(nlvl==27))

mean(perf)
std(perf)
[m ind]=min(nrep./nlvl)
% unique(StimCode)
nrep(16)
nrep-nrep2

nNT./nlvl
nTA./nlvl

