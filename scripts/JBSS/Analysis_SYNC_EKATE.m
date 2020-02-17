clear all
Sources={}


% PROPRE
sourcesUtilesPasPropres={[11	7	3	14];
    [2	8	14	10	];
    [2	6	3	]		;
    [7	6	11	2	1	13];
    [4	7	6	1	9	3];
    [9	2	3	1]	;
    [2	6	13	]	;
    [4	8	11	1]	;
    [2	9	5	]	;
    [14	3	2	12]	;
    [5	3	10	2]	;
    [5	3	10	2]	;
    [5	3	10	2]	;
    [7	1	2	3	8]	;
    [7	4	11	5]	;
    [2	5	15	1]	;
    [5	7	1]	;
    [5	6	1]	;
    [1	3	9]		};
BSS={'JBSS','BSS'}
for TYPEBSS=1:2
    for Test_Training_Set_ratio=[0.1 0.2 0.3 0.4 0.5];
        setBalance=1;
        clear COH
        for nbusers=1:19
            disp(['LOAD user' num2str(nbusers)])
            clear results EpochClass
            if TYPEBSS==1
                disp('JBSS')
                Directory='D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\Gp'
                Directory=[Directory num2str(nbusers) '\']
                fullDir=[Directory 'results' num2str(nbusers) '5.mat']
            elseif TYPEBSS==2
                disp('BSS')
                Directory='D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\BSS\Gp'
                Directory=[Directory num2str(nbusers) '\']
                fullDir=[Directory 'resultsBSS' num2str(nbusers) '.mat']
            else
                disp('WARNING TYPEBSS INCORRECT')
            end
            load(fullDir)
            resultssave=results;
            for SelectedSources=[-1:0];
                for indTrial=1:100
                    clc
                    results=resultssave;
                    disp(['Group' num2str(nbusers)])
                    disp(['SelectedSources' num2str(SelectedSources)])
                    disp(['indTrial' num2str(indTrial)])
                    
                    clear Xtr Ytr Xte Yte tmp tmp2
                    EpochClass{nbusers}=results.OJoB.sources(1).trialinfo;
                    RND=[];
                    for indPlayer=1:2
                        disp(['indPlayer' num2str(indPlayer)])
                        tmp=results.OJoB.sources(indPlayer).trial;
                        results.OJoB.sources(indPlayer).trialOLD=results.OJoB.sources(indPlayer).trial;
                        tmp2=[cell2mat3D(tmp)];
                        EpochClass{nbusers}(EpochClass{nbusers}==2)=0;
                        [Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(tmp2,EpochClass{nbusers},Test_Training_Set_ratio,RND,[],setBalance);
                        
                        
                        switch SelectedSources
                            case -1
                                sUt=sourcesUtilesPasPropres{nbusers};
                            case 0
                                sUt=1:size(Xtr,1);
                            otherwise
                                sUt=SelectedSources;
                                
                        end
                        Sources=results.OJoB.sources(indPlayer).trial;
                        SourcesEpochs{nbusers,indPlayer}=cell2mat3D(Sources);
                        %SourcesEpochs_usefull{nbusers,indPlayer}=Xtr(sUt,:,:);
                        
                        results.OJoB.sources(indPlayer).trial= mat3D2cell(Xtr(sUt,:,:));
                        
                        
                    end
                    
                    if SelectedSources==0
                        %compute average sources
                        
                        COH.TA{nbusers,indTrial}(:,:,indPlayer)=mean(SourcesEpochs{nbusers,indPlayer}(:,:,Ytr==1),3);
                        COH.NT{nbusers,indTrial}(:,:,indPlayer)=mean(SourcesEpochs{nbusers,indPlayer}(:,:,Ytr==0),3);
                        COH.EpochClass{indTrial}=Ytr;
                    end
                    
                    results.OJoB.sources(1).trialinfo=Ytr;
                    results.OJoB.sources.trial;
                    %% SALEEEE
                    
                    cfg=[];
                    cfg.method     = 'global';
                    cfg.nfft       = 512;
                    cfg.freqRange  = [0 128/2];
                    cfg.conditions = [0;1];
                    cfg.window=chebwin(0.6*)
                    Cohtmp=coherence(cfg, results.OJoB.sources);
                    COH.Sources{SelectedSources+2}=sUt;
                    COHuser(:,:,:,:,indTrial)=single(Cohtmp.all);
                    [COH.rhoTA{SelectedSources+2}(nbusers,:,indTrial) COH.rho_instTA{SelectedSources+2}(nbusers,:,indTrial) COH.rho_lagTA{SelectedSources+2}(nbusers,:,indTrial) ] = dependence2_tmp(squeeze(Cohtmp.C(:,:,:,1)));
                    [COH.rhoNT{SelectedSources+2}(nbusers,:,indTrial) COH.rho_instNT{SelectedSources+2}(nbusers,:,indTrial) COH.rho_lagNT{SelectedSources+2}(nbusers,:,indTrial) ]= dependence2_tmp(squeeze(Cohtmp.C(:,:,:,2)));
                    %[rho,rho_instant,rho_lag]
                    %{
                    figure(1)
                    
                    plot(COH.rhoTA{SelectedSources+2}(nbusers,:,indTrial))
                    hold all;
                    plot(COH.rhoNT{SelectedSources+2}(nbusers,:,indTrial))
                    hold off
                    %}
                end
                COH.coherence{nbusers,SelectedSources+2}=mean(COHuser,5);
                COHuser=[];
            end
            %{
        figure(2)
        plot(mean(COH.rhoTA{1}(:,:,indTrial),1))
        hold all;
        plot(mean(COH.rhoNT{1}(:,:,indTrial),1))
        hold off
            %}
        end
        %%
        %{
    figure(3)
    subplot(211)
    plot(mean(mean(COH.depTA,3),1))
    hold all;
    plot(mean(mean(COH.depNT,3),1))
    hold off
    subplot(212)
    plot(mean(mean(COH.depTA_allsources,3),1))
    hold all;
    plot(mean(mean(COH.depNT_allsources,3),1))
    hold off
        %}
        save(['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\' BSS{TYPEBSS} '\COH2_' num2str(Test_Training_Set_ratio*10) '.mat'],'COH')
    end
end
%% RESULTS ANALYSIS (load)
clear all
load(['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\COH2_' num2str(1) '.mat'],'COH')

%%
groups=5
Freqs=freqs_val(128,size(COH.rhoTA{1}(groups,:,:),2));
figure(3)
subplot(211)
plot(Freqs,mean(mean((COH.rhoTA{1}(groups,:,:)),3),1))
hold all;
plot(Freqs,mean(mean((COH.rhoNT{1}(groups,:,:)),3),1))
hold off
legend('TA','NT')
subplot(212)
plot(Freqs,mean(mean(COH.rhoTA{2}(groups,:,:),3),1))
hold all;
plot(Freqs,mean(mean(COH.rhoNT{2}(groups,:,:),3),1))
hold off
legend('TA','NT')