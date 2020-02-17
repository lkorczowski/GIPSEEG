clear all
Sources={}


% PROPRE
sourcesUtilesPasPropres={[11	7	3	14];%1
    [2	8	14	10	];%2
    [2	6	3	]		;%3
    [7	6	11	2	1	13];%4
    [4	7];%5
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
BSS={'JBSS'}
for TYPEBSS=1
    for Test_Training_Set_ratio=[0.1];
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
            %%
            for SelectedSources=[-1:2];
                %%
                for indTrial=1
                    clc
                    results=resultssave;
                    disp(['Group' num2str(nbusers)])
                    disp(['SelectedSources' num2str(SelectedSources)])
                    disp(['indTrial' num2str(indTrial)])
                    
                    clear Xtr Ytr Xte Yte tmp tmp2 Sources
                    EpochClass{nbusers}=results.OJoB.sources(1).trialinfo;
                    RND=[];
                    for indPlayer=1:2
                        disp(['indPlayer' num2str(indPlayer)])
                        tmp=results.OJoB.sources(indPlayer).trial;
                        results.OJoB.sources(indPlayer).trialOLD=results.OJoB.sources(indPlayer).trial;
                        tmp2=[cell2mat3D(tmp)];
                        EpochClass{nbusers}(EpochClass{nbusers}==2)=0;
                        RND
                        [Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(tmp2,EpochClass{nbusers},Test_Training_Set_ratio,RND,[],setBalance);
                        
                        
                        switch SelectedSources
                            case -1
                                sUt=sourcesUtilesPasPropres{nbusers};
                            case 0
                                sUt=1:size(Xtr,1);
                            otherwise
                                sUt=sourcesUtilesPasPropres{nbusers}(SelectedSources);
                                
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
                    nbSamples=size(results.OJoB.sources(1).trial{1,1},2);
                    cfg=[];
                    cfg.method     = 'global';
                    cfg.nfft       = 512;
                    cfg.freqRange  = [0 128/2];
                    cfg.conditions = [0;1];
                    cfg.window     = ones(nbSamples,1)%;hamming(floor(min(nbSamples)));
                    [Cohtmp COH.Xfft{SelectedSources+2,nbusers,indTrial}]=coherence(cfg, results.OJoB.sources);
                    COH.Sources{SelectedSources+2}=sUt;
                    COHuser(:,:,:,:,indTrial)=single(Cohtmp.all);
                    [COH.rhoTA{SelectedSources+2}(nbusers,:,indTrial) COH.rho_instTA{SelectedSources+2}(nbusers,:,indTrial) COH.rho_lagTA{SelectedSources+2}(nbusers,:,indTrial) ] = dependence2_tmp(squeeze(Cohtmp.C(:,:,:,cfg.conditions==1)));
                    [COH.rhoNT{SelectedSources+2}(nbusers,:,indTrial) COH.rho_instNT{SelectedSources+2}(nbusers,:,indTrial) COH.rho_lagNT{SelectedSources+2}(nbusers,:,indTrial) ]= dependence2_tmp(squeeze(Cohtmp.C(:,:,:,cfg.conditions==0)));
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
        COH.SourcesUtiles=sourcesUtilesPasPropres;
        save(['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\' BSS{TYPEBSS} '\COH5_' num2str(Test_Training_Set_ratio*10) '.mat'],'COH')
    end
end
%% RESULTS ANALYSIS (load)
clear all
load(['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\COH5_' num2str(1) '.mat'],'COH')

%%
groups=5
Freqs=freqs_val(128,size(COH.rhoTA{1}(groups,:,:),2));
figure(3)

SelSources=3
subplot(411)
plot(Freqs,mean(mean((COH.rhoTA{SelSources}(groups,:,:)),3),1))
hold all;
plot(Freqs,mean(mean((COH.rhoNT{SelSources}(groups,:,:)),3),1))
hold off
axis([0 30 0 1])
legend('TA','NT')
if length(groups)==1
    if SelSources==1
title(['Source ' num2str(COH.SourcesUtiles{groups})])
    elseif SelSources==2
       title(['Source all'])  
    else
        title(['Source ' num2str(COH.SourcesUtiles{groups}(SelSources-2))])

    end
end

SelSources=4
subplot(412)
plot(Freqs,mean(mean(COH.rhoTA{SelSources}(groups,:,:),3),1))
hold all;
plot(Freqs,mean(mean(COH.rhoNT{SelSources}(groups,:,:),3),1))
hold off
axis([0 30 0 1])

legend('TA','NT')
if length(groups)==1
    if SelSources==1
title(['Source ' num2str(COH.SourcesUtiles{groups})])
    elseif SelSources==2
       title(['Source all'])  
    else
        title(['Source ' num2str(COH.SourcesUtiles{groups}(SelSources-2))])

    end
end

SelSources=1
subplot(413)
plot(Freqs,mean(mean(COH.rhoTA{SelSources}(groups,:,:),3),1))

hold all;
plot(Freqs,mean(mean(COH.rhoNT{SelSources}(groups,:,:),3),1))
hold off
axis([0 30 0 1])

legend('TA','NT')
if length(groups)==1
    if SelSources==1
title(['Source ' num2str(COH.SourcesUtiles{groups})])
    elseif SelSources==2
       title(['Source all'])  
    else
        title(['Source ' num2str(COH.SourcesUtiles{groups}(SelSources-2))])

    end
end

SelSources=2
subplot(414)
plot(Freqs,mean(mean(COH.rhoTA{SelSources}(groups,:,:),3),1))
hold all;
plot(Freqs,mean(mean(COH.rhoNT{SelSources}(groups,:,:),3),1))
hold off
axis([0 30 0 1])

legend('TA','NT')
if length(groups)==1
    if SelSources==1
title(['Source ' num2str(COH.SourcesUtiles{groups})])
    elseif SelSources==2
       title(['Source all'])  
    else
        title(['Source ' num2str(COH.SourcesUtiles{groups}(SelSources-2))])

    end
end
saveas(gcf,['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\COH4_' num2str(1) '.tiff'])

figure(1)
SelSources=2
IndFreq=7
size(COH.coherence{groups,SelSources})
spTA=COH.coherence{groups,SelSources}(:,:,IndFreq,2);
spNT=COH.coherence{groups,SelSources}(:,:,IndFreq,1);
subplot(121)
imagesc(abs(spTA-eye(size(spTA))))
title('TA')
colorbar('SouthOutside')

subplot(122)
imagesc(abs(spNT-eye(size(spNT))))
title('NT')
colorbar('SouthOutside')
saveas(gcf,['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\COH4_' num2str(1) '_cosp.tiff'])


%%
SelSources=3
    allFFT=[COH.Xfft{SelSources,groups,:}];
fftTA=cell2mat3D(allFFT{2}); %TA
fftNT=cell2mat3D(allFFT{1}); %NT
fftTAmean=abs(mean((fftTA),3));
fftNTmean=abs(mean((fftNT),3));

fftTAreal=mean(real(fftTA),3);
fftNTreal=mean(real(fftNT),3);
fftTAimag=mean(imag(fftTA),3);
fftNTimag=mean(imag(fftNT),3);

figure(10)
subplot(311)
plot(Freqs,fftTAmean)
axis([0 30 0 1])

ylabel('$$\sum_{k}|X_k|$$','Interpreter','latex')
if length(groups)==1
    if SelSources==1
title(['TA Source ' num2str(COH.SourcesUtiles{groups})],'FontSize',24)
    elseif SelSources==2
       title(['TA Source all'])  
    else
        title(['TA Source ' num2str(COH.SourcesUtiles{groups}(SelSources-2))],'FontSize',24)

    end
end
subplot(312)
plot(Freqs,abs(fftTAreal))
ylabel('$$|(\sum_{k}re(X_k)|$$','Interpreter','latex')
axis([0 30 0 1])

subplot(313)
plot(Freqs,abs(fftTAimag))
axis([0 30 0 1])

ylabel('$$|(\sum_{k}im(X_k)|$$','Interpreter','latex')

saveas(gcf,['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\COH4_' num2str(1) '_cosp_source_' num2str(COH.SourcesUtiles{groups}(SelSources-2)) '.tiff'])

figure(11)
subplot(311)
plot(Freqs,fftNTmean)
axis([0 30 0 1])

ylabel('$$\sum_{k}|X_k|$$','Interpreter','latex')
if length(groups)==1
    if SelSources==1
title(['NT Source ' num2str(COH.SourcesUtiles{groups})],'FontSize',24)
    elseif SelSources==2
       title(['NT Source all'])  
    else
        title(['NT Source ' num2str(COH.SourcesUtiles{groups}(SelSources-2))],'FontSize',24)

    end
end
subplot(312)
plot(Freqs,abs(fftNTreal))
ylabel('$$|(\sum_{k}re(X_k)|$$','Interpreter','latex')
axis([0 30 0 1])

subplot(313)
plot(Freqs,abs(fftNTimag))
axis([0 30 0 1])

ylabel('$$|(\sum_{k}im(X_k)|$$','Interpreter','latex')

saveas(gcf,['D:\data\Hyperscanning\EKATE\Groups\Results\FLORENT\JBSS\COH4_' num2str(1) '_cosp_source_' num2str(COH.SourcesUtiles{groups}(SelSources-2)) 'NT.tiff'])


