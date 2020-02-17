% SCRIPT FOR Approximate Joint Diagonalization Comparison
% @FLorent Bouchard, 2018

%% INPUTS AND PARAMETERS
% ################ FOLDERS ##################
clear all
% INPUT/OUTPUT FOLDER (including .\data_original\ and .\unmix\)
MAIN_FOLDER='D:\DATA\ssvep_dat\';
OUT_FOLDER=[MAIN_FOLDER 'figures\'];if ~exist(OUT_FOLDER);mkdir(OUT_FOLDER);end
% ################# PARAMETERS ##################
ssvep_freqs=[13,17,21];
fmax=256; %max frequency (Hz)
fres=0.5; %frequency resolution (Hz)
LB=1;%low cut-off freq (Hz)
HB=43;%hight cut-off freq (Hz)
FILT_param=[LB,HB, 4, NaN]; %[LB, HB, order, decimation(not yet computed)]
event_code=[33024,33025,33027,33026]; %resting state,  13Hz, 17Hz, 21Hz
event_table={'rest','13 Hz','17 Hz','21 Hz'};
L_epoch=3; %length windows in seconds
L_welch=2; %length welch estimator windows in seconds
l_welch=1.75; %length welch estimator step in seconds
freqs=[LB:fres:HB];%freqs for welch estimator
offset=2.5;% WARNING : using 2.5seconds for experiment specific sync for event_pos


%% MAIN LOOP (subject, session, method)
close all;
SUBJECTS=[1:7];SESSIONS=[1:2];
for indsub=SUBJECTS; %subject count
    for indses=SESSIONS; %session count
        
        % ################## loading ########################"
        load([MAIN_FOLDER 'data_original\subject0' ...
            num2str(indsub) '_session0' num2str(indses) '.mat'],...
            'raw','event_pos','event_type', 'chnames','fs');
        %load unmixing matrices B_all
        load([MAIN_FOLDER 'unmix\unmix_subject0' ...
            num2str(indsub) '_session0' num2str(indses) '_all_div.mat'],...
            'B_all');
        
        % ################## conversion ########################"
        N_methods = size(B_all,3);%number of methods
        fs=double(fs);
        event_pos=event_pos(ismember(event_type,event_code));%keep only useful event_pos
        event_type=event_type(ismember(event_type,event_code));%keep only useful event_type
        for i=1:size(chnames,1);ElectrodesName{i}=chnames(i,:);end

        for indmet=1:N_methods;%methods count
            
            OUT_FOLDER=[MAIN_FOLDER 'figures\s' num2str(indsub) '-' num2str(indses) '\m' num2str(indmet) '\'];
            if ~exist(OUT_FOLDER);mkdir(OUT_FOLDER);end
            
            Triggers=zeros(size(raw,2),1);Triggers(event_pos)=1; %build triggers
            FILT_param(4)=round(double(fs)/fmax); %compute decimation factor
            
            % filtering (electrode & source domain)
            [xfilt, xfs, xtrig]=preprocessingEEG(raw', fs,FILT_param,Triggers); %Band-pass filt
            X=xfilt';
            B=B_all(:,:,indmet); %extracting method-specific spatial unmixing matrix
            S=B*X; %compute sources
            
            % epoching (WARNING : OFFSET IS SET TO 2.5 seconds due to experiment-specific triggers' sync)
            Xk=epoch_p300(X,xtrig,L_epoch*xfs,xfs*offset); %electrode domain epoching
            Sk=epoch_p300(S,xtrig,L_epoch*xfs,xfs*2.5); %source domain epoching
            
            % compute power spectral density with welch-estimator
            % NOTE : while Fourrier Transform is linear (i.e. Sfk=B*Xfk), I'm not sure for welch.
            % Therefore I compute welch independently for electrode and source
            % domain (just in case but slower)
            Xfk=zeros(size(Xk,1),length(freqs),size(Xk,3)); % electrodes spectrum
            Sfk=zeros(size(Sk,1),length(freqs),size(Sk,3)); %source spectrum
            for k=1:size(Xk,3);
                Xfk(:,:,k)=pwelch(Xk(:,:,k)',L_welch*xfs,l_welch*xfs,freqs,xfs)';
                Sfk(:,:,k)=pwelch(Sk(:,:,k)',L_welch*xfs,l_welch*xfs,freqs,xfs)';
            end
            
            % Retroprojection of sources in electrodes domain (to avoid scaling issue)
            A=pinv(B); %estimated mixing matrix
            Xtilde=zeros([size(X),size(A,2)]); % backprojected sources [N x Ttot x P] (P: nb sources)
            Xtilde_k=zeros([size(Xk) size(A,2)]); % backprojected epochs [N x T x K x P]
            Xtilde_fk=zeros(size(Xk,1),length(freqs),size(Xk,3)); % backprojected spectrum [N x T x K x P]
            
            for indA=1:size(A,2)
                Atilde=A;
                Atilde(:,1:end ~= indA)=0;%indA-th source only
                Xtilde(:,:,indA)=Atilde*S;
                Xtilde_k(:,:,:,indA)=epoch_p300(Xtilde(:,:,indA),xtrig,L_epoch*xfs,xfs*offset);
                for k=1:size(Xtilde_k,3);
                    Xtilde_fk(:,:,k,indA)=pwelch(Xtilde_k(:,:,k,indA)',L_welch*xfs,l_welch*xfs,freqs,xfs)';
                end
            end
            
            %% ############## PLOT MEAN SPATIAL SPECTRUM ###############
            %electrode domain
            figure;
            plot(freqs,([...mean(mean(Xf(1,:,event_type==33024),3),1)',...
                mean(mean(Xfk(:,:,event_type==33025),3),1)',... %13 Hz
                mean(mean(Xfk(:,:,event_type==33027),3),1)',...%17 Hz
                mean(mean(Xfk(:,:,event_type==33026),3),1)']),...%32 Hz
                'linewidth',2);
            legend(event_table(2:end));
            vline([6.5 26],'b--');vline(13,'b');
            vline([8.5 34],'r--');vline(17,'r');
            vline([10.5 42],'y--');vline(21,'y');
            print(gcf,[OUT_FOLDER 's' num2str(indsub) '-' num2str(indses)],'-dpng','-r300')
            %eegplot(X)
            
            %source domain (nb figures = P)
            for indA=1:size(A,2)
                figure;
                plot(freqs,([...mean(mean(Xf(1,:,event_type==33024),3),1)',...
                    mean(mean(Xtilde_fk(:,:,event_type==33025,indA),3),1)',... %13 Hz
                    mean(mean(Xtilde_fk(:,:,event_type==33027,indA),3),1)',... %17 Hz
                    mean(mean(Xtilde_fk(:,:,event_type==33026,indA),3),1)']),... %21Hz
                    'linewidth',2);
                legend(event_table(2:end));
                vline([6.5 26],'b--');vline(13,'b');
                vline([8.5 34],'r--');vline(17,'r');
                vline([10.5 42],'y--');vline(21,'y');
                print(gcf,[OUT_FOLDER 's' num2str(indsub) '-' num2str(indses) '_source' num2str(indA)],'-dpng','-r300')
                close
            end
            
            %% ############## TOPOPLOT ###############
            % WARNING : eeglab REQUIRED AND findElectrodes from LK_toolbox
            
            figure
            nb_columns=4; % <-------------- HARDCODED FIGURE WIDTH
            nb_rows=2; % <-------------- HARDCODED FIGURE HEIGHT
            FontSize=12;
            TemporalCells=0;%only required when ploting temporal signals, so set to 0
            CompLABELS=Generate_Components_Numbers(1:size(A,2))';
            [~, localization]=findElectrodes(ElectrodesName, ElectrodesName);
            for i=1:size(A,2)
                SpatialCells=(1:nb_columns*nb_rows);
                SpatialCells=SpatialCells(~ismember(SpatialCells,TemporalCells));
                subplot(nb_rows,nb_columns,[SpatialCells(i)]);
                topoplot((A(:,i)),localization,'plotrad',0.7,'headrad',0.5,'colormap',colormap(jet(20)));title(CompLABELS{i},'FontSize',FontSize);
            end
            set(gcf, 'PaperPosition', [0 0 15 7.5]);
            hb=colorbar;
            set(hb, 'position', [0.92 0.3 0.02 0.6],'FontSize',FontSize);
            print(gcf,[OUT_FOLDER 's' num2str(indsub) '-' num2str(indses) '_source_topolots'],'-dpng','-r300')
            
            %% ############## SCORES SAVE ###############
            for i=1:size(Xtilde_fk,4)
                [scoreA(:,i),scoreB(:,:,i)]=frequency_ratio(Xtilde_fk(:,:,:,i));
                scoreB_mean(:,i)=mean(scoreB(:,:,i),2);
            end
            
            for indA=1:size(A,2)
                figure;
                plot(freqs,([...mean(mean(Xf(1,:,event_type==33024),3),1)',...
                    mean(scoreB(:,event_type==33025,indA),2),... %13 Hz
                    mean(scoreB(:,event_type==33027,indA),2),... %17 Hz
                    mean(scoreB(:,event_type==33026,indA),2)]),... %21Hz
                    'linewidth',2);
                legend(event_table(2:end));
                vline([6.5 26],'b--');vline(13,'b');
                vline([8.5 34],'r--');vline(17,'r');
                vline([10.5 42],'y--');vline(21,'y');
                print(gcf,[OUT_FOLDER 's' num2str(indsub) '-' num2str(indses) '_source' num2str(indA) '_f-ratio'],'-dpng','-r300')
                close
            end
            %%
            close all
            figure
            for indA=1:size(A,2)
                ratio13(indA)=mean(scoreB(ismember(freqs, 13),event_type==33025,indA),2);%13 Hz
                ratio17(indA)=mean(scoreB(ismember(freqs, 17),event_type==33027,indA),2); %17 Hz
                ratio21(indA)= mean(scoreB(ismember(freqs, 21),event_type==33026,indA),2); %21 Hz
            end
            subplot(131);plot(sort(ratio13,'descend'),'b*-');title(event_table(2));ylim([0 1])
            subplot(132);plot(sort(ratio17,'descend'),'r*-');title(event_table(3));ylim([0 1])
            subplot(133);plot(sort(ratio21,'descend'),'y*-');title(event_table(4));ylim([0 1])
            
            print(gcf,[OUT_FOLDER 's' num2str(indsub) '-' num2str(indses) '_source_f-ratio'],'-dpng','-r300')
        end
    end
end
