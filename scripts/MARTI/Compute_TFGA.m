%% Compute Grand Average for different measures
% Use output from TF_bitmaps3.m TF.TA and TF.NT , 32x10x256 
cd D:\Stage\Matlab\
addpath(genpath('D:\Stage\Matlab\eeglab13_4_4b' ))  
addpath(genpath('D:\Stage\Matlab\Pilot_multi'),genpath('D:\Stage\Matlab\eeglab13_4_4b'),genpath('D:\Stage\Matlab\biosig\doc'),genpath('D:\Stage\Matlab\biosig\t200_FileAccess') )
%%

Directory='D:\Stage\Matlab\DATA\';
outputfolder='D:\Stage\Matlab\DATA_4AR\2015 BRAIN INVADERS\MARTI_GroupsName.mat';
load(outputfolder) %load GroupsName

%pre-allocate fields with zeros
tCon=zeros(32,10, 256);
tMAmp=zeros(32,10, 256);
tMdir=zeros(32,10, 256);
ntCon=zeros(32,10, 256);
ntMAmp=zeros(32,10, 256);
ntMdir=zeros(32,10, 256);
TA.Con=zeros(4,32,10, 256);
TA.MAmp=zeros(4,32,10, 256);
TA.Mdir=zeros(4,32,10, 256);
NT.Con=zeros(4,32,10, 256);
NT.MAmp=zeros(4,32,10, 256);
NT.Mdir=zeros(4,32,10, 256);
nsubjects=0;

 G2ex=[ 9 10 11 15];

for condition=1:4        
    for indUser= 1:21;
        if sum(indUser==G2ex)>0
            continue
        else
            
        for player= 1:2

      load([Directory 'TFOutput4\G' num2str(indUser) '_p' num2str(player) 'c' num2str(condition) '.mat'] )
        
      %Add new data to TA struct
      tCon=tCon+TF.TA;
      tMAmp=tMAmp+ abs(TF.TA);
      tMdir=tMdir+ angle(TF.TA);
      
      ntCon=ntCon+TF.NT;
      ntMAmp=ntMAmp+ abs(TF.NT);
      ntMdir=ntMdir+ angle(TF.NT);
      
          
                nsubjects=nsubjects+1;

        end

        end

    end
     % Write to TA-struct
      TA.Con(condition,:,:,:)=tCon;
      TA.MAmp(condition,:,:,:)=tMAmp;
      TA.Mdir(condition,:,:,:)=tMdir;
         % Write to NT-struct
      NT.Con(condition,:,:,:)=ntCon;
      NT.MAmp(condition,:,:,:)=ntMAmp;
      NT.Mdir(condition,:,:,:)=ntMdir;
      
      % Make empty again
      tCon=zeros(32,10, 256);
tMAmp=zeros(32,10, 256);
tMdir=zeros(32,10, 256);
ntCon=zeros(32,10, 256);
ntMAmp=zeros(32,10, 256);
ntMdir=zeros(32,10, 256);
end
nsubjects=nsubjects/4;
        
 %%       
 nsubjects=40;
 
    Con(1,:,:,:,:)=abs(TA.Con/nsubjects);
    Con(2,:,:,:,:)=abs(NT.Con/nsubjects);
 Mdir(1,:,:,:,:)=abs(TA.Mdir/nsubjects);
   Mdir(2,:,:,:,:)=abs(NT.Mdir/nsubjects);
 MAmp(1,:,:,:,:)=abs(TA.MAmp/nsubjects);
   MAmp(2,:,:,:,:)=abs(NT.MAmp/nsubjects);
   
    MDir2(1,:,:,:,:)=angle(TA.Con/nsubjects);
    MDir2(2,:,:,:,:)=angle(NT.Con/nsubjects);

%% And plot

 %% Plot phase concentration, emphasis on electrodes differences
 load ([Directory 'Mapping.mat'], 'Mapping' )

     % wavelet parameters
    min_freq = 2;
    max_freq = 20;
    num_frex = 10;
    clim =[0  50 ]
     %chans=[1 1 3 2 2 4 5 3  6 7 13 9 14 10 15 17 22 23 24 20 26 27 28 29 30] ;
     chans=[5 6 9 10  13 15 14 23 18 19 22 24] ;

        lc=length(chans);
for condi=1:4
    for tar=1%:2        
        if tar==1
            tit=[num2str(condi) 'TA'];
        else
            tit=[num2str(condi) 'NT'];
        end
            hfig=figure('Name', tit)
            set(hfig, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.7 0 0.3 1])

            for chani=1:lc
                            subplot(6,2,chani)

        % pz cz centered around cz
        % for each freq normalize for max to be one and min to be 0.
        sig2plot=squeeze(Con(tar,condi,chans(chani),:,:));
    %     norms2p=zeros(num_frex,size(sig2plot,2));
    %         norms2ptemp=zeros(num_frex,size(sig2plot,2));
    % 
    %     for freqi=1:10
    %       mini=min(sig2plot(freqi,:));
    %       
    %       norms2ptemp(freqi,:)=(sig2plot(freqi,:)-mini);
    %       maxi=max(norms2ptemp(freqi,:));
    %        norms2p(freqi,:)=norms2ptemp(freqi,:)/maxi;
    % 
    % 
    %       
    %     end

        % other wavelet parameters
        frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
        time = -1:1/128:1; 
        half_of_wavelet_size = (length(time)-1)/2;
        ytickskip = 2:2:num_frex; % indices into frequencies for which freqs to plot

            imagesc(time,frequencies,sig2plot, clim) %, clim
            set(gca,'ytick',frequencies(ytickskip),'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-1 1])
            xlabel('Time (ms)'), ylabel('Frequency (Hz)')
             title([ Mapping(chans(chani))])

            end
    end
end

 %% Plot Phase concentration, emphasis on condition differences
 load ([Directory 'Mapping.mat'], 'Mapping' )

     % wavelet parameters
    min_freq = 2;
    max_freq = 20;
    num_frex = 10;
    clim =[0  140 ]
     %chans=[1 1 3 2 2 4 5 3  6 7 13 9 14 10 15 17 22 23 24 20 26 27 28 29 30] ;
     %chans=[5 6 9 10  13 15 14 23 18 19 22 24] ;
 % chans=[ 9 10  13 14 15  18 19 ] ;
chans=28;
        lc=length(chans);
        
        for chani=1:lc 
       hfig=figure('Name', [Mapping{chans(chani)}])           
       set(hfig, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.7 0.2 0.3 0.8])

         for tar=1:2        
        

        
for condi=1:4
       subplot(4,2,(condi-1)*2+tar)

        % pz cz centered around cz
        % for each freq normalize for max to be one and min to be 0.
        sig2plot=squeeze(Con(tar,condi,chans(chani),:,:));
    %     norms2p=zeros(num_frex,size(sig2plot,2));
    %         norms2ptemp=zeros(num_frex,size(sig2plot,2));
    % 
    %     for freqi=1:10
    %       mini=min(sig2plot(freqi,:));
    %       
    %       norms2ptemp(freqi,:)=(sig2plot(freqi,:)-mini);
    %       maxi=max(norms2ptemp(freqi,:));
    %        norms2p(freqi,:)=norms2ptemp(freqi,:)/maxi;
    % 
    % 
    %       
    %     end

        % other wavelet parameters
        frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
        time = -1:1/128:1; 
        half_of_wavelet_size = (length(time)-1)/2;
        ytickskip = 2:2:num_frex; % indices into frequencies for which freqs to plot

            imagesc(time,frequencies,sig2plot, clim) %, clim
            set(gca,'ytick',frequencies(ytickskip),'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-1 1])
            xlabel('Time (ms)'), ylabel('Frequency (Hz)')
             title([num2str(condi)])
colormap(gray)
            end
    end
        end

        
         %% Plot mean direction, emphasis on electrodes differences
 load ([Directory 'Mapping.mat'], 'Mapping' )

     % wavelet parameters
    min_freq = 2;
    max_freq = 20;
    num_frex = 10;
    clim =[0  5 ]
     %chans=[1 1 3 2 2 4 5 3  6 7 13 9 14 10 15 17 22 23 24 20 26 27 28 29 30] ;
     chans=[5 6 9 10  13 15 14 23 18 19 22 24] ;

        lc=length(chans);
for condi=1:4
    for tar=1%:2        
        if tar==1
            tit=[num2str(condi) 'TA'];
        else
            tit=[num2str(condi) 'NT'];
        end
            hfig=figure('Name', tit)
            set(hfig, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.7 0 0.3 1])

            for chani=1:lc
                            subplot(6,2,chani)

        % pz cz centered around cz
        % for each freq normalize for max to be one and min to be 0.
        sig2plot=squeeze(Mdir(tar,condi,chans(chani),:,:));
    %     norms2p=zeros(num_frex,size(sig2plot,2));
    %         norms2ptemp=zeros(num_frex,size(sig2plot,2));
    % 
    %     for freqi=1:10
    %       mini=min(sig2plot(freqi,:));
    %       
    %       norms2ptemp(freqi,:)=(sig2plot(freqi,:)-mini);
    %       maxi=max(norms2ptemp(freqi,:));
    %        norms2p(freqi,:)=norms2ptemp(freqi,:)/maxi;
    % 
    % 
    %       
    %     end

        % other wavelet parameters
        frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
        time = -1:1/128:1; 
        half_of_wavelet_size = (length(time)-1)/2;
        ytickskip = 2:2:num_frex; % indices into frequencies for which freqs to plot

            imagesc(time,frequencies,sig2plot,clim) %, clim
            set(gca,'ytick',frequencies(ytickskip),'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-1 1])
            xlabel('Time (ms)'), ylabel('Frequency (Hz)')
             title([ Mapping(chans(chani))])

            end
    end
end

%%         %% Plot mean direction, emphasis on electrodes differences
 load ([Directory 'Mapping.mat'], 'Mapping' )

     % wavelet parameters
    min_freq = 2;
    max_freq = 20;
    num_frex = 10;
    %clim =[0  140 ]
     %chans=[1 1 3 2 2 4 5 3  6 7 13 9 14 10 15 17 22 23 24 20 26 27 28 29 30] ;
     chans=[5 6 9 10  13 15 14 23 18 19 22 24] ;

        lc=length(chans);
for condi=1:4
    for tar=1%:2        
        if tar==1
            tit=[num2str(condi) 'TA'];
        else
            tit=[num2str(condi) 'NT'];
        end
            hfig=figure('Name', tit)
            set(hfig, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.7 0 0.3 1])

            for chani=1:lc
                            subplot(6,2,chani)

        % pz cz centered around cz
        % for each freq normalize for max to be one and min to be 0.
        sig2plot=squeeze(MAmp(tar,condi,chans(chani),:,:));
    %     norms2p=zeros(num_frex,size(sig2plot,2));
    %         norms2ptemp=zeros(num_frex,size(sig2plot,2));
    % 
    %     for freqi=1:10
    %       mini=min(sig2plot(freqi,:));
    %       
    %       norms2ptemp(freqi,:)=(sig2plot(freqi,:)-mini);
    %       maxi=max(norms2ptemp(freqi,:));
    %        norms2p(freqi,:)=norms2ptemp(freqi,:)/maxi;
    % 
    % 
    %       
    %     end

        % other wavelet parameters
        frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
        time = -1:1/128:1; 
        half_of_wavelet_size = (length(time)-1)/2;
        ytickskip = 2:2:num_frex; % indices into frequencies for which freqs to plot

            imagesc(time,frequencies,sig2plot) %, clim
            set(gca,'ytick',frequencies(ytickskip),'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-1 1])
            xlabel('Time (ms)'), ylabel('Frequency (Hz)')
             title([ Mapping(chans(chani))])

            end
    end
end