clear all
Directory= 'D:\data\erwan\'; %change path if needed
%Directory ='D:\Users\Oui-Oui\Documents\Mes Documents GIPSA\MATLAB\';
load( [Directory 'FullErwanData_RAW.mat'])
%%
BASEDIR='D:\data\Data-New-Marco-All\'
BASEDIR2='E:\Data-New-Marco-All4\'

Users=Generate_Users_Numbers(1:24)
f1=1
N=4
decimationN=4;
userNB=4
for f2=[20]
    
    
    for userNB=1:24
        user=Users{userNB};
        USERDIR = fullfile(BASEDIR,user)
        USERDIR2 = fullfile(BASEDIR2,user)
        
        listdir = dir(USERDIR);
        Nsession = length(listdir)-2;
        Type={'online','training'}
        e=[];
        for k=1:Nsession
            for Phase=0:1
                for mod=1:2
                    try
                        tic
                        close all
                        FULLDIR = fullfile(USERDIR,['Session' num2str(k)],['Phase' num2str(Phase)],Type{mod});
                        FULLDIR2 = fullfile(USERDIR2,['Session' num2str(k)],['Phase' num2str(Phase)],Type{mod});
                        
                        %Directory='D:\data\Data-New-Marco-All\01\Session1\Phase0\Online\';
                        %filenameEEG=['\EEG_' num2str(f1), '_', num2str(f2)  'Hz.txt'];
                        %filenameSTIM=['\STIMnew.txt'];
                        
                        EEG=data(userNB).session{k}.phase{(Phase+1)}.(Type{mod});
                                                fs=EEG.Fs;

                        %E=EEG.s;
                        %Flash=EEG.Flash;% useless

                                                
                        %{
                        [FFT freq]=Spectra(E(fs*120:end,:),fs);
                        figure(2);semilogy(freq,mean(FFT,2))
                        
                        %E=[zeros(fs*10,size(E,2));E];
                        
                        figure(1)
                        subplot(411)
                        time=(1:size(E,1))/fs;
                              %                  E=E-convn(E,ones(fs/2,1)/fs*2,'same');
                        Econv=convn(E,ones(fs,1)/fs,'same');
                        plot(time,Econv(1:end,2)-E(1:end,2))
                        plot(time,E(1:end,2));hold all;
                        plot(time,Econv(1:end,2));plot(time,Econv(1:end,2)-E(1:end,2))
                        hold off
                                                %E=E-Econv;
                        
                                                time=(1:size(E,1))/fs;

                        subplot(412)
                        plot(time,E(1:end,2));hold all;
                        drawTag(Flash,fs)
                        
                        hold on;               plot(time,E(1:end,2));hold off
                        axis([0 150 -500 500]);xlabel('time (s)');ylabel('amplitude (µV');title('RAW data')
                        legend('EEG signal','flash')
                        %f1=5
                        %}
                        [EEG.s EEG.Fs EEG.Flash]=preprocessingEEG(EEG.s,EEG.Fs,[f1 f2 N decimationN],EEG.Flash);
                        %{
                        subplot(413)
                        time=(1:size(Ef,1))/fsf;

                        
                        drawTag(flashf,fsf)
                        hold on;               plot(time,Ef(1:end,2)');hold off
                        axis([0 150 -500 500]);xlabel('time (s)');ylabel('amplitude (µV');title('BandPass 1-12Hz')
                        legend('EEG signal','flash')
                           %}
                        
                        %indFlash=find(flashf); % OLD MANUAL (useless)
                        %stim=[indFlash,EEG.Y]; % OLD MANUAL (useless)
                        
                        %[FFT freq]=Spectra(Ef(fsf*60:end,:),fsf);
                        %figure(2);hold all;semilogy(freq,mean(FFT,2));hold off;legend('RAW','Filtered');axis([0 55 0 10]);title('Power Spectrum');xlabel('Freq (Hz)');ylabel('PSD µV ^2/Hz')
                        %figure(1);subplot(313);plot(time,conv(Ef(1:end,2),[1 1 1 1 1 1]/6,'same'))
                        toc
                        %%
                        
                        %%
                        %eegplot(Ef(:)')
                        %{
                        %%%%%% MANUAL OLD VERSION
                        dlmwrite([FULLDIR filenameEEG],Ef(:)','delimiter',' ', 'precision','%.6f')
                        writeStim( [FULLDIR filenameSTIM], stim )
                        dlmwrite([FULLDIR2 filenameEEG],Ef(:)','delimiter',' ', 'precision','%.6f')
                        writeStim( [FULLDIR2 filenameSTIM], stim )
                        %}
                        %%%%%% AUTONOMOUS VERSION
                        EEG_mat2txt([FULLDIR2 '\'], EEG)
                        toc
                    catch error
                        disp(['ERROR with user ' num2str(userNB)])
                        continue
                    end
                end
            end
        end
    end
end
