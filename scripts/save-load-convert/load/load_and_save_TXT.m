BASEDIR='D:\data\Data-New-Marco-All\'

Users=Generate_Users_Numbers(1:24)
 
for userNB=[1:24]
    user=Users{userNB};
USERDIR = fullfile(BASEDIR,user)

listdir = dir(USERDIR);
Nsession = length(listdir)-2;
Type={'Online','Training'}
e=[];
for k=1:Nsession
for phase=0:1
for mod=1:2
        try
            tic
    FULLDIR = fullfile(USERDIR,listdir(k+2).name,['Phase' num2str(phase)],Type{mod});
        File = fullfile(FULLDIR,'EEG.txt')


%Directory='D:\data\Data-New-Marco-All\01\Session1\Phase0\Online\';
fs=128;
formatSpec = '%s';
[Text]=fopen([FULLDIR '\EEG.txt'])
EEG = textscan(Text,formatSpec,...
    'Delimiter', ' ', ...
    'CollectOutput', true);

tmp=cellfun(@str2num,EEG{:},'UniformOutput',false);
Null=cellfun(@isempty,tmp);
tmp{~Null};
E=[tmp{~Null}];
%%
%[FFT,Freq]=Spectra(E,fs);
%figure;subplot(211);plot(Freq,log(FFT))
  wo=50/(fs/2);bw=wo/35;
        [b,a]=iirnotch(wo,bw);
        s = filtfilt(b,a,E);
        %[FFT,Freq]=Spectra(s,fs);
        %xlabel('Frequency (Hz)');title('Before Notch')
        %subplot(212);plot(Freq,log(FFT))
        %xlabel('Frequency (Hz)');title('After Notch')
        
        %%

dlmwrite([FULLDIR '\EEGnotch.txt'],s,'delimiter',' ', 'precision','%.6f')
toc
        catch error
            e(k,phase,mod)=error
        end
end
end
end
end