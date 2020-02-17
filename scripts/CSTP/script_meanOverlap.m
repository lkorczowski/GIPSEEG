clear all
load EEG_ex

%meanOverlap(EEG.Channels',EEG.Trigger,128)

StimPos=find(EEG.Trigger);
ERP1=epoch_p300(EEG.Channels',EEG.Trigger,128);
tic
for indSection=1:length(StimPos)/6 %for every section (row / column)
Trigger{indSection}=zeros(size(EEG.Trigger));
Trigger{indSection}(StimPos((indSection-1)*6+1:(indSection)*6))=1;%find the position of the section's trials
ERP(:,:,indSection)=meanOverlap(EEG.Channels',Trigger{indSection},128);%average 
end
toc
StimPos=find(EEG.Trigger);
StimPos=StimPos(EEG.EpochClass==1);
Trigger2=zeros(size(EEG.Trigger));
Trigger2(StimPos)=1;%find the position of the section's trials
ERP2(:,:,indSection)=meanOverlap(EEG.Channels',Trigger2,128);%average 


subplot(121);
plotEEG(mean(ERP1,3));
subplot(122)
plotEEG(mean(ERP,3));


subplot(121);
plotEEG(mean(ERP1(:,:,EEG.EpochClass==1),3));
subplot(122);
plotEEG(mean(ERP2,3));
