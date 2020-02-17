function StimCodeC=triggers_BI2015b_extractONLINE(StimCode)

%% Time relative performance
Conditions=[107 108 109 110]; %COSY CONS CMSY CMNS
Accurary=[105,106,111,112]; %P1 P2 BOTH NONE
StimCodeC=struct;
StimCodeC.StimCode=StimCode;
StimCodeC.COND_stim=ismember(StimCode,Conditions);
StimCodeC.COND_pos=find(StimCodeC.COND_stim);

StimCodeC.StimCode_COND=zeros(length(StimCodeC.StimCode),1);
for indC=1:length(StimCodeC.COND_pos)
    StimCodeC.StimCode_COND(StimCodeC.COND_pos(indC):end)= StimCodeC.StimCode(StimCodeC.COND_pos(indC));
end
StimCodeC.ACC_stim=ismember(StimCode,Accurary);
StimCodeC.ACC_pos=find(StimCodeC.ACC_stim);
StimCodeC.ACC_COND=StimCodeC.StimCode_COND(StimCodeC.ACC_pos);
StimCodeC.nP1=nnz(StimCode==105);
StimCodeC.nP2=nnz(StimCode==106);
StimCodeC.nP12=nnz(StimCode==111);
StimCodeC.nP0=nnz(StimCode==112);
StimCodeC.ACC_online=StimCode(StimCodeC.ACC_stim);
StimCodeC.ACC_table=zeros(length(StimCodeC.ACC_online),3);
for ind=1:length(StimCodeC.ACC_online)
    if StimCodeC.ACC_online(ind)==105 | StimCodeC.ACC_online(ind)==111
        StimCodeC.ACC_table(ind,1)=1;
    end
    if StimCodeC.ACC_online(ind)==106 | StimCodeC.ACC_online(ind)==111
        StimCodeC.ACC_table(ind,2)=1;
    end
    if StimCodeC.ACC_online(ind)==112
        StimCodeC.ACC_table(ind,3)=1;
    end
end
StimCodeC.ACC_table(:,4)=StimCodeC.ACC_COND;

