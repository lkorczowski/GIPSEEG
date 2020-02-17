%script_results_analysis_ROCvsP1dist

%% compute Riemmanien distance between P1
clear all
close all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
load([Directory 'ROC71.mat'])
AUsers=Generate_Users_Numbers(1:71);
nusers=2;
allcombi=nchoosek(AUsers,nusers);
[P1 r m]=zscore(P1,2);
for Col=1:length(allcombi) 
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        P1_riem_dist(Col)=distance(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'),'riemman');
       	Cref=blkdiag(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'));
        C=cov([P1(:,:,IND(1));P1(:,:,IND(2))]');
        P1toP2_riem_dist(Col)=distance(Cref,C,'riemman'); % SYNC_riem
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end
%save([Directory 'P1_riem_dist.mat'],'P1_riem_dist','allcombi','P1toP2_riem_dist')
%% compute analysize
close all
for i=1:length(AUsers)
[IND,INDc]=find(ismember(allcombi,AUsers{i}));
AN(i,1)=mean(P1toP2_riem_dist(IND)); %mean of each individu
AN(i,2)=var(P1toP2_riem_dist(IND)); %variance of the individual performance
%AN(Col,3)=AUC(Col); %group performance

end
figure
subplot(211);plot(AN(:,1)');subplot(212);plot(AN(:,2)');
figure
plot(P1toP2_riem_dist,'*')

%%
close all
ROC=results.AUC;
maxROC=results.maxAUC;
stdROC=results.StdAUC;

allcombiROC=results.allcombi;

figure;plot(AN(:,1),ROC,'*')
xlabel('Mean Geo Dist Synchro')
ylabel('Area under Curve ROC')
hold on
clear Table
for K = 1 : length(allcombiROC) 
text(AN(K,1)*1.001, ROC(K), allcombiROC{K}); 
Table(K,:)=[str2num(allcombiROC{K}) ROC(K) AN(K,1) AN(K,2) maxROC(K) stdROC(K)]; 
end
hold off


%% save results
csvwrite([Directory 'TableROC_SYNC.csv'],Table)


%%

%%
clear all
%close all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
load([Directory 'ROC71.mat'])
%AUsers={'01','02','03','07','08','09','10','11','12','13','14','16','18','19',...
 %   '20','21','22','23','24','25','28','30','35','36','37','38','39',...
  %  '43','44','45','46','47','49','51','52','53','55','56','58','59','61',...
   % '63','64','66'};
   AUsers=Generate_Users_Numbers(1:71)
   subjects={'09','17','26','27','34','58',}
nusers=2;
allcombi=nchoosek(subjects,nusers);
[P1 r m]=zscore(P1,2);
for Col=1:size(allcombi,1) 
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        P1_riem_dist(Col)=distance(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'),'riemann');
       	Cref=blkdiag(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'));
        C=cov([P1(:,:,IND(1));P1(:,:,IND(2))]');
        P1toP2_riem_dist(Col,1)=distance(Cref,C,'fro');%fro norm
        P1toP2_riem_dist2(Col,1)=SYNC_Riem(supressMean(P1(:,:,IND(1))),supressMean(P1(:,:,IND(2))));
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end

%figure;plot(P1toP2_riem_dist,P1toP2_riem_dist2,'*')
%%
clear List listIND
[sortSYNC INDs]=sort(P1toP2_riem_dist,'descend');
combiSYNC=allcombi(INDs,:);
medianSYNC=median(P1toP2_riem_dist)
iter=0;
Alt=1;
listIND=1;
while ~isempty(combiSYNC) && iter<1e4
    if Alt
        tmp=combiSYNC(1,:);
        List(listIND,:)=[tmp , sortSYNC(1), 'high']; 
        listIND=listIND+1;
        [IND trash]=find(ismember(combiSYNC,tmp));
        newList=setdiff((1:size(combiSYNC,1)),unique(IND));
        combiSYNC=combiSYNC(newList,:);
        sortSYNC=sortSYNC(newList);
    else
        tmp=combiSYNC(end,:);
        List(listIND,:)=[tmp , sortSYNC(end), 'low']; 
        listIND=listIND+1;
        [IND trash]=find(ismember(combiSYNC,tmp));
        newList=setdiff((1:size(combiSYNC,1)),unique(IND));
        combiSYNC=combiSYNC(newList,:);
        sortSYNC=sortSYNC(newList);
    end
    
    iter=iter+1;
    Alt=~Alt;
end
%%
for i=1:length(AUsers)
[IND,INDc]=find(ismember(allcombi,AUsers{i}));
AN(i,1)=mean(P1toP2_riem_dist(IND)); %mean of each individu
AN(i,2)=var(P1toP2_riem_dist(IND)); %variance of the individual performance
%AN(Col,3)=AUC(Col); %group performance

end

%%
clear all
%close all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
load([Directory 'ROC71.mat'])
%AUsers={'01','02','03','07','08','09','10','11','12','13','14','16','18','19',...
 %   '20','21','22','23','24','25','28','30','35','36','37','38','39',...
  %  '43','44','45','46','47','49','51','52','53','55','56','58','59','61',...
   % '63','64','66'};
   AUsers=Generate_Users_Numbers(1:71)
   subjects={'09','17','26','27','34','58',}
nusers=2;
allcombi=nchoosek(subjects,nusers);
[P1 r m]=zscore(P1,2);
allcombi={'43' '55';...
'28' '45';...
'12' '49';...
'13' '56';...
'21' '39';...
'18' '63';...
'10' '30';...
'24' '66';...
'01' '19';...
'08' '16';...
'46' '64';...
'23' '51';...
'02' '47';...
'20' '38';...
'25' '37';...
'03' '07';...
'36' '61';...
'11' '59';...
'27' '34'}
for Col=1:size(allcombi,1) 
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        P1_riem_dist(Col)=distance(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'),'riemann');
       	Cref=blkdiag(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'));
        C=cov([P1(:,:,IND(1));P1(:,:,IND(2))]');
        P1toP2_riem_dist(Col,1)=distance(Cref,C,'fro');%fro norm
        P1toP2_riem_dist2(Col,1)=SYNC_Riem(supressMean(P1(:,:,IND(1))),supressMean(P1(:,:,IND(2))));
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end

for i=1:length(P1toP2_riem_dist)
  tmp=allcombi(i,:);
  List(i,:)=[tmp , P1toP2_riem_dist(i)]; 
        
end
