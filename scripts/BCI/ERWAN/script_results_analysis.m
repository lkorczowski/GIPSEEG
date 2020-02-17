%% load results
AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
allcombi=results.allcombi;
AUC1=results.AUC1;
AUC=results.AUC;
%% compute analysize
for Col=1:size(allcombi,1)
Users=allcombi(Col,:);
IND=find(ismember(AUsers,allcombi(Col,:)));
AN(Col,1)=mean(AUC1(IND)); %mean of each individu
AN(Col,2)=var(AUC1(IND)); %variance of the individual performance
AN(Col,3)=AUC(Col); %group performance

end
%% analysis ROC
close all
ANc=AN(any(AN(:,3),2),:);
figure
subplot(221)
plot((ANc(:,1)),(ANc(:,3)),'.')
xlabel('Mean of individual AUC ROC')
ylabel('Group AUC ROC')


subplot(222)

plot(sqrt(ANc(:,2)),ANc(:,3)-ANc(:,1),'.')
xlabel('Deviation of individual AUC ROC')
ylabel('Increase in performance (group - indiv) AUC')


subplot(223)

[values, AX] = hist3([(ANc(:,1)),ANc(:,3)],[60 60]);
imagesc(AX{1},AX{2},values');
colorbar
xlabel('Mean of individual AUC ROC')
ylabel('Group AUC ROC')
%axis equal
axis xy
% set(gca,'XScale','log');


subplot(224)

%[values, AX]=hist3([sqrt(ANc(:,2)),log((ANc(:,3)-ANc(:,1)).^2+1)],[50 50]);
[values, AX]=hist3([sqrt(ANc(:,2)),ANc(:,3)-ANc(:,1)],[50 50]);
%-ANc(:,1)
imagesc(AX{1},AX{2},values');
colorbar;
xlabel('Deviation of individual AUC ROC')
ylabel('Increase in performance (group - indiv) AUC')
%set(gca,'YDir','normal')
axis xy

%
%% LOAD FILES WITH MAT
clear all;
close all
AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
Users = {'01','02','03','04'}; % selected players

load('D:\data\ALLdata.mat')
    Xtr=ALLdata.Xtr;
    Xte=ALLdata.Xte;
    Ytr=ALLdata.Ytr;
    Yte=ALLdata.Yte;
 nusers=length(Users);
 allcombi=nchoosek(AUsers,24);
auc=[];
se=1; %selection session for users data 
ph=2; %selection phase for users data

 allcombi=nchoosek(AUsers,nusers);
%%generate multiplayers set
IND=find(ismember(AUsers,Users))

%%generate multiplayers set

[Xout, Yout]=Multisubject_generator(Xtr(IND),Ytr(IND));
[Xtraining, Ytraining]=Randomize_trials(Xout, Yout);

[Xout, Yout]=Multisubject_generator(Xte(IND),Yte(IND));
[Xtest, Ytest]=Randomize_trials(Xout, Yout);

%% pca and save V
close all
clear all
load('D:\data\summary.mat')

P1all=mean(P1(:,:,:),3);
P1all=P1all(randperm(size(P1all,1)),:);
%{
plot(P1all')
COV=cov(P1all');
[V, D]=eig(COV);
[Y I]=sort(diag(D),'descend')
plot(cumsum(Y)./sum(Y))
IND=find((cumsum(Y)./sum(Y))>0.998,1)

Vred=V(:,I(1:IND));

plot((Vred'*P1all)')
%}
V=mypca(P1all,4,1)
plot((V*P1all)')
save('PCA4t.mat','V')
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

