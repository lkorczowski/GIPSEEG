%% load the data EKAT 
clear all
Directory= 'D:\data\Hyperscanning\EKATE\Groups\'
load([Directory 'ALLgroups.mat'])
%% compute the P300

BImapping=[1,2,5,3,6,12,14,16,21,22,23,24,25,27,28,29];
BImapping=[BImapping [1,2,5,3,6,12,14,16,21,22,23,24,25,27,28,29]+32];

for i=1:length(ALLdata.Xte)
P{i}=Epoch_average(ALLdata.Xte{i}(:,:,:),ALLdata.Yte{i},2,'multiP1');
Pclean{i}=Epoch_average(ALLdata.Xte{i}(BImapping,:,~ALLdata.isBad{i}),ALLdata.Yte{i}(~ALLdata.isBad{i}),2,'multiP1');
end


%% plot and save P300 in Pz and Cz
close all
for group=19
P1(:,:,group)=mean(ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}(:,:,:)==1)),3);
P0(:,:,group)=mean(ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}(:,:,:)==0)),3);
end
P1GA=mean(P1,3);%grand average
P0GA=mean(P0,3);%grand average
    winElec=[14,21,22,23,24,25,27,28,29]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
figure
hold all;
plotEEG(P1GA(winElec,:))
plotEEG(P0GA(winElec,:))
hold off

figure
hold all;
h1=plotEEG(P1GA([14 21],:),0.2,128,{'Cz' 'Pz'},40,[.1 0.5 1])
h0=plotEEG(P0GA([14 21],:),0.2,128,{'Cz' 'Pz'},40,[1 0.5 .1])

hold off
xlabel('Time (s)')
legend([h1(1) h0(1)],{'K^+', 'K^-'},'location','southeastoutside')
legend('boxoff')
          %  set(gcf, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
            set(gcf, 'color', [1 1 1])
saveas(gcf,'D:\Mes Documents GIPSA\Présentations & Communications\images et logo\P300_TA_NT.png')

%% plot and save P300 covariance matrices
close all
for group=1:19
    X_TA{group}=ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==1));
    X_NT{group}=ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==0));
P1(:,:,group)=mean(ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==1)),3);
P0(:,:,group)=mean(ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==0)),3);
end
map=[1,0,0;
    1,.1,.1;
    1,.2,.2;
    1,.3,.3;
    1,.4,.4;
    1,.5,.5;
    1,.6,.6;
    1,.7,.7;
    1,.8,.8;
    1,.9,.9;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1];
P1GA=mean(P1,3);%grand average
P0GA=mean(P0,3);%grand average
    figure
    group=11
    trial=34
    Xk=[P1(1:32,:,group)/norm(P1(1:32,:,group),'fro') ;X_TA{group}(1:32,:,trial)/norm(X_TA{group}(1:32,:,trial),'fro') ;...
        P1(33:64,:,group)/norm(P1(33:64,:,group),'fro') ;X_TA{group}(33:64,:,trial)/norm(X_TA{group}(32:64,:,trial),'fro') ];
    COVk=cov(Xk');
    COVk(1:5,1:5)
    clims = [1e-5 2.e-4];
    imagesc(abs(COVk),clims)
    colormap(flipud(map))
%set(gca,'YDir','inverse')
            set(gcf, 'color', [1 1 1])
        set(gca,'XtickLabel',[],'YtickLabel',[])
hold on;plot([64.5,64.5],[0,128],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([0,128],[64.5,64.5],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([32.5,32.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[32.5,32.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([96.5,96.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[96.5,96.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
saveas(gcf,'D:\Mes Documents GIPSA\Présentations & Communications\images et logo\COV_TA.png')
figure
trial=33
    Xk=[P1(1:32,:,group)/norm(P1(1:32,:,group),'fro') ;X_NT{group}(1:32,:,trial)/norm(X_NT{group}(1:32,:,trial),'fro') ;...
        P1(33:64,:,group)/norm(P1(33:64,:,group),'fro') ;X_NT{group}(33:64,:,trial)/norm(X_NT{group}(32:64,:,trial),'fro') ];
    COVk=cov(Xk');
        COVk(1:5,1:5)

       imagesc(abs(COVk),clims)
    colormap(flipud(map))
%set(gca,'YDir','inverse')
            set(gcf, 'color', [1 1 1])
                    set(gca,'XtickLabel',[],'YtickLabel',[])
hold on;plot([64.5,64.5],[0,128],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([0,128],[64.5,64.5],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([32.5,32.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[32.5,32.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([96.5,96.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[96.5,96.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;

saveas(gcf,'D:\Mes Documents GIPSA\Présentations & Communications\images et logo\COV_NT.png')
%%
%% plot and save P300 MEAN covariance matrices
close all
for group=1:19
    X_TA{group}=ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==1));
    X_NT{group}=ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==0));
P1(:,:,group)=mean(ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==1)),3);
P0(:,:,group)=mean(ALLdata.Xte{group}(:,:,find(ALLdata.Yte{group}==0)),3);
end
map=[1,0,0;
    1,.1,.1;
    1,.2,.2;
    1,.3,.3;
    1,.4,.4;
    1,.5,.5;
    1,.6,.6;
    1,.7,.7;
    1,.8,.8;
    1,.9,.9;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1;
    1,1,1];
P1GA=mean(P1,3);%grand average
P0GA=mean(P0,3);%grand average
    figure
    group=11
    for trial=1:50
    
    Xk(:,:,trial)=[P1(1:32,:,group)/norm(P1(1:32,:,group),'fro') ;X_TA{group}(1:32,:,trial)/norm(X_TA{group}(1:32,:,trial) ,'fro') ;...
        P1(33:64,:,group)/norm(P1(33:64,:,group),'fro') ;X_TA{group}(33:64,:,trial)/norm(X_TA{group}(33:64,:,trial),'fro')];
    end
    Xk_cov=covariances(Xk);
    
    COVk=mean_covariances(Xk_cov,'ld');
    COVk(1:5,1:5)
    clims = [1e-6 1e-4];
    imagesc(abs(COVk),clims)
    colormap(flipud(map))
%set(gca,'YDir','inverse')
            set(gcf, 'color', [1 1 1])
        set(gca,'XtickLabel',[],'YtickLabel',[])
hold on;plot([64.5,64.5],[0,128],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([0,128],[64.5,64.5],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([32.5,32.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[32.5,32.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([96.5,96.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[96.5,96.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
saveas(gcf,'D:\Mes Documents GIPSA\Présentations & Communications\images et logo\COVmean_TA.png')
%
figure
trial=33

     for trial=1:50
    
    Xk(:,:,trial)=[P1(1:32,:,group)/norm(P1(1:32,:,group),'fro') ;X_NT{group}(1:32,:,trial)/norm(X_NT{group}(1:32,:,trial) ,'fro') ;...
        P1(33:64,:,group)/norm(P1(33:64,:,group),'fro') ;X_NT{group}(33:64,:,trial)/norm(X_NT{group}(33:64,:,trial),'fro')];
    end
    Xk_cov=covariances(Xk);
    
    COVk=mean_covariances(Xk_cov,'ld');
    COVk(1:5,1:5)
    clims = [1e-6 1e-4];
    imagesc(abs(COVk),clims)
    colormap(flipud(map))
%set(gca,'YDir','inverse')
            set(gcf, 'color', [1 1 1])
                    set(gca,'XtickLabel',[],'YtickLabel',[])
hold on;plot([64.5,64.5],[0,128],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([0,128],[64.5,64.5],'Color',[0 0 0],'LineWidth',1);hold off;
hold on;plot([32.5,32.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[32.5,32.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([96.5,96.5],[0,128],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;
hold on;plot([0,128],[96.5,96.5],'--','Color',[.8 .8 .8],'LineWidth',1);hold off;

saveas(gcf,'D:\Mes Documents GIPSA\Présentations & Communications\images et logo\COVmean_NT.png')