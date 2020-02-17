%% LOAD THE DATA FOR TEST
close all
clear all
Directory='D:\data\Hyperscanning\BI-multiplayers\Groups\'
        load([Directory 'ALLgroups.mat'])
        load( [Directory 'results_AUC_STATS_all_VS_intra_vsMAPPING_P0_25.mat'],'BImap')
%% PREPARE EPOCHS & PLOT IT
X=ALLdata.Xte{1}(BImap{1},:,:);
Y=ALLdata.Yte{1};

        P=[10,20]
        [Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(X,Y,P,[],[],1)
        
        
        figure
        subplot(121)
        plotEEG(mean(X(:,:,Y==0),3),5,128)
        subplot(122)
        plotEEG(mean(X(:,:,Y==1),3),5,128)
        
        figure
        subplot(121)
        plotEEG(mean(Xtr(:,:,Ytr==0),3),5,128)
        subplot(122)
        plotEEG(mean(Xtr(:,:,Ytr==1),3),5,128)
       
        
        %% riemannian classification test
        cov_method={};
        P.nusers=2
        P.P300_ref_orientation='multiP1'
        P.GroupsName='01'
        P.Stats='ALL'
        P.method_mean='ld'
        P.method_dist=[]
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,P.nusers,P.P300_ref_orientation,cov_method);
%eegplot([P1(:,:,1);P1(:,:,2)],'srate',128,'winlength',1)
%return
COVte = covariances_p300_hyper(Xtest,P1,P.nusers,P.P300_ref_orientation,cov_method);

disp(['Classification initialized on Group ' P.GroupsName])
%%% Classification par MDM riemannienne
[Yestimated, distances,C,COV,EIG] = mdm_hyper(COVte,COVtr,Ytraining,P.nusers,P.Stats,P.P300_ref_orientation,P.method_mean,P.method_dist,0);
        ConfM=confusionmat(Ytest,Yestimated)
        count(Yestimated==Ytest)/length(Ytest)
        figure
        plotroc(Ytest',-diff(distances'))
        
        %% CSP + LDA 1 from Delorme
        % data epoching, class-grouping and CSP
        %training CSP
        Class=unique(Y);
        nof=3
        wnd=1:128
for k = 1:length(Class)
    EpochsClass=Xtraining(:,:,Ytraining==Class(k));
    %EPO{k} = mean(EpochsClass,3);
    EPO{k} = EpochsClass(:,:);
end

[V,D] = eig(cov(EPO{2}'),cov(EPO{1}')+cov(EPO{2}'));
S = V(:,[1:nof end-nof+1:end]);
 
% log-variance feature extraction and LDA
for k = 1:2
    feat{k} = squeeze(log(var(reshape(EPO{k}'*S, length(wnd),[],2*nof))));
end



w = ((mean(feat{2})-mean(feat{1}))/(cov(feat{1})+cov(feat{2})))';%hyperplan equation
b = (mean(feat{1})+mean(feat{2}))*w/2;%biais

T=eye(size(Xtest,2));
for x=1:size(Xtest,3)
    Features(:,:,x)=Xtest(:,:,x)'*S;
y(x) = log(var(T*Features(:,:,x)))*w - b;
end

figure
subplot(121)
plotEEG(mean(Features(:,:,Ytest==0),3)')
subplot(122)

plotEEG(mean(Features(:,:,Ytest==1),3)')

        Yestimated=double(y>0)'
        ConfM=confusionmat(Ytest,Yestimated)
        count(Yestimated==Ytest)/length(Ytest)
        figure
        plotroc(Ytest',y)
        
        
             %% CSP + LDA from another source
        close all
        clear Features
        Class=unique(Y);
        nof=6
        wnd=1:128
        X1=Xtraining(:,:,Ytraining==1);X0=Xtraining(:,:,Ytraining==0);
[unmixing] = csp(X0(:,:), X1(:,:),nof)
%[unmixing] = csp(mean(Xtraining(:,:,Ytraining==0),3),mean(Xtraining(:,:,Ytraining==1),3),3)
 EPO{1}=X0(:,:); EPO{2}=X1(:,:);
% log-variance feature extraction and LDA
for k = 1:2
    feat{k} = squeeze(log(var(reshape(unmixing*EPO{k}, length(wnd),[],2*nof))));
end

w = ((mean(feat{2})-mean(feat{1}))/(cov(feat{1})+cov(feat{2})))';%hyperplan equation
b = (mean(feat{1})+mean(feat{2}))*w/2;%biais

T=eye(size(Xtest,2));
for x=1:size(Xtest,3)
    Features(:,:,x)=unmixing*Xtest(:,:,x);
y(x) = log(var(Features(:,:,x)*T,[],2))'*w - b;
end

figure
subplot(121)
plotEEG(mean(Features(:,:,Ytest==0),3))
subplot(122)

plotEEG(mean(Features(:,:,Ytest==1),3))

        Yestimated=double(y>0)'
        ConfM=confusionmat(Ytest,Yestimated)
        count(Yestimated==Ytest)/length(Ytest)
                [PerfX,PerfY,~,AUC,OPTROCPT] = perfcurve(Ytest',y',1);

        figure
        
        plotroc(Ytest',y)
        
        %% SWLDA
        clear Xtr_c Xte_c
        for indK=1:size(Xtr,3)
        for indC=1:size(Xtr,1)
        Xtr_c(indC,:,indK) = decimate(Xtr(indC,:,indK),4,1);
        end
        end
         for indK=1:size(Xte,3)
        for indC=1:size(Xte,1)
                Xte_c(indC,:,indK) = decimate(Xte(indC,:,indK),4,1);
        end
         end
             [weights, intercept] = train_swlda(Xtr,Ytr);
              y1 = test_swlda(Xte,weights,intercept);
              count((sign(y1)>0)==Yte)/40;
              auc = area_under_curve(y1',Yte')
              [weights, intercept] = train_swlda(Xtr_c,Ytr);
              y1 = test_swlda(Xte_c,weights,intercept);
              count((sign(y1)>0)==Yte)/40;
              auc = area_under_curve(y1',Yte')
              %% reglda
              y=reglda(Xte,Xtr,Ytr)