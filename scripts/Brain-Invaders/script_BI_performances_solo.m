%% LOAD FILES WITH MAT
clear all;
close all
% select the players in the analysis
nb_players=71;
AUsers=Generate_Users_Numbers(1:nb_players);
%AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
%Users = {'11','15','18','21'}; % selected players
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
FigureDirec=[Directory 'Figure\'];
load([Directory 'ALLsolo.mat'])
%% prepare data
nusers=2;
allcombi=nchoosek(AUsers,nusers);
%
auc=[];
P=0.25;
%% loop generate data
for Trial=1:6 %generate and test data sets many times in order to 
for Col=1:length(allcombi) 
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        Xte=ALLdata.Xte(IND);Yte=ALLdata.Yte(IND);
        [Xout, Yout,RND]=Multisubject_generator(ALLdata.Xte(IND),ALLdata.Yte(IND),RND);
        [X, Y]=Randomize_trials(Xout, Yout);
        [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X,Y,P);
        %% classification
        
        
        %elec=[1,2,6,7,8,9,10,11,12,13,14,15,16];
               elec=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
        % Estimation des matrices de covariance spÃ©ciale P300
        [COVtr, P1] = covariances_p300(Xtraining(elec,:,:),Ytraining);
        COVte = covariances_p300(Xtest(elec,:,:),P1);
        
        %condiWoutRegu=cond(COVte(:,:,1))
        %condiWRegu=cond(COVte(:,:,1)+10e-5*eye(size(COVte,1)))
        
        %%% Classification par MDM riemannienne
        [Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining);
        
        Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
        disp(['Performance classification using closest reference :' num2str(Perfclassif) '%']);
        
        % matrice de confusion
        disp('matrice de confusion');
        disp(confusionmat(Ytest,Yestimated));
        
        % courbe roc visualisée
        %{
        scores = -diff(distances')';
        [TPR,FPR,TH] =roc(Ytest',scores');
        % Area under curve (higer is better, max=1);
        [X,Y,~,AUC,OPTROCPT] = perfcurve(Ytest',scores',1);
        AUC
        figure;plot(X,Y);
        xlabel('False positive rate'); ylabel('True positive rate')
        title('ROC for classification by logistic regression')
        %}
        
        % courbe roc sauvergardée
        scores = -diff(distances')';
        [X,Y,~,AUC(Col),OPTROCPT] = perfcurve(Ytest',scores',1);
        disp(['Area Under Curve : ' num2str(AUC(Col))]);
        figure;f1=plotroc(Ytest',scores');
        %legend('random',cat(2,allcombi{Col,:}))
        legend('random',char(allcombi(Col,:)))
        text(0.5,0.5,['AUC=' num2str(AUC(Col),4)],'FontSize',12)
        saveas(f1,[FigureDirec 'ROC' cat(2,allcombi{Col,:}) 'trial' num2str(Trial+1) '.jpeg'])
        close
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end
%
results.AUC=AUC;
results.allcombi=allcombi;
save([Directory 'ROC' num2str(nb_players) 'trial' num2str(Trial+1) '.mat'],'results')
%}
end
%% 
clear AUC
for Trial=1:6
    load([Directory 'ROC' num2str(nb_players) 'trial' num2str(Trial) '.mat'],'results');
    AUC(Trial,:)=results.AUC;
end
clear results
results.AUC=mean(AUC,1);
results.maxAUC=max(AUC,[],1);
results.StdAUC=std(AUC,2);
results.allcombi=allcombi;

save([Directory 'ROC' num2str(nb_players) '.mat'],'results')

%%
close all
figure
subplot(411)
errorbar(results.AUC,results.StdAUC)
subplot(412)
hist(results.AUC,15)
subplot(413)
hist(results.AUC(results.StdAUC<0.04),15)
subplot(414)
hist(results.maxAUC,15)
figure
hist(results.StdAUC,15)