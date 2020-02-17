function [AUC, scores,distances, ConfM, Perfclassif,Ytest,RND, P1, SHUF]=mdm_chain_hyper(ALLdata,Parameters,P)
%INPUTS
%struct ALLDATA (Xte (eeg epochs), Yte (epochs labels), isBad (artifacts
%labels)

%in struct P : the parameters Parameters
%Method_mean
%method_dist
%nusers
%BImap
%Stats
%GroupsName
%PourTest
%EIG
%P300_ref_orientation
%RND
%% PREPARATION FOR ALL THE CLASSIFIERS
if P.nusers==2
    BImapping=[P.BImap P.BImap+32];
end
i=find(strcmp(Parameters.GroupsName,P.GroupsName));
Ratio=P.TestSetRatio;
if ~isfield(P,'samples')
    P.samples=1:size(ALLdata.Xte{i(1)},2);
end
X=ALLdata.Xte{i}(BImapping,P.samples,:);Y=ALLdata.Yte{i};

if P.ShuffleCouple==1
    [X SHUF]=Shuffle_Set(X,Y,P.SHUF);
else
    SHUF=P.SHUF;
end
%CHECK IF IT IS WITH ONLY a commonP1
if strcmp(P.Stats,'common')
    P.P300_ref_orientation='commonP1';
end


if isfield(ALLdata,'isBad')
NoArte=~ALLdata.isBad{i};%index of epochs without artifacts
else
    NoArte=1:size(X,3);
end
[Xtraining Ytraining Xtest Ytest RND]=Generate_Training_Test_Set(X(:,:,NoArte),Y(NoArte),Ratio,P.RND,P.NbMax,P.SetBalance);
%      [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X,Y,P);


if strcmp(P.Stats,'SWLDA-hyper')
    [weights, intercept] = train_swlda(Xtraining,Ytraining);
    scores = test_swlda(Xtest,weights,intercept);
    Yestimated=double(scores>0);
    Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
    % matrice de confusion
    ConfM=confusionmat(Ytest,Yestimated);
    % courbe roc
    %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
    distances=[];
    % Area under curve (higer is better, max=1);
    %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
    [PerfX,PerfY,~,AUC,OPTROCPT] = perfcurve(Ytest',scores',1);
    fprintf(['AUC Score for EIG(' num2str(P.EIG) '), Stats(' P.Stats '), P300ref(' P.P300_ref_orientation ') : ' num2str(AUC) '\n\n'])
    P1=[];
elseif strcmp(P.Stats,'SWLDA-multi')
    nb_elec=size(Xtraining,1)/P.nusers;
    for indUser=1:P.nusers
        [weights, intercept] = train_swlda(Xtraining((1:nb_elec)+nb_elec*(indUser-1),:,:),Ytraining);
        distances(:,indUser) = test_swlda(Xtest((1:nb_elec)+nb_elec*(indUser-1),:,:),weights,intercept);
    end
    scores=mean(distances,2);
    Yestimated=double(scores>0);
    Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
    % matrice de confusion
    ConfM=confusionmat(Ytest,Yestimated);
    % courbe roc
    %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
    % Area under curve (higer is better, max=1);
    %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
    [PerfX,PerfY,~,AUC,OPTROCPT] = perfcurve(Ytest',scores',1);
    fprintf(['AUC Score for EIG(' num2str(P.EIG) '), Stats(' P.Stats '), P300ref(' P.P300_ref_orientation ') : ' num2str(AUC) '\n\n'])
    P1=[];
elseif strcmp(P.Stats,'SWLDA-solo')
    nb_elec=size(Xtraining,1)/P.nusers;
    for indUser=1:P.nusers
        [weights, intercept] = train_swlda(Xtraining((1:nb_elec)+nb_elec*(indUser-1),:,:),Ytraining);
        distances(:,indUser) = test_swlda(Xtest((1:nb_elec)+nb_elec*(indUser-1),:,:),weights,intercept);
    
    scores(:,indUser)=distances(:,indUser);
    Yestimated(:,indUser)=double(scores(:,indUser)>0);
    Perfclassif(indUser)=length(find(Yestimated(:,indUser)==Ytest))/length(Ytest);
    % matrice de confusion
    ConfM{indUser}=confusionmat(Ytest,Yestimated(:,indUser));
    % courbe roc
    %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
    % Area under curve (higer is better, max=1);
    %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
    [PerfX,PerfY,~,AUC(indUser),OPTROCPT] = perfcurve(Ytest',scores(:,indUser)',1);
   
        fprintf(['AUC Score with player' num2str(indUser) ' for EIG(' num2str(P.EIG) '), Stats(' P.Stats '), P300ref(' P.P300_ref_orientation ') : ' num2str(AUC(indUser)) '\n'])
    end
    fprintf('\n')
    P1=[];
else
    %%classification
    % Estimation des matrices de covariance spéciale P300
    %cov_method='ledoit-wolf';
    cov_method={};
    [COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,P.nusers,P.P300_ref_orientation,cov_method);
    %eegplot([P1(:,:,1);P1(:,:,2)],'srate',128,'winlength',1)
    %return
    COVte = covariances_p300_hyper(Xtest,P1,P.nusers,P.P300_ref_orientation,cov_method);
    
    disp(['Classification initialized on Group ' P.GroupsName])
    %%% Classification par MDM riemannienne
    [Yestimated, distances,C,COV,EIG] = mdm_hyper(COVte,COVtr,Ytraining,P.nusers,P.Stats,P.P300_ref_orientation,P.method_mean,P.method_dist,P.EIG);
    
    if strcmp(P.Stats,'MDM') %%for the solo, it will classify indepently both players
        for indUser=1:length(Yestimated) %Yestimated is a cell with the 2 scores
            Perfclassif(indUser)=length(find(Yestimated{indUser}==Ytest))/length(Ytest);
            % matrice de confusion for each player
            ConfM{indUser}=confusionmat(Ytest,Yestimated{indUser});
            % courbe roc
            %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
            scores(:,indUser)= -diff(distances{indUser}')';
            %[TPR,FPR,TH] =roc(Ytest',scores{indNB}'); %use perfcurve
            
            % Area under curve (higer is better, max=1);
            %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
            [PerfX,PerfY,~,AUC(indUser),OPTROCPT] = perfcurve(Ytest',scores(:,indUser)',1);
            fprintf(['AUC Score with player' num2str(indUser) ' for EIG(' num2str(P.EIG) '), Stats(' P.Stats '), P300ref(' P.P300_ref_orientation ') : ' num2str(AUC(indUser)) '\n'])
        end
        fprintf('\n')
        
    else
        Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
        % matrice de confusion
        ConfM=confusionmat(Ytest,Yestimated);
        % courbe roc
        %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
        scores= -diff(distances')';
        [TPR,FPR,TH] =roc(Ytest',scores');
        
        % Area under curve (higer is better, max=1);
        %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
        [PerfX,PerfY,~,AUC,OPTROCPT] = perfcurve(Ytest',scores',1);
        fprintf(['AUC Score for EIG(' num2str(P.EIG) '), Stats(' P.Stats '), P300ref(' P.P300_ref_orientation ') : ' num2str(AUC) '\n\n'])
    end
end