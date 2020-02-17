function [AUC, scores,distances, ConfM, Perfclassif,Ytest]=mdm_chain_hyper_chrono(ALLdata,Parameters,P)
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
%%
       if P.nusers==2
        BImapping=[P.BImap P.BImap+32];
       end
            i=find(strcmp(Parameters.GroupsName,P.GroupsName));
            Ratio=P.TestSetRatio;
            X=ALLdata.Xte{i}(BImapping,:,:);Y=ALLdata.Yte{i};
            
            NoArte=~ALLdata.isBad{i};%index of epochs without artifacts
            [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X(:,:,NoArte),Y(NoArte),Ratio);
            %      [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X,Y,P);

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
                
                if strcmp(P.Stats,'solo')
                    for indNB=1:length(Yestimated)
                         Perfclassif{indNB}=length(find(Yestimated{indNB}==Ytest))/length(Ytest);
                % matrice de confusion
                ConfM=confusionmat(Ytest,Yestimated{indNB});
                % courbe roc
                                %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
                scores{indNB}= -diff(distances{indNB}')';
                [TPR,FPR,TH] =roc(Ytest',scores{indNB}');
                
                % Area under curve (higer is better, max=1);
                                %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
                [PerfX,PerfY,~,AUC{indNB},OPTROCPT] = perfcurve(Ytest',scores{indNB}',1);
                fprintf(['AUC Score with player' num2str(indNB) ' for EIG(' num2str(P.EIG) '), Stats(' P.Stats '), P300ref(' P.P300_ref_orientation ') : ' num2str(AUC{indNB}) '\n'])
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