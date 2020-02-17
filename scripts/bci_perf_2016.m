Res=[]
clear AUC STD
for indPl=1:length(allParam.users)  %players (figure)
    
    
    
    % new figure
    
    for testR=1:length(allParam.TestSetRatio) %training set ratio (factor)
        
        [indicesP, namesP, cvec]=struc2res(Pall,Param,{1:2,indPl,testR,2}) %struc2res extracts specific results from the structure of parameters Pall
        for indI=1:length(indicesP) % for all combinaison
            AUC{indI,testR}=[Pall(indicesP{indI}).AUC];
            STD{indI,testR}=[Pall(indicesP{indI}).STD];
        end
        
    end
    subplot(1,3,indPl)
    errorbar(repmat([allParam.TestSetRatio{:}]',1,size(AUC,1)),cellfun(@mean,AUC'),cellfun(@mean,STD')/3);legend({namesP.method_mean},'location','south')
    xlabel('training set (%)')
    ylabel('AUC')
    title(['users' num2str(Param.users{indPl})])
    axis([0 0.6 0.75 1])
    
    
    
end 