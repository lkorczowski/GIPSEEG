%this script is made to extract some function from the toolbox and put it
%in a external file.
%it has also the option to automatically add need information
restoredefaultpath
run('D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\installer.m')
%functions_to_copy={'ACSTP','CSTP', 'applyCSTP', 'EnsembleAverage', 'meanOverlap', 'epoch_p300', 'WeightsEstimation',...
%    'CorrectLatency', 'ConvergenceLatencies', 'LatencyEstimation', 'best_Pz', 'CSTPinfo', 'installer','script_ACSTP','preprocessingEEG','removeEpochArtifact',...
%	'printCurrentState','plotEEG','global_field_power','count','applyWeights','script_ACSTP_test','FindScaling','findElectrodes','script_ACSTP_testARNAUD',...
 %   'EEG_cnt2mat','Position2Trigger'}
 
 functions_to_copy={'ACSTP', 'installer','script_ACSTP_test', 'ACSTPshow','pop_acstp','eegplugin_acstp','ACSTPprint','script_test_pop_acstp'}
	
	%functions_to_copy={'MARTI_gdf2mat','Generate_Users_Numbers','load_EEG_data'}
%% MANUAL MODIFICATIONS (it opens all the function and you can replace by ctrl+F)
%{
for indFct=1:length(functions_to_copy)
        edit(functions_to_copy{indFct})
end
%}

%% AUTOMATIC MODIFICATIONS


output_directory='D:/Mes Documents GIPSA/MATLAB/SVN/ACSTP/';
mkdir(output_directory)
mkdir(output_directory,'function')
mkdir(output_directory,'script')

foldersup='';
for indFct=1:length(functions_to_copy)
    clear line_to_modify history_found A line_history line_author line_work relatedwork_found author_found history_found
    loc=which(functions_to_copy{indFct});
    if strcmp(functions_to_copy{indFct},'installer') 
        foldersup='';
    elseif strcmp(functions_to_copy{indFct}(1:4),'scri')
        foldersup='script/';
    else
       foldersup='function/'; 
    end
    %copyfile(loc,'D:/24/');
    
    % Read txt into cell A
    fid = fopen(loc,'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    history_found=[];
    while ischar(tline)
        if ~isempty(tline)
            if length(tline)>10
        history_found(i)=strcmp(tline(1:11),'% *** Histo');
        author_found(i)=strcmp(tline(1:11),'% *** Autho');
        relatedwork_found(i)=strcmp(tline(1:11),'% *** Relat');
            end
        is_comment(i)=strcmp(tline(1),'%');
        end
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    
    end
    fclose(fid);
    line_history=find(history_found,1);
    line_author=find(author_found,1);
    line_work=find(relatedwork_found,1);
    line_to_modify=[line_history,line_author,line_work];
    %%
    
    if length(line_to_modify)<3
        disp([functions_to_copy{indFct} ' has no history']);
        edit(functions_to_copy{indFct})
    else
        Date=date;
    A{line_history}=['% *** History: ' Date];
     A{line_author}=['% *** Author: Louis KORCZOWSKI, GIPSA-Lab, ' Date(end-3:end)];
      A{line_work}=['% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)']; 

    end
    % Change cell A
    %A{69} = sprintf('%d',99);
    % Write cell A into txt
    if ~exist(output_directory)
        mkdir(output_directory);
        disp(['Folder ' output_directory ' created'])
    end
    fid = fopen([output_directory foldersup functions_to_copy{indFct} '.m'], 'w');
    for p = 1:numel(A)
        p;
        if A{p+1} == -1
            fprintf(fid,'%s', A{p});
            break
        else
            fprintf(fid,'%s\n', A{p});
        end
    end
    
end
restoredefaultpath
disp('copy over, path restored')
clear all