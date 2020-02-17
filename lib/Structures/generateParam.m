function P=generateParam(Parameters,Index)
% Parameters [struct] with N entry
% Index with 2 columns :
% 1st column : entry name
% 2nd column : entry index
% if entry exists in Parameters but not explitely set in Index,
% the first parameter of the entry is set (default)
%
% EXEMPLE :
% >>P.method={'algebric','geometric'}
% >>P.number=[1,5,6,7];
%
% >>generateParam(P,{'method',2})
% returns :
% ans.method: 'geometric'
% ans.number: 1
%
% >>generateParam(P,{'method',1;'number',4})
% returns :
% ans.method: 'algebric'
% ans.number: 7
for n = fieldnames(Parameters)'
    
        name = n{1};
        tmp=Parameters.(name);
        if isempty(tmp)
                    P.(name) = tmp;
        elseif iscell(tmp)
            I=find(strcmp(Index(:,1),name));
            if isempty(I)
                                P.(name) = tmp{1};

            else
                P.(name) = tmp{Index{I,2}};
            end
        else
            
            I=find(strcmp(Index(:,1),name));
                if isempty(I)
                                    P.(name) = tmp(1);

                else
                P.(name) = tmp(Index{I,2});
                end
                
        end
end


end