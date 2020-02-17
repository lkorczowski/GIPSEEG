

function [indicesP, namesP, cvec]=struc2res(Pall,Param,cargs)

indicesP={};
nameP=fieldnames(Pall);
nameArgs=fieldnames(Param);
nArgs = length(nameArgs);
nameArgs=fieldnames(Param);
cvec=[];
for indO = 1: nArgs
    inpName = nameArgs{indO};
    largs(indO)=length(Param.(inpName));
    if nargin<3
    cargs{indO}=1:largs(indO);
    end
    if indO==1
        cvec=cargs{indO};
    else
        cvec=combvec(cvec,cargs{indO});
    end
end

for indC=1:size(cvec,2)
    indicesP{indC}=ones(size(Pall));
    for indO = 1: nArgs
        
        inpName = nameArgs{indO};
        if ischar(Param.(inpName){1})
            indicesP{indC}=indicesP{indC} & strcmp(Param.(inpName){cvec(indO,indC)},{Pall.(inpName)});
            namesP(indC).(inpName)=Param.(inpName){cvec(indO,indC)};
        else
            indicesP{indC}=indicesP{indC} & cmpcell({Pall.(inpName)},Param.(inpName){cvec(indO,indC)});
            namesP(indC).(inpName)=Param.(inpName){cvec(indO,indC)};
        end
        
    end
end

