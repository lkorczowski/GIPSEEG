function RND2=schuffle_indices_classwise(RND,Y,NbClasses)

Class=unique(Y);
NbC=length(Class);
if NbC~=length(NbClasses)
    warning('error in number of classes in schuffle_indices_classwise.m')
    RND2=[];
    return
end
ChosenTrials=[];

   
for z=1:NbC
    if NbClasses(z)>0
        ChosenTrials=[ChosenTrials;find((Y(RND))==Class(z),NbClasses(z))];
    end
end
    RND2=RND(ChosenTrials);
    
