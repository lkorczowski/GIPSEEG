function RND2=schuffle_indices_classwise(RND,Y,NbClasses)

Class=unique(Y);
NbC=length(Class);
if NbC~=length(NbClasses)
    disp('error in number of classes in schuffle_indices_classwise.m')
    return
end
ChosenTrials=[];

   
for z=1:NbC
   ChosenTrials=[ChosenTrials;find((Y(RND))==Class(z),NbClasses(z))];
end
    RND2=RND(ChosenTrials);
    
