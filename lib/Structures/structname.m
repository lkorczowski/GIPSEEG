function name=structname(X)
Ind=1;
for n = fieldnames(X)'
        name{Ind} = n{1};

        Ind=Ind+1;
end
end