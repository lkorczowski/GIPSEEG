function AUsers=Generate_Components_Numbers(Numbers)
for i=Numbers
    if i<10
        AUsers{i}=['c0' num2str(i)];
    else
        AUsers{i}=['c' num2str(i)];
    end
end
AUsers=AUsers(find(~cellfun(@isempty,AUsers)));