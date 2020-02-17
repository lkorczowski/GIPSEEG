function AUsers=Generate_Users_Numbers(Numbers)
for i=Numbers
    if i<10
        AUsers{i}=['0' num2str(i)];
    else
        AUsers{i}=['' num2str(i)];
    end
end
AUsers=AUsers(find(~cellfun(@isempty,AUsers)));