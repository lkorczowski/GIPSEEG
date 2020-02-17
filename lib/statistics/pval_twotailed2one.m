function pc=pval_two2one(r,p)

for i=1:length(r)
    if r(i)<0
        pc(i)=p(i)/2;
    else
        pc(i)=1-(1-p(i))/2;
    end
end