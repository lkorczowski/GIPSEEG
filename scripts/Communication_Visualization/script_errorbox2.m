

[Y errY]=gp_median_quantile({bi2012a.acc bi2013a.acc bi2014a.acc bi2014b.acc bi2015a.acc bi2015b.acc})

figure
Colors=cellfun(@(x) x/10-.1,{[5 10 10],[10. 3. 3.],[3.33 3.33 1]/5},'Uniform',0)

hbar=barwitherr(errY,Y)
for indMethods=1:length(hbar)
    set(hbar(indMethods),'FaceColor',Colors{indMethods});
end