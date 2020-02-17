clear all
Directory='D:\Mes Documents GIPSA\Présentations & Communications\2016-05 BCI Meeting 2016\'
load([Directory 'data\results.mat'])
close all
Paper=[0 0 25 15];
PaperFactor=1;
FontSize=24
XLegend1={'A','B','C','D','E','F'};

[Y errY]=gp_median_quantile({bi2012a.acc bi2013a.acc bi2014a.acc bi2014b.acc bi2015a.acc(:) bi2015b.acc(:,3:4)})

figure
Colors=cellfun(@(x) x/10-.1,{[5 10 10],[10. 3. 3.],[3.33 3.33 1.33]*2},'Uniform',0)

hbar=barwitherr(errY,Y)
for indMethods=1:length(hbar)
    set(hbar(indMethods),'FaceColor',Colors{indMethods});
end
    set(gcf, 'PaperPosition', Paper*PaperFactor,'units','normalized','outerposition',[0.45 0.1 0.5 .9])
                set(gca,'xtick',1:length(XLegend1),'xticklabel',XLegend1,'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');

                print(gcf,[Directory 'data/summary'],'-dtiff','-r450')