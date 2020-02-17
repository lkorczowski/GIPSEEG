fig = figure;
Fboot
subplot1(8,2)%,'Gap', [.00 .00], 'XTickL', 'Margin', 'YTickL', 'Margin');
ticklab = {'1', '1-2', '2-5', '5-10', '10-50', '50-100'};
w = {'NT','TA'};
s = {'1','2','3','4','5','6','7','8'};
for indUser=1:8
for indBoot=1:6
BarPlot(:,indBoot,:,indUser)=cat(1,(mean(mean(Fboot(indUser,indBoot).EA_rmse,1),2)),(mean(mean(Fboot(indUser,indBoot).EAcstp_rmse,1),2)));
end
end
for wix=1:2 % class (columns)
    for six=1:8 % users (rows)
%         tmp = [perf.Cond1.OJoB(:,six,wix,6)+perf.Cond2.OJoB(:,six,wix,5)+perf.Cond2.OJoB(:,six,wix,3),perf.Cond2.NOJoB(:,six,wix,6)+perf.Cond2.NOJoB(:,six,wix,5)+perf.Cond2.NOJoB(:,six,wix,3),perf.Cond2.JNoLAD(:,six,wix,6)+perf.Cond2.JNoLAD(:,six,wix,5)+perf.Cond2.JNoLAD(:,six,wix,3)];
        tmp = [BarPlot(:,:,wix,six)]';
        subplot1(wix+(six-1)*2);
        b=bar(tmp);
        b(1).FaceColor = [0.2,0.2,0.2];
        
        b(2).FaceColor = [1,1,1];
        b(3).FaceColor = [0.7,0.7,0.7];
        axis([0.5 6.5 0 100]);
        if six == 7 && wix ==1
            ax = gca;
            ax.XTick = 1:6;
            ax.XTickLabel = ticklab;
            ax.XTickLabelRotation = 90;
            y=ylabel('succès (%)');
            set(y, 'Units', 'Normalized', 'Position', [-0.45, 0.7, 0]);
            xlabel('conditionnements');
        else
            ax = gca;
            ax.XTick = 1:6;
            ax.XTickLabel = {};
            ax.YTickLabel = {};
        end
        if six==1
            ax = gca;
            ax.XAxisLocation = 'top';
            xlabel(w(wix),'Interpreter','latex');
            if wix==3
                xlabel({'\makebox[4in][c]{poids inter-sujets $w_{inter}$}','\makebox[4in][c]{0,7}'},'Interpreter','latex');
            end
        end
        if wix==5
            ax = gca;
            ax.YAxisLocation = 'right';
            ylabel(s(six),'Interpreter','latex');
            if six==4
                ylabel({'2','SNR'},'Interpreter','latex');
            end
        end
    end
end

%# set size of figure's "drawing" area on screen
set(gcf, 'Units','centimeters', 'Position',[0 0 55 60]);
% 
% %# set size on printed paper
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 55 60])
% %# WYSIWYG mode: you need to adjust your screen's DPI (*)
% set(gcf, 'PaperPositionMode','auto');
% 
print(fig,'perfAll','-dpng');