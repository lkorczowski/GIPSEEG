% working examples for Non-Orthogonal Bilinear Approximate Joint
% Diagonalization by Gauss Planar Transformation (to be submitted)
% *** Authors : L. Korczowski, F. Bouchard, C. Jutten, M. Congedo

% Comparison CAJD versus AJD

clear all;
close all;

%% define global parameters
% <<<<<<< HEAD
nTest = 250;
N     = 8;
T     = 128
K     = 100;
F     = 20;
condA = [2 10];
allCondE = {[2 10]};
allSNR   = {0,10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10^0};
outputDIR = 'D:\Mes Documents GIPSA\MATLAB\figures\NoBAD_eusipco\'
outputDIR='D:\GoogleDrive\MATLAB\figures\thesis\chap6\'
outputFILE=[outputDIR 'CAJD_sim']
NEWSIM=1;%boolean, if 1 new sim is done (erasing the old one), if 0 only plot (faster)
%% perform Test on CAJD model
if NEWSIM
for tix=1:nTest
    disp(tix);
    for eix=1:length(allCondE)
        for six=1:length(allSNR)
            %% define current test parameters
            condE  = allCondE{eix};
            SNR    = allSNR(six);
            
            %% simulate data
            opSim = struct('N',N,'T',T,'K',K,'F',F,'condA',condA,'condE',condE,'SNR',SNR);
            [X,Cf,A,E] = simBAJD_dat(opSim);
            cond(E);
            % get all covariance matrices for AJD
            [C,Ct] = convertDat_BAJD(X,Cf);
%             [U,S,V]=svd(mean(X,3)',0);
%             V=V(:,1:N)
            %% define BSS parameters
            B0      = eye(N);
            D0      = orth(randn(T,N));
            epsilon = 1e-18;
            itMax   = 100;
%             B0      = inv(V)';%eye(N);
%             D0      = pinv(U)';%orth(randn(T,N));
%             epsilon = 1e-18;
%             itMax   = 50;
            
            %% perform CAJD
            opCAJD = struct('B0',B0,'D0',D0,'eps',epsilon,'itMax',itMax,'A',A,'E',E);
            [B_CAJD,D_CAJD,S_CAJD,C_CAJD,info_CAJD] = gp_CAJD_GPT(X,Cf,opCAJD);
            
            %% perform AJD
            opAJD = struct('B0',B0,'eps',epsilon,'itMax',itMax,'A',A);
            [B_AJD,C_AJD,info_AJD] = gp_AJD_GPT(C,opAJD);
            %% perform BAJD
            opBAJD = struct('B0',B0,'D0',D0,'eps',epsilon,'itMax',itMax,'A',A,'E',E);
            [B_BAJD,D_BAJD,S_BAJD,info_BAJD]= gp_BAJD_GPT(X,opBAJD);
            
            %% save results
            % results(tix,eix,six) = struct('CAJD',[],'AJD',[]);
            % results(tix,eix,six).dat = struct('X',X,'Cf',Cf,'A',A,'E',E);
            % results(tix,eix,six).CAJD = struct('conv',[info_CAJD.conv],'B',B_CAJD,'D',D_CAJD,'critOFFC',[info_CAJD.critOFFC],'critOFFS',[info_CAJD.critOFFS],'critA', [info_CAJD.critA],'critE', [info_CAJD.critE]);
            % results(tix,eix,six).AJD = struct('conv',[info_AJD.conv],'B',B_AJD,'critOFF',[info_AJD.critOFF],'critA',[info_AJD.critA]);
            
            results.critA_CAJD{tix,eix,six} = [info_CAJD.critA];
            results.critE_CAJD{tix,eix,six} = [info_CAJD.critE];
            results.critA_BAJD{tix,eix,six} = [info_BAJD.critA];
            results.critE_BAJD{tix,eix,six} = [info_BAJD.critE];
            results.critA_AJD{tix,eix,six} =  [info_AJD.critA];
        end
    end
end
save(outputFILE,'results');
end
%% LOAD DATA 
load(outputFILE)

%% CONVERGENCE plot results for the Amari Criterion on B 
close all

critA_CAJD = results.critA_CAJD;
critE_CAJD = results.critE_CAJD;
critA_BAJD = results.critA_BAJD;
critE_BAJD = results.critE_BAJD;
critA_AJD = results.critA_AJD;

limx = [1 39.5];
limy1 = [ -16 1];
limy2 = [ -16 1];
FontSize=12;
Paper=[0 0 15 10];
LineWidth=1.5;
MarkersSize=6;
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
Legends={'AJD_{inf}','CAJD_{inf}','BAJD_{inf}','AJD_{100}','CAJD_{100}','BAJD_{100}'}
% figure;
% for tix=1:nTest
%     plot(log10(critA_AJD{tix,1,1}),'b');
%     hold on;
%     plot(log10(critA_CAJD{tix,1,1}),'r');
% end
% xlim(limx)
% ylim(limy)
%
% xlabel('number of sweeps');
% ylabel('$I_{M-A}(B)$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
% legend('AJD_{orth}','CAJD_{orth}');
figure
subplot1(1,2,'YTickL','All')
for tix=1:nTest
    subplot1(1)
    plot(log10(critA_AJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{1});
    hold on;
    plot(log10(critA_CAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critA_BAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    
    ylabel('$I_{{M-A}}~~(\mathbf{B})$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
    legend(Legends{1:3},'location','southwest');
    xlim(limx)
    ylim(limy1)
    xlabel('number of sweeps');
    title('(a)')
    
    subplot1(3)
    plot(log10(critA_AJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{1});
    hold on;
    plot(log10(critA_CAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critA_BAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    legend(Legends{4:6},'location','southwest');
    set(gca,'yaxislocation','right');
    %     set(gca,'yticklabel')
    ylim(limy2)
    xlabel('number of sweeps');
    title('(b)')
    
end
xlim(limx)
set(gcf,'paperposition',Paper)

% legend(Legends,'location','southwest');
set(gcf,'color',[1 1 1])
print(gcf, [outputDIR '\fig\_Simu\MoAm'],'-dpng','-r450')


%% CONVERGENCE plot results for the Amari Criterion on D
close all
critA_CAJD = results.critA_CAJD;
critE_CAJD = results.critE_CAJD;
critA_BAJD = results.critA_BAJD;
critE_BAJD = results.critE_BAJD;
critA_AJD = results.critA_AJD;

limx = [1 39.5];
limy1 = [ -16 1];
limy2 = [ -10 1];
FontSize=12;
Paper=[0 0 15 10];
LineWidth=1.5;
MarkersSize=6;
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
Legends={'AJD_{inf}','CAJD_{inf}','BAJD_{inf}','AJD_{100}','CAJD_{100}','BAJD_{100}'};
% figure;
% for tix=1:nTest
%     plot(log10(critA_AJD{tix,1,1}),'b');
%     hold on;
%     plot(log10(critA_CAJD{tix,1,1}),'r');
% end
% xlim(limx)
% ylim(limy)
%
% xlabel('number of sweeps');
% ylabel('$I_{M-A}(B)$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
% legend('AJD_{orth}','CAJD_{orth}');
figure
subplot1(1,2,'YTickL','All')
for tix=1:nTest
    subplot1(1)
    %     plot(log10(critA_AJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{1});
    hold on;
    plot(log10(critE_CAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critE_BAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    
    ylabel('$I_{{M-A}}~~(\mathbf{D})$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
    legend(Legends{2:3},'location','southwest');
    xlim(limx)
    ylim(limy1)
    xlabel('number of sweeps');
    title('(a)')
    
    subplot1(3)
    hold on;
    plot(log10(critE_CAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critE_BAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    legend(Legends{5:6},'location','southwest');
    set(gca,'yaxislocation','right');
    ylim(limy2)
    xlabel('number of sweeps');
    title('(b)')
    
end
xlim(limx)
set(gcf,'paperposition',Paper)

% legend(Legends,'location','southwest');
set(gcf,'color',[1 1 1])
print(gcf, [outputDIR '\fig\_Simu\MoAmE'],'-dpng','-r450')
%% FINAL CONVERGENCE Amari Criterion on B 
close all
clear A_CAJD A_BAJD A_AJD E_CAJD E_BAJD

% plot parameters
limx = [1 39.5];
limy1 = [ -16 1];
limy2 = [ -10 1];
FontSize=12;
Paper=[0 0 15 10];
LineWidth=1.5;
MarkersSize=6;
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
Legends={'CAJD','BAJD','AJD'};

critA_CAJD = results.critA_CAJD;
critE_CAJD = results.critE_CAJD;
critA_BAJD = results.critA_BAJD;
critE_BAJD = results.critE_BAJD;
critA_AJD = results.critA_AJD;

% catch the final value
for tix=1:nTest;for tiy=1;for tiz=1:length(allSNR);
    A_CAJD(tix,tiy,tiz)=critA_CAJD{tix,tiy,tiz}(end);
    E_CAJD(tix,tiy,tiz)=critE_CAJD{tix,tiy,tiz}(end);
    A_BAJD(tix,tiy,tiz)=critA_BAJD{tix,tiy,tiz}(end);
    E_BAJD(tix,tiy,tiz)=critE_BAJD{tix,tiy,tiz}(end);
    A_AJD(tix,tiy,tiz)=critA_AJD{tix,tiy,tiz}(end);
end;end;end;

A_median=10*log10([squeeze(median(A_CAJD)),squeeze(median(A_BAJD)),squeeze(median(A_AJD))]);
E_median=10*log10([squeeze(median(E_CAJD)),squeeze(median(E_BAJD))]);
A_p10=10*log10([squeeze(prctile(A_CAJD,10)),squeeze(prctile(A_BAJD,10)),squeeze(prctile(A_AJD,10))]);
A_p90=10*log10([squeeze(prctile(A_CAJD,90)),squeeze(prctile(A_BAJD,90)),squeeze(prctile(A_AJD,90))]);
E_p10=10*log10([squeeze(prctile(E_CAJD,10)),squeeze(prctile(E_BAJD,10))]);
E_p90=10*log10([squeeze(prctile(E_CAJD,90)),squeeze(prctile(E_BAJD,90))]);

figure
semilogx([allSNR{:}],A_median,'^-','Markers',MarkersSize,'Linewidth',LineWidth);legend(Legends);hold on; ax = gca;
ax.ColorOrderIndex = 1;semilogx([allSNR{:}],A_p10,'--');ax = gca;
ax.ColorOrderIndex = 1;semilogx([allSNR{:}],A_p90,'--');hold off
    ylabel('$I_{{M-M}}~~(\mathbf{B})$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
xlabel('Noise variance $\sigma$', 'interpreter','latex','fontsize',FontSize);

print(gcf, [outputDIR '\fig\_Simu\A_MoMaFinal'],'-dpng','-r450')

figure
semilogx([allSNR{:}],E_median,'^-','Markers',MarkersSize,'Linewidth',LineWidth);legend(Legends);hold on; ax = gca;
ax.ColorOrderIndex = 1;semilogx([allSNR{:}],E_p10,'--');ax = gca;
ax.ColorOrderIndex = 1;semilogx([allSNR{:}],E_p90,'--');hold off
    ylabel('$I_{{M-M}}~~(\mathbf{D})$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
    xlabel('Noise variance $\sigma$', 'interpreter','latex','fontsize',FontSize);

    print(gcf, [outputDIR '\fig\_Simu\E_MoMaFinal'],'-dpng','-r450')
