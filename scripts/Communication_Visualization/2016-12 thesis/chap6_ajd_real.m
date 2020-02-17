% (0) working examples for Mining the Bilinear Structure of Data by
% Approximate Joint Diagonalization (to be submitted)
% *** Authors : L. Korczowski, F. Bouchard, C. Jutten, M. Congedo

% Comparison AJD, BAJD and CAJD with Gauss Planar Transformation (GPT)
% algorithm

% TUTO
% run LK_toolbox installer
% run gipseep toolbox installer
% modify data path and output
if 1
clear all; close all
outputDIR='D:\GoogleDrive\MATLAB\figures\thesis\chap6\'
load('F:\data\Hyperscanning\EKATE\Groups\ALLgroups.mat')
load('BImap.mat');load('BImapping.mat')
subjects=[1:10]
winElec=[7 9 10 11 12 13 14 15 16];
winTime=[floor((0.05)*128):ceil((0.550)*128)]; %(14)* time window (in sample) used for latency calculation and Pz selection
InitMethods=2; %1: B:eye/D:randn 2: acstp
Xmethods=3; %1: {Xk}, 2: {Xta, Xnt}, 3: {Xbar}bootstrap
%% (1) Compute results for all subjects and SAVE
%
tic
    options.freqs=[1 20]
    options.fs=128
for tix=subjects
    close all
    X=ALLdata.Xte{tix}(BImap{1},:,:);
    Y=ALLdata.Yte{tix};
    N     = size(X,1);
    T     = size(X,2);
    K     = size(X,3);
    F     = 100;

    
    %%preparation and compute ACSTP
    EEG=struct;ACSTPoptions=struct;
    EEG.Channels=transpose(X(:,:))
    tmp=zeros(size(EEG.Channels,1),1);tmp(1:size(X,2):end)=1;
    tmpTA=zeros(size(EEG.Channels,1),1);tmpTA(1:size(X,2):end)=Y;
    EEG.Trigger=tmp;
    EEG.TriggerTA=tmpTA;
    EEG.EpochClass=Y;
    EEG.Fs=128;
    EEG.ElectrodesName=LABELS;
    ACSTPoptions=struct;
    ACSTPoptions.Epoch_size=128;ACSTPoptions.LatencyCorr_max=0;ACSTPoptions.Mask_Electrodes=winElec;
    ACSTPoptions.Mask_Time=winTime;%ACSTPoptions.SubspaceDim=N-2;
    [Xhat ACSTPstruct]=ACSTP(EEG, ACSTPoptions);
    %     ACSTPshow(EEG,ACSTPoptions,ACSTPstruct,1);
    %     set(gcf,'color',[1 1 1],'paperposition',[0 0 35 25])
    %     print(gcf, [outputDIR 'fig\ACSTP_' num2str(tix)],'-dtiff','-r450')
    %     ACSTPshow(EEG,ACSTPoptions,ACSTPstruct,3);
    %     set(gcf,'color',[1 1 1],'paperposition',[0 0 35 25])
    %     print(gcf, [outputDIR 'fig\ACSTPcomp_' num2str(tix)],'-dtiff','-r450')
    
    %% define global parameters
    close all
    
    % compute cospectra (from cross-spectra)
    Cf=[];%zeros(N,N);
    for indN1=1:N
        for indN2=1:N
            if indN1>=indN2 %upper triangle of the hermitian matrix
                %cross spectra computation (slow)
                [tmp Freqs]=cpsd(X(indN1,:)',X(indN2,:)',T,0,2*T,options.fs);
                %keep only wanted frequencies (set "options.freqs)
                indicesF=options.freqs(1)<Freqs & options.freqs(2)>Freqs;
                % CRITICAL STEP : convert cross-spectra to another SOS metric
                % real part : phase-aligned coherency 
                % imaginary part : lagged coherency
                % angle : lag
                % absolute : coherence
                Cf(indN1,indN2,:)=real(tmp(indicesF));
                if indN1~=indN2
                    %build lower triangle of the hermitian matrix
                    Cf(indN2,indN1,:)=conj(Cf(indN1,indN2,:));
                end
                
            end
            
        end
    end
    for indtmp=1:37;tmptmp(indtmp)=distance_riemann(Cf(:,:,indtmp),diag(diag(Cf(:,:,indtmp))));
    tmptmp2(indtmp)=norm(Cf(:,:,indtmp)-diag(diag(Cf(:,:,indtmp))),'fro')^2/norm(diag(Cf(:,:,indtmp)),'fro')^2;
    end
    figure;plot([tmptmp;tmptmp2]');legend('lol','lol2');
    max(tmptmp)/min(tmptmp)
    max(tmptmp2)/min(tmptmp2)
end
    
    Ct=zeros(N,N,K);
    for indK=1:K
        Ct(:,:,indK)=cov(X(:,:,indK)');
    end
    %% Define the target matrices and compute AJD, BAJD and CAJD and save
    weightTA=6;
    for six=Xmethods %define Xbar
        if six==1 % target matrices are the AEA (TA, NT) of Xk
            Xbar=cat(3,mean(X(:,:,Y==1),3)*weightTA,mean(X(:,:,Y==0),3));
        elseif six==2 % target matrices are Xk
            Xbar=X;
        else % target matrices are AEA (TA, NT) of Xk with bootstrap
            Xbar=gp_mean_bootstrap(X,[1 40],50,Y);
        end
        
        % [Xbar Ybar] = gp_mean_bootstrap(X,[5 10],20,Y);
        Cft = convertDat_BAJD(Xbar,Cf);
        
        
        %% perform Test
        disp(tix);
        for eix=InitMethods %define initiatilization
            
            %% define BSS parameters
            if eix==1
                B0      = eye(N);
                D0      = orth(randn(T,N));
            else
                B0 =  ACSTPstruct.fullAs{2}(1:N,1:N);
                D0 =  ACSTPstruct.fullAt{2}(1:T,1:N);
            end
            epsilon = 1e-18;
            itMax   = 1000;
            
            %% perform GPT-CAJD
            opCAJD = struct('B0',B0,'D0',D0,'eps',epsilon,'itMax',itMax);
            [B_CAJD,D_CAJD,S_CAJD,C_CAJD,info_CAJD] = gp_CAJD_GPT(Xbar,Cf,opCAJD);
            
            %% perform GPT-AJD
            opAJD = struct('B0',B0,'eps',epsilon,'itMax',itMax);
            [B_AJD,C_AJD,info_AJD] = gp_AJD_GPT(Cft,opAJD);
            
            %% perform GPT-BAJD
            opBAJD = struct('B0',B0,'D0',D0,'eps',epsilon,'itMax',itMax);
            [B_BAJD,D_BAJD,S_BAJD,info_BAJD] = gp_BAJD_GPT(Xbar,opBAJD);
            
            %% save results
            %             results(tix,eix,six) = struct('dat',[],'CAJD',[],'AJD',[],'BAJD',[],'ACSTP',[],'EEG',[]);
            results(tix,eix,six).dat = struct('X',X,'Cf',Cf,'Xbar',Xbar);
            results(tix,eix,six).CAJD = struct('C',C_CAJD,'S',S_CAJD,'conv',[info_CAJD.conv],'B',B_CAJD,'D',D_CAJD,'critOFFC',[info_CAJD.critOFFC],'critOFFS',[info_CAJD.critOFFS],'critA', [info_CAJD.critA],'critE', [info_CAJD.critE]);
            results(tix,eix,six).AJD = struct('C',C_AJD,'conv',[info_AJD.conv],'B',B_AJD,'critOFF',[info_AJD.critOFF],'critA',[info_AJD.critA]);
            results(tix,eix,six).BAJD = struct('S',S_BAJD,'conv',[info_BAJD.conv],'B',B_BAJD,'D',D_BAJD,'critOFF',[info_BAJD.critOFF],'critA', [info_BAJD.critA],'critE', [info_BAJD.critE]);
            
            results(tix,eix,six).ACSTP=ACSTPstruct;
            results(tix,eix,six).EEG=EEG;
        end
    end
end
toc
save([outputDIR 'CAJDreal'],'results')

%% (2) CRITERIA plot ALL Off-Diag Criteria for AJD/CAJD and BAJD/CAJD. Good to compare different initialization and target matrices

initMethods={'$\mathbf{B}_0=\mathbf{I}_N$, \mathbf{D}_0=randn$','$\mathbf{B}_0=\mathbf{B}_{acstp}, \mathbf{D}_0=\mathbf{D}_{acstp}$'};
XMethods={'$$\bar{\mathbf{X}}_{TA},\bar{\mathbf{X}}_{TA}$$','$\mathbf{X}_k$','$bootstrap \bar{\mathbf{X}}$'};
close all
% plot the convergence criteria for the spatial constract function, i.e.
% related to the diagonalization of C (only). available for AJD and CAJD
figure;
subplot1(2,3,'Gap',[0 0])
FontSize=12
subjects=12
AXIS=[0 100 -20 1]
for eix=1:2 %row
    for six=1:3 % columns
        for tix=subjects
            subplot1([eix six])
            if ~isempty(results(tix,eix,six).CAJD)
                plot(log10(results(tix,eix,six).CAJD.critOFFC),'b');
                hold on;
                plot(log10(results(tix,eix,six).AJD.critOFF),'r');
            end
            axis(AXIS)
            if six==1,ylabel('$conv$ (dB)','interpreter','latex') ;end
            if eix==1,title(XMethods{six},'interpreter','latex') ;end
            if eix==2,xlabel('$niter$','interpreter','latex') ;end
            
        end
        hold off
    end
end
Legends={'CAJD','AJD','BAJD'};%{'NoBAD','NoLAD','NoBAD bilinear'};
legend(Legends,'location','best');
set(gcf,'color',[1 1 1],'paperposition',[0 0 45 25])
set(gca,'fontsize',18)
print(gcf, [ outputDIR 'convergenceB' num2str(tix) ],'-dpng','-r450')
% plot the convergence criteria for the spatio-temporal constract function, i.e.
% related to the diagonalization of S (only). available for BAJD and CAJD
figure
subplot1(2,3,'Gap',[0 0])
AXIS=[0 100 -3 4]

for eix=1:2 %row
    for six=1:3 % columns
        for tix=subjects
            subplot1([eix six])
            if ~isempty(results(tix,eix,six).CAJD)
                plot(log10(results(tix,eix,six).CAJD.critOFFS),'b');
                hold on;
                plot(log10(results(tix,eix,six).BAJD.critOFF),'r');
                
            end
            axis(AXIS)
            if six==1,ylabel('$conv$ (dB)','interpreter','latex') ;end
            if eix==1,title(XMethods{six},'interpreter','latex') ;end
            if eix==2,xlabel('$niter$','interpreter','latex') ;end
            
        end
        hold off
    end
end
Legends={'CAJD','BAJD'};%{'NoBAD','NoLAD','NoBAD bilinear'};
legend(Legends,'location','best');
set(gcf,'color',[1 1 1],'paperposition',[0 0 45 25])
set(gca,'fontsize',18)
print(gcf, [ outputDIR 'convergenceD' num2str(tix) ],'-dpng','-r450')

%% (3) CRITERIA select specific Off-Diag Criteria for AJD/CAJD and BAJD/CAJD. For good figures
initMethods={'$\mathbf{B}_0=\mathbf{I}_N$, \mathbf{D}_0=randn$','$\mathbf{B}_0=\mathbf{B}_{acstp}, \mathbf{D}_0=\mathbf{D}_{acstp}$'};
XMethods={'$$\bar{\mathbf{X}}_{TA},\bar{\mathbf{X}}_{TA}$$','$\mathbf{X}_k$','$bootstrap \bar{\mathbf{X}}$'};
close all
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
% plot the convergence criteria for the spatial constract function, i.e.
% related to the diagonalization of C (only). available for AJD and CAJD
figure;
subplot1(1,2,'Gap',[0.01 0],'YTickL','All')
FontSize=12
subjects=[1:19]
XLIM=[0 99]
YLIM=[-8 1]
Init2plot=2
XMethod2plot=3
FontSize=12
Paper=[0 0 20 12]
LineWidth=1.5
MarkersSize=6
for eix=Init2plot %row
    for six=XMethod2plot % columns
        for tix=subjects
            subplot1([1])
            if ~isempty(results(tix,eix,six).CAJD)
                plot(log10(results(tix,eix,six).CAJD.critOFFC),'color',Colors{2},'Linewidth',LineWidth,'Markers',MarkersSize);
                hold on;
                plot(log10(results(tix,eix,six).AJD.critOFF),'color',Colors{1},'Linewidth',LineWidth,'Markers',MarkersSize);
            end
            xlim(XLIM);            ylim(YLIM)
            
            ylabel('$I_{n-d}~~(\mathbf{B})$ (dB)','interpreter','latex') ;
            title('(a)','interpreter','latex') ;
            xlabel('$niter$','interpreter','latex') ;
            
        end
        hold off
    end
end
Legends={'CAJD','AJD'};%{'NoBAD','NoLAD','NoBAD bilinear'};
legend(Legends,'location','southwest');
% plot the convergence criteria for the spatio-temporal constract function, i.e.
% related to the diagonalization of S (only). available for BAJD and CAJD
YLIM=[ -3 4]
for eix=Init2plot %row
    for six=XMethod2plot % columns
        for tix=subjects
            subplot1([2])
            if ~isempty(results(tix,eix,six).CAJD)
                plot(log10(results(tix,eix,six).CAJD.critOFFS),'color',Colors{2},'Linewidth',LineWidth,'Markers',MarkersSize);
                hold on;
                plot(log10(results(tix,eix,six).BAJD.critOFF),'color',Colors{3},'Linewidth',LineWidth,'Markers',MarkersSize);
                
            end
            xlim(XLIM);            ylim(YLIM)
            if six==1,ylabel('$conv$ (dB)','interpreter','latex') ;end
            title('(b)','interpreter','latex') ;
            if eix==2,xlabel('$niter$','interpreter','latex') ;end
            ylabel('$I_{n-d}~~(\mathbf{D})$ (dB)','interpreter','latex') ;
            
        end
        hold off
    end
end
set(gca,'yaxislocation','right');
Legends={'CAJD','BAJD'};%{'NoBAD','NoLAD','NoBAD bilinear'};
legend(Legends,'location','southwest');
set(gcf,'color',[1 1 1],'paperposition',Paper)
if ~(length(subjects)==1)
    print(gcf, [ outputDIR 'critofs'],'-dpng','-r450')
else
    print(gcf, [ outputDIR 'critof' num2str(tix) ],'-dpng','-r450')
end

%% (4) CRITERIA histogram for all the SNR
% load all parameters
close all
clear y errY Table TableSTD TableQuMin TableQuMax tmpTable
allN=[1 2]
allHist=allN
stdfactor=1;
FontSize=16;
PaperFactor=1
Paper=[0 0 20 6]
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
Legends={'AJD','CAJD','BAJD'};
XLegend1={'(a)','(b)'}%'I_{n-d}(\mathbf{B})' 'I_{n-d}(\mathbf{B},\mathbf{D})'}%Legends;%{'','40','20','16','12','8','4'}';
MinMax=[0 length(XLegend1)+1 -45 20];
% Colors=cellfun(@(x) x/10-.1,{[5 10 10],[10. 3. 3.],[3.33 3.33 3]},'Uniform',0)
Xlabels='$n$'
% XLegend2={'','200','100','80','60','40','20'}';
pval_err=0.2
yLegends={'$I_{{n-d}}$ ~(dB)'}
outputDIR='D:\Mes Documents GIPSA\MATLAB\figures\NoBAD_eusipco\'
outputDIR2='D:\Mes Documents GIPSA\Présentations & Communications\2016-01 BSS versus BAJD ERP\Figures\'

 outputFIG=[outputDIR2 'CritOfquantile']
% outputFIG=[outputDIR 'CritOfSTD']

subjects=1:19;
Init2plot=2;
XMethod2plot=3;
Histos=1:2
% prepare the data such as
% y : are the columns of the histogram
% errY : are the std of the histogram
% with CrOf and MoAm:
% 1method      3SNR                      7trials
% 3       1     6      1      1     1     600
subplot1(1,2)
iter2plot=300
for indRow=1
    for indHisto=Histos%:length(allHist)
        %critical part, the input for the histogram
        clear tmpTable
        for tix=subjects
            if indHisto==1
                            tmpTable(2,tix)=10*log10(squeeze([results(tix,Init2plot,XMethod2plot).CAJD.critOFFC(iter2plot)])); %make such as [nb methods x nb trials]

                tmpTable(1,tix)=10*log10(squeeze([results(tix,Init2plot,XMethod2plot).AJD.critOFF(iter2plot)])); %make such as [nb methods x nb trials]
                tmpTable(3,tix)=NaN; %make such as [nb methods x nb trials]
                             [~,Tpval(indRow,indHisto),~,statsT(indRow,indHisto)]=ttest2(tmpTable(1,:),tmpTable(2,:),0.05);
                
            else
                tmpTable(1,tix)=NaN; %make such as [nb methods x nb trials]
                            tmpTable(2,tix)=10*log10(squeeze([results(tix,Init2plot,XMethod2plot).CAJD.critOFFS(iter2plot)])); %make such as [nb methods x nb trials]

                tmpTable(3,tix)=10*log10(squeeze([results(tix,Init2plot,XMethod2plot).BAJD.critOFF(iter2plot)])); %make such as [nb methods x nb trials]
                                             [~,Tpval(indRow,indHisto),~,statsT(indRow,indHisto)]=ttest2(tmpTable(3,:),tmpTable(2,:),0.05);

            end
        end
        
        Table(indHisto,:,indRow)=median(tmpTable,2);%stack the observations
        %         NT(indHisto,:,indRow)=mean(mean(tmpNT1,2),1);
        %         TAcstp(indHisto,indMethod,indRow)=mean(mean(tmpTA2,2),1); % now
        %         included in the previous ligne
        %         NTcstp(indHisto,indMethod,indRow)=mean(mean(tmpNT2,2),1);
        %std
        TableSTD(indHisto,:,indRow)=std(tmpTable,[],2)*stdfactor;
        TableQuMin(indHisto,:,indRow)=quantile(tmpTable,pval_err/2,2);
        TableQuMax(indHisto,:,indRow)=quantile(tmpTable,1-pval_err/2,2);
        % need to clean up a bit for the statistics (here comparison
        % between only two methods)
        %         for all the methods
        %             [~,NTpval(indRow,indHisto)]=ttest(tmpNT1(:),tmpNT2(:),0.05);
        %             TApval(TApval>0.05)=nan;
        %             NTpval(NTpval>0.05)=nan;
        %         end
    end
end

disp('preparation DONE')
%xticks = linspace(0,7, numel(XLegend1));
MaxRow=1;MaxColumn=1;
for indFigure=1:1
    close all
    hrmse=figure;
    subplot1(MaxRow,MaxColumn,'Gap',[0 0.005]); %change for [Max Row, max Column]
    
    
    ind=1;
    for indRow=(1:MaxRow)+(indFigure-1)*8
        %         for indHisto=Histos
        clear errY
        %         subplot1(indHisto)
        y=[Table(:,:,indRow)];
        %bar([NT(1,:)' NTcstp(1,:)'])
%         errY(:,:,1) = [TableSTD(:,:,indRow)];   % the std
        % OR
                   errY(:,:,1) = [TableQuMin(:,:,indRow)-y]; errY(:,:,2) = [TableQuMax(:,:,indRow)-y]; %the quantile
        
        %          remplacementVal=(MinMax(end)-errY*stdfactor)*0.5;
        %          y((errY*stdfactor+y)>=MinMax(end)*.8)=remplacementVal((errY*stdfactor+y)>=MinMax(end)*0.8);
        hbar=barwitherr(errY, y);    % Plot with errorbars
        axis(MinMax)
        %         sigstar({[0.8,1.2], [0.8,1.2]+1,[0.8,1.2]+2,[0.8,1.2]+3,[0.8,1.2]+4,[0.8,1.2]+5},TApval(indRow,:))
        for indMethods=1:length(hbar)
            set(hbar(indMethods),'FaceColor',Colors{indMethods});
        end
        ylabel(yLegends{indRow},'FontSize',FontSize,'fontname','times new roman','FontAngle','italic','interpreter','latex');
        
        if ind==MaxRow
            set(gca,'xtick',1:length(XLegend1),'xticklabel',XLegend1,'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
            rotateXLabels(gca,0)
            %             xlabel(Xlabels,'interpreter','latex','fontsize',FontSize)
            %             xlabh = get(gca,'XLabel');
            % set(xlabh,'Position',get(xlabh,'Position')-[1 10 1])
            legend(Legends,'location','southeast')
        end
        ind=ind+1;
    end
end
set(hrmse, 'PaperPosition', Paper*PaperFactor,'units','normalized','outerposition',[0.45 0.1 0.5 .9])

set(gcf, 'color', [1 1 1])
% colormap('gray')
%     annotation('textbox', [0.05,0.001,0.5,.05],'units','normalized',...
%         'String',['CrOf - std/' num2str(1/stdfactor)],'fontsize',(FontSize),'fontname','arial new roman','EdgeColor','none','fontangle','italic');
%     annotation('textbox', [0.05,0.001,0.5,.025],'units','normalized',...
%         'String',['°p-val<0.05 *p-val<0.01'],'fontsize',(FontSize),'fontname','arial new roman','EdgeColor','none','fontangle','italic');

print(hrmse,[outputFIG],'-dpng','-r450')
print(hrmse,[outputFIG],'-dpng','-r450')

%% (5) CRITERIA select specific convergence for AJD/CAJD and BAJD/CAJD. For good figures
initMethods={'$\mathbf{B}_0=\mathbf{I}_N$, \mathbf{D}_0=randn$','$\mathbf{B}_0=\mathbf{B}_{acstp}, \mathbf{D}_0=\mathbf{D}_{acstp}$'};
XMethods={'$$\bar{\mathbf{X}}_{TA},\bar{\mathbf{X}}_{TA}$$','$\mathbf{X}_k$','$bootstrap \bar{\mathbf{X}}$'};
close all
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
% plot the convergence criteria for the spatial constract function, i.e.
% related to the diagonalization of C (only). available for AJD and CAJD
figure;
FontSize=12
subjects=[1 6 11]
XLIM=[0 999]
YLIM=[-8 5]
Init2plot=2
XMethod2plot=3
FontSize=12
Paper=[0 0 15 10]
LineWidth=1.5
MarkersSize=6

for eix=Init2plot %row
    for six=XMethod2plot % columns
        for tix=subjects
            if ~isempty(results(tix,eix,six).CAJD)
                plot(log10(results(tix,eix,six).CAJD.conv),'color',Colors{2},'Linewidth',LineWidth,'Markers',MarkersSize);
                hold on;
                plot(log10(results(tix,eix,six).AJD.conv),'color',Colors{1},'Linewidth',LineWidth,'Markers',MarkersSize);
                plot(log10(results(tix,eix,six).BAJD.conv),'color',Colors{3},'Linewidth',LineWidth,'Markers',MarkersSize);
                
            end
            xlim(XLIM);            ylim(YLIM)
            
            ylabel('$conv$ (dB)','interpreter','latex') ;
            title(XMethods{six},'interpreter','latex') ;
            xlabel('$niter$','interpreter','latex') ;
            
        end
        hold off
    end
end
Legends={'CAJD','AJD','BAJD'};%{'NoBAD','NoLAD','NoBAD bilinear'};
legend(Legends,'location','southwest');
% plot the convergence criteria for the spatio-temporal constract function, i.e.
% related to the diagonalization of S (only). available for BAJD and CAJD
set(gcf,'color',[1 1 1],'paperposition',Paper)
if ~(length(subjects)==1)
    print(gcf, [ outputDIR 'Convergences'],'-dpng','-r450')
else
    print(gcf, [ outputDIR 'Convergence' num2str(tix) ],'-dpng','-r450')
end

%% (6) SOURCES figure Topoplot for a given subject CAJD
indInit=2
indX=3
close all

for subject=2

if subject<10
    figDIR=([outputDIR 'fig\ekate0' num2str(subject) '\']);
else
    figDIR=([outputDIR 'fig\ekate' num2str(subject) '\']);
end
mkdir(figDIR)
X=ALLdata.Xte{subject}(BImap{1},:,:);
Y=ALLdata.Yte{subject};
ACSTPstruct=results(subject,indInit,indX).ACSTP;
EEG=results(subject,indInit,indX).EEG;
D=results(subject,indInit,indX).CAJD.D;
B=results(subject,indInit,indX).CAJD.B;
K=size(X,3);
N=size(X,1);
T=size(X,2);
S_hat=[];
[~, tmpo]=findElectrodes(EEG.ElectrodesName, EEG.ElectrodesName);
clear localization
localization=tmpo;

for indK=1:K
    S_hat(:,:,indK)=B'*X(:,:,indK);
end

% plot temporal observation of X (mean, mean acstp, mean filtered, pinv(D))
figure
subplot1(1,5,'Gap',[0.02 0.01])
subplot1(1);plotEEG(mean(X(:,:,Y==1),3),4,EEG.Fs,EEG.ElectrodesName); title('${\bar{X}}_{TA}$','interpreter','latex')
subplot1(2);plotEEG((ACSTPstruct.fullAt{2}(:,1:16))');  title('$D_{ACSTP}$','interpreter','latex');set(gca,'fontsize',12)
subplot1(3);plotEEG(normEEG(mean(S_hat(:,:,Y==1),3)),2); title('$B^T {\bar{X}}_{TA}$','interpreter','latex');set(gca,'Yticklabel',[])
subplot1(4);plotEEG(normEEG(mean(S_hat(:,:,Y==0),3)),2); title('$B^T {\bar{X}}_{NT}$','interpreter','latex');set(gca,'Yticklabel',[])
subplot1(5);plotEEG(normEEG(pinv(D))); title('$D^{+}$','interpreter','latex');set(gca,'Yticklabel',[])

set(gcf,'color',[1 1 1],'paperposition',[0 0 45 25])
set(gcf,'position',[0.05 0.1 0.45 0.8],'unit','normalized');    set(gcf,'position',[0.05 0.1 0.45 0.8],'unit','normalized')



print(gcf, [figDIR 'S' num2str(subject) '_averages_X' num2str(indX) ],'-dtiff','-r450')

% topoplot of the CAJD components
figure
figure('name','CAJD')
gp_source_topoplot([ figDIR 'S' num2str(subject) 'comp_X' num2str(indX) '_CAJD'],EEG,results(subject,indInit,indX).CAJD.B,results(subject,indInit,indX).CAJD.D)

% topoplot of the BAJD components
figure('name','BAJD')
gp_source_topoplot([ figDIR 'S' num2str(subject) 'comp_X' num2str(indX) '_BAJD'],EEG,results(subject,indInit,indX).BAJD.B,results(subject,indInit,indX).BAJD.D)

% topoplot of the AJD components
figure('name','AJD')
gp_source_topoplot([ figDIR 'S' num2str(subject) 'comp_X' num2str(indX) '_AJD'],EEG,results(subject,indInit,indX).AJD.B)

% plot spatial covariances of the sources
figure('name','Mean Covariance Comparison')
subplot(131)
imagesc((abs(normEEG(X(:,:))*normEEG(X(:,:))' )));title('Ck')
subplot(132)
imagesc((abs(normEEG((results(subject,indInit,indX).AJD.B'*X(:,:)))*normEEG((results(subject,indInit,indX).AJD.B'*X(:,:)))')))
;title('B^T Xk B (AJD)')
subplot(133)
imagesc((abs(normEEG((results(subject,indInit,indX).CAJD.B'*X(:,:)))*normEEG((results(subject,indInit,indX).CAJD.B'*X(:,:)))')))
title('B^T Xk B (CAJD)')
set(gcf,'color',[1 1 1],'paperposition',[0 0 45 15])
print(gcf, [ figDIR 'S' num2str(subject) 'autocorr_X' num2str(indX)],'-dpng','-r450')
end

%% (7) SOURCES compute the source projection into the electrodes domain
% close all
close all
scale=3
subject=1
indInit=2
indX=3
clear Shatk1 Shatk0 normF
X=ALLdata.Xte{subject}(BImap{1},:,:);
Y=ALLdata.Yte{subject};
ACSTPstruct=results(subject,indInit,indX).ACSTP;
EEG=results(subject,indInit,indX).EEG;
D=results(subject,indInit,indX).CAJD.D;
B=results(subject,indInit,indX).CAJD.B;

Y1=(Y==1);
X1=X(:,:,Y1);
Y0=(Y==0);
X0=X(:,:,Y0);

CompCAJD=13
CompBAJD=7
CompAJD=16

%          Ks=randi(100,1,5)
[Xhat_cajd1 Xhat_cajd1b]=gp_source_projection(X1,CompCAJD,results(subject,indInit,indX).CAJD.B,results(subject,indInit,indX).CAJD.D);
[Xhat_bajd1 Xhat_bajd1b]=gp_source_projection(X1,CompBAJD,results(subject,indInit,indX).BAJD.B,results(subject,indInit,indX).BAJD.D);
[Xhat_ajd1]=gp_source_projection(X1,CompAJD,results(subject,indInit,indX).AJD.B);

[Xhat_cajd0 Xhat_cajd0b]=gp_source_projection(X0,CompCAJD,results(subject,indInit,indX).CAJD.B,results(subject,indInit,indX).CAJD.D);
[Xhat_bajd0 Xhat_bajd0b]=gp_source_projection(X0,CompBAJD,results(subject,indInit,indX).BAJD.B,results(subject,indInit,indX).BAJD.D);
[Xhat_ajd0]=gp_source_projection(X0,CompAJD,results(subject,indInit,indX).AJD.B);

%% (8a) SOURCES plot and save the plot of the projected sources (7required)
scale=3

    close all
    figure
    subplot1(1,3)
    Ks=[ 65,61,28,86,56,74,85,64,32,67];
    for indK=Ks%size(X,3)
        subplot1(1)
        plotEEG(X1(:,:,indK)/5,scale,EEG.Fs,EEG.ElectrodesName,[],[0.3 0.3 0.3]);
                hold on

        subplot1(2)
        plotEEG(Xhat_cajd1(:,:,indK),scale,EEG.Fs,EEG.ElectrodesName,[],[0.3 0.3 0.3]);
                hold on

        subplot1(3)
        plotEEG(Xhat_bajd1(:,:,indK),scale,EEG.Fs,EEG.ElectrodesName,[],[0.3 0.3 0.3]);
        hold on
        title('TA')
        set(gca,'xticklabel',[])
%         subplot1(3)
%         plotEEG((Shatk0(:,:,indK)),scale,EEG.Fs,EEG.ElectrodesName);
%         hold on
%         title('NT')
%         set(gca,'yticklabel',[])
        
        
    end
        subplot1(1);plotEEG(mean(X1,3)/5,scale,EEG.Fs,EEG.ElectrodesName,[],[1 0 0],[],[],2)
%         x=[0.3,.3];
% y=[-scale,scale];
% plot(x,y,'--')
    subplot1(2);plotEEG(mean(Xhat_cajd1,3),scale,EEG.Fs,EEG.ElectrodesName,[],[1 0 0],[],[],2)
%     plot(x,y,'--')
    subplot1(3);plotEEG(mean(Xhat_bajd1,3),scale,EEG.Fs,EEG.ElectrodesName,[],[1 0 0],[],[],2)

%     subplot1(3);plotEEG(mean(Shatk0,3),scale,EEG.Fs,EEG.ElectrodesName,[],[0 0 1],[],[],2)
%     plot(x,y,'--')

    set(gca,'yticklabel',[])
    
    dim = [.05 .005 .5 .1];
    str=['reconstructed source at the single trial level of the source=' num2str(CompCAJD) ' for NT and NT. Butterfly plot for K=' num2str(length(Ks))]
    set(gcf,'color',[1 1 1],'paperposition',[0 0 20 25])
    annotation('textbox',dim,'String',str,'unit','normalized','fontsize',12,'linewidth',0,'edgecolor',[1 1 1]);
    set(gcf,'position',[0.5 0.1 0.4 0.8],'unit','normalized');    set(gcf,'position',[0.5 0.1 0.4 0.8],'unit','normalized')
    
    print(gcf, [figDIR 'S' num2str(subject) '_X' num2str(indX)  'proj_s' num2str(CompCAJD)],'-dpng','-r450')
    
    % subplot(121); plot(mean(X1(7,:,:),3)','k','linewidth',1.5);
    % subplot(122); plotEEG((mean(Shat(Compo,:,:),3)),'k','linewidth',1.5);
    
    hold off

%% (8b) SOURCES TA variability with BAJD model  (7required)
close all
FontSize=14
% Zbarz=Shatk1(:,:,Ks)
scale=6
allChan=[4,7,11]
yticklabels = -scale/2:scale:scale/2;
yticks = linspace(-(1.01)*(scale), (1.01)*(scale), numel(yticklabels));
Fs=128;
offset=0;
xticklabels = (offset/Fs):.25:(Fs+offset)/Fs;
xticks = linspace(0, (1.01)*(Fs)/Fs, numel(xticklabels));

subplot1(length(allChan),3)
for indC=1:length(allChan)
    subplot1([indC 1])
    [h1 h2 h3]=plotEEGvariability(X1,'channel',allChan(indC))
    if indC==1 title('${\mathbf{X}_k}$ TA','FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale*3 scale*3]);
    ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    set(get(gca,'YLabel'),'Rotation',0) %put the ylabel horizontal
    h = get(gca,'ylabel');
    pos = get(h,'position');
    xlimits = get(gca,'xlim');
    pos(1) = xlimits(1) -0.1 * (xlimits(2) - xlimits(1));
    set(h,'position',pos)
    
end

for indC=1:length(allChan)
    subplot1([indC 2])
    [h1 h2 h3]=plotEEGvariability(Xhat_cajd1,'channel',allChan(indC))
    if indC==1 title(['${\hat{\mathbf{X}}_{k,' num2str(CompCAJD) '}}$ ~CAJD TA'],'FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale scale]);
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    
end

for indC=1:length(allChan)
    subplot1([indC 3])
    [h1 h2 h3]=plotEEGvariability(Xhat_bajd1,'channel',allChan(indC))
    if indC==1 title(['${\hat{\mathbf{X}}_{k,' num2str(CompBAJD) '}}$ ~BAJD TA'],'FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale scale]);
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    
end

print(gcf, [figDIR 'S' num2str(subject) '_X' num2str(indX)  'proj_s' num2str(CompCAJD) 'AEA TA'],'-dpng','-r450')

%% (8c) SOURCES TA variability with AJD model  (7required)
close all
FontSize=14
% Zbarz=Shatk1(:,:,Ks)
scale=10
allChan=[4,7,11]
yticklabels = -scale/2:scale:scale/2;
yticks = linspace(-(1.01)*(scale), (1.01)*(scale), numel(yticklabels));
Fs=128;
offset=0;
xticklabels = (offset/Fs):.25:(Fs+offset)/Fs;
xticks = linspace(0, (1.01)*(Fs)/Fs, numel(xticklabels));

subplot1(length(allChan),4)
for indC=1:length(allChan)
    subplot1([indC 1])
    [h1 h2 h3]=plotEEGvariability(X1,'channel',allChan(indC))
    if indC==1 title('${\mathbf{X}_k}$ TA','FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale*3 scale*3]);
    ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    set(get(gca,'YLabel'),'Rotation',0) %put the ylabel horizontal
    h = get(gca,'ylabel');
    pos = get(h,'position');
    xlimits = get(gca,'xlim');
    pos(1) = xlimits(1) -0.1 * (xlimits(2) - xlimits(1));
    set(h,'position',pos)
    
end

for indC=1:length(allChan)
    subplot1([indC 2])
    [h1 h2 h3]=plotEEGvariability(Xhat_cajd1b,'channel',allChan(indC))
    if indC==1 title(['${\hat{\mathbf{X}}_{k,' num2str(CompCAJD) '}}$ ~CAJD TA'],'FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale scale]);
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    
end

for indC=1:length(allChan)
    subplot1([indC 3])
    [h1 h2 h3]=plotEEGvariability(Xhat_bajd1b,'channel',allChan(indC))
    if indC==1 title(['${\hat{\mathbf{X}}_{k,' num2str(CompBAJD) '}}$ ~BAJD TA'],'FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale scale]);
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    
end

for indC=1:length(allChan)
    subplot1([indC 4])
    [h1 h2 h3]=plotEEGvariability(Xhat_ajd1,'channel',allChan(indC))
    if indC==1 title(['${\hat{\mathbf{X}}_{k,' num2str(CompAJD) '}}$ ~AJD TA'],'FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale scale]);
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    
end

print(gcf, [figDIR 'S' num2str(subject) '_X' num2str(indX)  'proj_s' num2str(CompCAJD) 'AEA TA 2'],'-dpng','-r450')


%% (9) SOURCES NT variability with BAJD model  (7required)
close all
FontSize=14
% Zbarz=Shatk1(:,:,Ks)
allChan=[4,7,11]
yticklabels = -scale/2:scale:scale/2;
yticks = linspace(-(1.01)*(scale), (1.01)*(scale), numel(yticklabels));
Fs=128;
offset=0;
xticklabels = (offset/Fs):.25:(Fs+offset)/Fs;
xticks = linspace(0, (1.01)*(Fs)/Fs, numel(xticklabels));

subplot1(length(allChan),2)
for indC=1:length(allChan)
    subplot1([indC 1])
    [h1 h2 h3]=plotEEGvariability(X0,'channel',allChan(indC))
    if indC==1 title('${\bar{X}}$ NT','FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale*5 scale*5]);
    ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    set(get(gca,'YLabel'),'Rotation',0) %put the ylabel horizontal
    h = get(gca,'ylabel');
    pos = get(h,'position');
    xlimits = get(gca,'xlim');
    pos(1) = xlimits(1) -0.1 * (xlimits(2) - xlimits(1));
    set(h,'position',pos)
    
end

for indC=1:length(allChan)
    subplot1([indC 2])
    [h1 h2 h3]=plotEEGvariability(Shatk0,'channel',allChan(indC))
    if indC==1 title('${\bar{X}}$ NT','FontSize',FontSize,'fontname','times new roman','Fontangle','italic','interpreter','latex');end
    axis([0 1 -scale scale]);
    set(gca,'Ytick',[])
    %     set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
    if indC==length(allChan) % if last line, put the xlabel
        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
    end
    
end

print(gcf, [figDIR 'S' num2str(subject) '_X' num2str(indX)  'proj_s' num2str(CompCAJD) 'AEA NT'],'-dpng','-r450')