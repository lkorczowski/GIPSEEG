%% RETROJECTION SPATIALE
close all
Type='AJD';
tix=1;%subject
EEG=results(tix,eix,six).EEG;
 X1=EEG.Epoch(:,:,EEG.EpochClass==1);X1bar=mean(X1,3);
  X0=EEG.Epoch(:,:,EEG.EpochClass==0);X0bar=mean(X0,3);
B=results(tix,eix,six).AJD.B;
A=pinv(B);
isARTE=eval(['results(tix,eix,six).' Type '.SACC_all']);

   Atilde1=A;
            Atilde1(1:end & (isARTE~=0),:)=0;
            Atilde2=A;
            Atilde2( ~(1:end & (isARTE~=0)),:)=0;
            X1bar_cleaned_ajd=Atilde1'*B'*Xbar1;

for indS=10
    Atilde1=A;
    Atilde1(~((1:end)==indS),:)=0;
    % clean setup
    %hold off;
    figure('name',[Type 'VAR' num2str(indS) ],'visible','off');subplot1(length(allChan),1,'Gap',[0.001, 0.000]);
    
    [Xhat1ajd]=gp_source_projection(X1,indS,B); %backprojection
    
    for indC=1:length(allChan)
        subplot1([indC 1])
        hold off
        plotEEGvariability(Xhat1ajd,'channel',allChan(indC),'fs',EEG.Fs,'offset',EEG.offset,'alpha',VARALPHA,'type',VARTYPE);
        if indC==1; title(Type,'FontSize',FontSize*1.33,'fontname','times new roman',...
                'Fontangle','italic','interpreter','latex');end
        axis([EEG.offset/EEG.Fs (size(EEG.Epoch,2)+EEG.offset-2)/EEG.Fs -scale scale]);
        ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
        set(gca,'Ytick',[])
        %set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
        if indC==length(allChan) % if last line, put the xlabel
            xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
            set(gca, 'XTick', xticks,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
        end
        set(gca,'YLabel',[]) %put the ylabel horizontal
        l=vline(0.3,'r');l.LineWidth=0.75;
        vline(0,'g:');l.LineWidth=0.75;
         hold all;   plotEEG(X1bar(indC,:))

    end
end

[rmse_AJD]=RMSE2(Xhat1ajd);%,X1bar_cleaned_ajd);
mean(rmse_AJD)

Type='CAJD';
B=results(tix,eix,six).CAJD.B;
A=pinv(B);
isARTE=eval(['results(tix,eix,six).' Type '.SACC_all']);
   Atilde1=A;
            Atilde1(1:end & (isARTE~=0),:)=0;
            Atilde2=A;
            Atilde2( ~(1:end & (isARTE~=0)),:)=0;
            X1bar_cleaned_cajd=Atilde1'*B'*Xbar1;

for indS=14
    Atilde1=A;
    Atilde1(~((1:end)==indS),:)=0;
    % clean setup
    %hold off;
    figure('name',[Type 'VAR' num2str(indS) ],'visible','off');subplot1(length(allChan),1,'Gap',[0.001, 0.000]);
    
    [Xhat1cajd]=gp_source_projection(X1,indS,B); %backprojection
    
    for indC=1:length(allChan)
        subplot1([indC 1])
        hold off
        plotEEGvariability(Xhat1cajd,'channel',allChan(indC),'fs',EEG.Fs,'offset',EEG.offset,'alpha',VARALPHA,'type',VARTYPE);
        if indC==1; title(Type,'FontSize',FontSize*1.33,'fontname','times new roman',...
                'Fontangle','italic','interpreter','latex');end
        axis([EEG.offset/EEG.Fs (size(EEG.Epoch,2)+EEG.offset-2)/EEG.Fs -scale scale]);
        ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
        set(gca,'Ytick',[])
        %set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
        if indC==length(allChan) % if last line, put the xlabel
            xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
            set(gca, 'XTick', xticks,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
        end
        set(gca,'YLabel',[]) %put the ylabel horizontal
        l=vline(0.3,'r');l.LineWidth=0.75;
        vline(0,'g:');l.LineWidth=0.75;
    end
end
[rmse_CAJD]=RMSE2(Xhat1cajd);%,X1bar_cleaned_cajd);
mean(rmse_CAJD)
hhist=figure('visible','on');
histogram(rmse_AJD);hold all;histogram(rmse_CAJD);
legend('AJD','CAJD')
ylabel('population')
xlabel('RMSE')
print(hhist,'D:\GoogleDrive\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\images\6-RMSE_AJD-CAJD','-dpdf')

[h,p,ci,stats]=ttest2(rmse_AJD,rmse_CAJD);

%% RETROPROECTION D
close all
Type='BAJD';
B=results(tix,eix,six).BAJD.B;
D=results(tix,eix,six).BAJD.D;
scale=10
A=pinv(B);
isARTE=eval(['results(tix,eix,six).' Type '.SACC_all']);
   Atilde1=A;
            Atilde1(1:end & (isARTE~=0),:)=0;
            Atilde2=A;
            Atilde2( ~(1:end & (isARTE~=0)),:)=0;
            X1bar_cleaned_cajd=Atilde1'*B'*Xbar1;
            
            STA=-EEG.offset+(floor(0.2*EEG.Fs):floor(0.45*EEG.Fs));

for indS=10
    Atilde1=A;
    Atilde1(~((1:end)==indS),:)=0;
    % clean setup
    %hold off;
    figure('name',[Type 'VAR' num2str(indS) ],'visible','off');subplot1(length(allChan),1,'Gap',[0.001, 0.000]);
    
    [Xhat1bajd]=gp_source_projection(X1,indS,B,D); %backprojection
        [Xhat0bajd]=gp_source_projection(X0,indS,B,D); %backprojection
    for indC=1:length(allChan)
        subplot1([indC 1])
        hold off
        plotEEGvariability(Xhat1bajd,'channel',allChan(indC),'fs',EEG.Fs,'offset',EEG.offset,'alpha',VARALPHA,'type',VARTYPE);
        if indC==1; title(Type,'FontSize',FontSize*1.33,'fontname','times new roman',...
                'Fontangle','italic','interpreter','latex');end
        axis([EEG.offset/EEG.Fs (size(EEG.Epoch,2)+EEG.offset-2)/EEG.Fs -scale scale]);
        ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
        set(gca,'Ytick',[])
        %set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
        if indC==length(allChan) % if last line, put the xlabel
            xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
            set(gca, 'XTick', xticks,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
        end
        set(gca,'YLabel',[]) %put the ylabel horizontal
        l=vline(0.3,'r');l.LineWidth=0.75;
        vline(0,'g:');l.LineWidth=0.75;
    end
end

[rmse_BAJD]=RMSE2(Xhat1bajd)%,X1bar_cleaned_cajd);
amp_BAJD=squeeze(mean(mean(Xhat1bajd(:,STA,:),1),2));
amp_BAJD0=squeeze(mean(mean(Xhat0bajd(:,STA,:),1),2));

Type='CAJD';
B=results(tix,eix,six).CAJD.B;
D=results(tix,eix,six).CAJD.D;
A=pinv(B);
isARTE=eval(['results(tix,eix,six).' Type '.SACC_all']);
   Atilde1=A;
            Atilde1(1:end & (isARTE~=0),:)=0;
            Atilde2=A;
            Atilde2( ~(1:end & (isARTE~=0)),:)=0;
            X1bar_cleaned_cajd=Atilde1'*B'*Xbar1;

for indS=14
    Atilde1=A;
    Atilde1(~((1:end)==indS),:)=0;
    % clean setup
    %hold off;
    figure('name',[Type 'VAR' num2str(indS) ],'visible','off');subplot1(length(allChan),1,'Gap',[0.001, 0.000]);
    
    [Xhat1cajd]=gp_source_projection(X1,indS,B,D); %backprojection
            [Xhat0cajd]=gp_source_projection(X0,indS,B,D); %backprojection

    for indC=1:length(allChan)
        subplot1([indC 1])
        hold off
        plotEEGvariability(Xhat1cajd,'channel',allChan(indC),'fs',EEG.Fs,'offset',EEG.offset,'alpha',VARALPHA,'type',VARTYPE);
        if indC==1; title(Type,'FontSize',FontSize*1.33,'fontname','times new roman',...
                'Fontangle','italic','interpreter','latex');end
        axis([EEG.offset/EEG.Fs (size(EEG.Epoch,2)+EEG.offset-2)/EEG.Fs -scale scale]);
        ylabel(LABELS(allChan(indC)),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
        set(gca,'Ytick',[])
        %set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
        if indC==length(allChan) % if last line, put the xlabel
            xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
            set(gca, 'XTick', xticks,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
        end
        set(gca,'YLabel',[]) %put the ylabel horizontal
        l=vline(0.3,'r');l.LineWidth=0.75;
        vline(0,'g:');l.LineWidth=0.75;
    end
end
[rmse_CAJD]=RMSE2(Xhat1cajd)%,X1bar_cleaned_cajd);
amp_CAJD=squeeze(mean(mean(Xhat1cajd(:,STA,:),1),2));
amp_CAJD0=squeeze(mean(mean(Xhat0cajd(:,STA,:),1),2));


hhist=figure('visible','on');
subplot(211)
histogram(amp_BAJD);hold all;histogram(amp_BAJD0);
xlim([-2 2])
legend('BAJD-TA','BAJD-NT')
ylabel('population')

subplot(212)
histogram(amp_CAJD);;hold all;histogram(amp_CAJD0);
xlim([-2 2])

legend('CAJD-TA','CAJD-NT')
ylabel('population')
xlabel('amplitude moyenne')
print(hhist,'D:\GoogleDrive\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\images\6-RMSE_BAJD-CAJD','-dpdf')

[h,p,ci,stats]=ttest2(amp_BAJD0,amp_BAJD);
[h,p,ci,stats]=ttest2(amp_CAJD0,amp_CAJD);

