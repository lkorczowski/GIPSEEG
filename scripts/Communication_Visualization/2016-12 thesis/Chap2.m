Dir='D:\Mes Documents GIPSA\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\images\'
Dir='C:\Users\Oui-Oui\Mes Documents GIPSA\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\images\'
% Dir='.\test\'
%% Exemple CSP
%%%%parameters%%%%%
Paper=[0 0 20 10];factor=0.25;FontS=12;Font='times new roman';
%%%%%%%%%
close all
Cov1=[100 70 ; 70 70];Mu1=[0 0];
Cov2=[30 30 ; 30 100];Mu2=[0 0];
R1 = mvnrnd(Mu1,Cov1,5000);
R2 =[mvnrnd(Mu2,Cov2,5000)];

figure
subplot1(1,2)
subplot1(1)

plot(R1(:,1),R1(:,2),'.b');
hold on;


plot(R2(:,1),R2(:,2),'.r')
hold off;
axis([-150 150 -150 150]*factor)
% xlabel('$\mathbf{x}_1$','interpreter','latex','fontsize',FontS,'fontname',Font)
% ylabel('$\mathbf{x}_2$','interpreter','latex','fontsize',FontS,'rotation',-0,'fontname',Font)

S1=cov(R1);S2=cov(R2);

% Compute CSP
[V,D]=eig(S1,S1+S2);
%[V,D]=eig(inv(S1)*S2);

% Apply CSP
F1=(R1*V);
F2=(R2*V);

subplot1(2)
plot(F1(:,1),F1(:,2),'.b');
hold on;


plot(F2(:,1),F2(:,2),'.r')
hold off;

%%%%%%%% SAVE %%%%%%%%%%%%
set(gca,'yaxislocation','right');
set(gca,'fontsize',FontS)
% xlabel('$\mathbf{w}_1^T \mathbf{X}$','interpreter','latex','fontsize',FontS,'fontname',Font)
% lh=ylabel('$\mathbf{w}_P^T \mathbf{X}$','interpreter','latex','fontsize',FontS,'rotation',-0,'fontname',Font)
%      p=get(lh,'position');
%      set(lh,'position',p+.1);

set(gcf,'color',[1 1 1],'paperposition',Paper)
% axis([-3.5 3.5 -7 7])
print(gcf, [ Dir '2_CSP'],'-dpng','-r450')
%%%%%%%% SAVE %%%%%%%%%%%%



%% Principe LDA
Paper=[0 0 20 10];factor=0.25;FontS=12;Font='times new roman';
% Dir='./figures/thesis/';
Dir=PATH;
%%%%%%%%%
close all
Cov1=[100 70 ; 70 70];Mu1=[5 10];
Cov2=[30 30 ; 30 100];Mu2=[-5 +50];
R1 = mvnrnd(Mu1,Cov1,5000);
R2 =[mvnrnd(Mu2,Cov2,5000)];

figure
subplot(3,3,[1 2 4 5])
plot(R1(:,1),R1(:,2),'.b');
hold on;

plot(R2(:,1),R2(:,2),'.r')
hold off;
set(gca,'YTick',[]);
set(gca,'XTick',[]);
% Y-axis projection
subplot(3,3,[3 6])
plot(1,R1(:,2),'.b');
hold on;
plot(1.1,R2(:,2),'.r');
hold off;
xlim([0 2]);set(gca,'XTick',[]);%set(gca,'YTick',[]);
 set(gca,'yaxislocation','right');

% X-axis projection
subplot(3,3,[7 8])
plot(R1(:,1),1,'.b');
hold on;
plot(R2(:,1),1.1,'.r');
hold off;
ylim([0 2]);%set(gca,'XTick',[]);
set(gca,'YTick',[]);


%calcul LDA

mu1=mean(R1);mu2=mean(R2);
cov1=cov(R1);cov2=cov(R1);

Sigv=cov1+cov2;
w=inv(Sigv)*[mu1-mu2]'

x= linspace(-25,40, 100);y= w(1)-w(2)*x+25;
% w*R1;
subplot(3,3,[1 2 4 5]);hold on;
plot(x,y);hold off

%%%%%%%% SAVE %%%%%%%%%%%%
set(gca,'fontsize',FontS)
% xlabel('$\mathbf{w}_1^T \mathbf{X}$','interpreter','latex','fontsize',FontS,'fontname',Font)
% lh=ylabel('$\mathbf{w}_P^T \mathbf{X}$','interpreter','latex','fontsize',FontS,'rotation',-0,'fontname',Font)
%      p=get(lh,'position');
%      set(lh,'position',p+.1);

set(gcf,'color',[1 1 1],'paperposition',Paper)
% axis([-3.5 3.5 -7 7])
print(gcf, [ Dir '2_LDA'],'-dpng','-r450')
%%%%%%%% SAVE %%%%%%%%%%%%
