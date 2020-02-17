rng('shuffle') %set schuffle to get a new random serie
close all
%generate 5 vectors
%2 small vectors

x= rand(2,1)*1-0.5;
y= rand(2,1)*1.2-0.5;
% 2 averages
x= [x; rand(2,1)*1.5-0.75];
y= [y;rand(2,1)*1.5];
%1 big
x= [x; rand(1,1)*.1-1];
y= [y;rand(1,1)*0-2];

z=x+i*y
%
figure
subplot(121)
L=0;% counter of lines
% random vector
compass(z,'-b')
h=findall(gca);
set(h(L+1:end),'LineWidth',1.5)
L=length(h);

MeanVec=mean(x)+i*mean(y);
hold on
compass(MeanVec,'-k')
h=findall(gca);
set(h(1:3),'LineWidth',2)
L=length(h);
xlabel('re^{i\phi}','fontsize',18,'fontname','times new roman','Fontangle','italic')

subplot(122)
% random vector
compass(z./abs(z),'-b')
h=findall(gca);
set(h,'LineWidth',1.5)

MeanVec=mean(z./abs(z))
hold on
compass(MeanVec,'-k')
h=findall(gca);
set(h(1:3),'LineWidth',2)
xlabel('e^{i\phi}','fontsize',18,'fontname','times new roman','Fontangle','italic')

print(gcf, ['.\compass_test1'],'-dtiff','-r450')





%% univariate measurement
figure
subplot(121)
L=0;% counter of lines
% random vector
compass(z,'-b')
h=findall(gca);
set(h(L+1:end),'LineWidth',1.5)
L=length(h);

phi=angle(z);
Pcon=exp(i*mean(phi))

hold on
compass(Pcon,'-k')
h=findall(gca);
set(h(1:3),'LineWidth',2)
L=length(h);
xlabel('Pcon','fontsize',18,'fontname','times new roman','Fontangle','italic')

print(gcf, ['.\compass_test2'],'-dtiff','-r450')


