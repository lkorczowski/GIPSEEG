%% exemple on the use of meanOverlap
clear all
clc
close all
N=3000;%size of the measured data
n=128;%size of the epoch
t=1:n;
s1=sin(t*100);%simulate typical response
s1=[s1;sin(t*20)];
RNDstep=100
Tags=[];
Weights=[];
indF=0;
NoiseVar=0.00;
while indF<N
    
    if isempty(Tags),Tags=randi(RNDstep);else, Tags=[Tags Tags(end)+randi(RNDstep)];end%flashs indices  
    indF=Tags(end)+n;
    if indF>N
        Tags(end)=[];
    else
        Weights=[Weights rand(1)];
    end
end
%Weights=Weights-mean(Weights)+1;
 Weights(10:20)=0.1;
 AvWeight=mean(Weights);
 Weights=Weights/mean(Weights(Weights~=0));  %normalizing mean

 s1n=s1+randn(size(s1))*NoiseVar;
figure
subplot1(5,1,'Gap',[0 0.05])
subplot1(1)

plot(s1');title('stereotypical response G_z')
Flash=zeros(1,N);%tag vector for the flash
Flash(Tags)=Weights;%add some overlapping flash
Top=toeplitz(zeros(1,n),Flash);
% Weights(10:1000)=0.01; %wrong weight estimation
% Weights=Weights/mean(Weights(Weights~=0));  %normalizing mean
E=s1n*Top;
subplot1(2)
plot(E');title(['Signal with overlapping flash (No Noise), K_z=' num2str(length(Weights))])
%use standard method to estimate
indF=find(Flash);
for i=1:length(indF)
    epoch(:,:,i)=E(:,indF(i):indF(i)+n-1)/Weights(i);
end
%average
subplot1(3)
s1mean=mean(epoch,3);%(1/sum(Weights))*
e=s1-s1mean;

plot(s1mean');hold on;plot(e','r--');legend('$\hat{s}_1$','$\hat{s}_2$','$\hat{e}_1$','$\hat{e}_2$');text(0.5, 0.5,['error=' num2str(trace(e*e'))])
title('Weighted AEA (1.8)')
set(legend,'Interpreter','latex')
Flash=zeros(1,N);%tag vector for the flash
Flash(Tags)=Weights;%add some overlapping flash
Top=toeplitz(zeros(1,n),Flash);
s1top=E*Top'*pinv(Top*Top');
s1top=(Top'\E')';

e=s1-s1top;
subplot1(4)
plot(s1top');hold on;plot(e','r--');legend('$\hat{s}_1$','$\hat{s}_2$','$\hat{e}_1$','$\hat{e}_2$');text(0.5, 0.5,['error=' num2str(trace(e*e'))])
set(legend,'Interpreter','latex')
title('Weighted Regression T\X')
tic
Weights(10)=1000000
[Emean Class]=meanOverlap(E,Flash,n,[],Weights);
toc
subplot1(5)
e=s1-Emean;
plot(Emean');hold on;plot(e','r--');legend('$\hat{s}_1$','$\hat{s}_2$','$\hat{e}_1$','$\hat{e}_2$');text(0.5, 0.5,['error=' num2str(trace(e*e'))])
set(legend,'Interpreter','latex')
title('Weighted Regression T\X + Weights correction (1.9)')

    set(gcf, 'color', [1 1 1])
set(gcf, 'PaperPosition', [0 0 20 16],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
print(gcf, ['.\tuto_meanOverlap'],'-dtiff','-r450')

