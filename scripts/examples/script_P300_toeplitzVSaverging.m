%% exemple on the use of meanOverlap
clear all
clc
close all
N=60000;%size of the measured data
n=1024;%size of the epoch
t=1:n;
s1=sin(t*100);%simulate typical response
s1=[s1;sin(t*20)];
RNDstep=60
Tags=[];
Weights=[];
indF=0;
LEGEND={'s1','s2','e1','e2'}
nbRows=6;
nbColumns=1;
indP=1;
while indF<N
    
    if isempty(Tags),Tags=randi(RNDstep);else, Tags=[Tags Tags(end)+randi(RNDstep)];end%flashs indices
    indF=Tags(end)+n;
    if indF>N
        Tags(end)=[];
    else
        Weights=[Weights 1];
    end
end

figure
%%%%%%%%%%%%%%%%%%%%%%%
subplot(nbRows,nbColumns,indP);indP=indP+1;
plot(s1');title('stereotypical response');legend(LEGEND)
Flash=zeros(1,N);%tag vector for the flash
Flash(Tags)=Weights;%add some overlapping flash
Top=toeplitz(zeros(1,n),Flash);
E=s1*Top;
E(:,500:505)=100;
Weights(find(Flash)>=(500-n) & find(Flash)<=505+n)=0.01%1
ylim([-2 2])

%%%%%%%%%%%%%%%%%%%%%%%
subplot(nbRows,nbColumns,indP);indP=indP+1;
plot(E');title(['Signal with overlapping, artefact(+100). K=' num2str(length(Weights))]);legend(LEGEND)
%use standard method to estimate
axis([400 600 -5 5])

%%%%%%%%%%%%%%%%%%%%%%% EAE without weights
indF=find(Flash);
for i=1:length(indF)
    epoch(:,:,i)=E(:,indF(i):indF(i)+n-1)*1;
end
subplot(nbRows,nbColumns,indP);indP=indP+1;
s1mean=mean(epoch,3);
e=s1-s1mean;
plot(s1mean');hold on;plot(e','r--');legend(LEGEND);text(0.5, 0.5,['error=' num2str(trace(e*e'))])
title('EAE without weights')
ylim([-2 2])
%%%%%%%%%%%%%%%%%%%%%%% EAE with estimated weights
indF=find(Flash);
for i=1:length(indF)
    epoch(:,:,i)=E(:,indF(i):indF(i)+n-1)*Weights(i);
end
subplot(nbRows,nbColumns,indP);indP=indP+1;
s1mean=mean(epoch,3);
e=s1-s1mean;
plot(s1mean');hold on;plot(e','r--');legend(LEGEND);text(0.5, 0.5,['error=' num2str(trace(e*e'))])
title('EAE with estimated weights')
ylim([-2 2])


%%%%%%%%%%%%%%%%%%%%%%% meanOver no weights
s1top=E*Top'*pinv(Top*Top');
[s1top]=meanOverlap(E,Flash,n);

e=s1-s1top;
subplot(nbRows,nbColumns,indP);indP=indP+1;
plot(s1top');hold on;plot(e','r--');legend(LEGEND);text(0.5, 0.5,['error=' num2str(trace(e*e'))])
tic
title('meanOverlap without weights')
ylim([-2 2])

%%%%%%%%%%%%%%%%%%%%%%% meanOver weights
tic
[Emean Class]=meanOverlap(E,Flash,n,[],Weights);
toc
subplot(nbRows,nbColumns,indP);indP=indP+1;
e=s1-Emean;
plot(Emean');hold on;plot(e','r--');legend(LEGEND);text(0.5, 0.5,['error=' num2str(trace(e*e'))])
title('meanOverlap with estimated weights')
ylim([-2 2])

print(['.\figures\CSTP\empirical\simulated'],'-dtiff','-r450')

%%
    