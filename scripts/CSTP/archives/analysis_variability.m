%analysis variability
%Zbarz (3)

%40
close all
Mapping={'Fp1';%1
'Fp2';%2
'F3';%3
'AFz';%4
'F4';%5
'T7';%6
'Cz';%7
'T8';%8
'P7';%9
'P3';%10
'Pz';%11
'P4';%12
'P8';%13
'O1';%14
'Oz';%15
'O2'}%16
ratios=[ 40 20 10 5]

electrode=1
class=2
figure
Scale=12
for j=1:size(Save,3)
h = text(-0.25, 0.5, 'row 2');
set(h, 'rotation', 90)

hold on
for i=1:size(Save,2)
P(:,:,i)=Save{10,i,j}(:,:,class);
Zbarz_noweight(:,:,i)=Save{11,i,j}(:,:,class);
Zbarz_weight(:,:,i)=Save{3,i,j}{1}(:,:,class);
Zbarz(:,:,i)=Save{3,i,j}{end}(:,:,class);
end
hold off
meanZbarz=mean(Zbarz,3);
varZbarz=var(Zbarz,[],3);
hold on;
subplot(4,size(Save,3),1+size(Save,3)*(j-1))
plotEEGvariability(P,electrode);
axis([0 128 -Scale Scale])
if j==1; title('X_{bar}'); end

ylabel([num2str(ratios(j)) ' sweep'])

subplot(4,size(Save,3),2+size(Save,3)*(j-1))
plotEEGvariability(Zbarz_noweight,electrode);
axis([0 128 -Scale Scale])
if j==1; title('X^{barhat}_{noweight}'); end


subplot(4,size(Save,3),3+size(Save,3)*(j-1))
plotEEGvariability(Zbarz_weight,electrode);
axis([0 128 -Scale Scale])

if j==1; title('X^{barhat}_{\sigma}'); end


subplot(4,size(Save,3),4+size(Save,3)*(j-1))
[h1 h2 h3]=plotEEGvariability(Zbarz,electrode);
axis([0 128 -Scale Scale])
if j==1; title('X^{barhat}_{\sigma,\epsilon}'); end

if j==1; legend([h1 h2 h3],{'all','mean','2std'},'location','NorthEastOutside'); end


end
text(-250.5,-020,['Subject' num2str(Subjects) ' class ' num2str(class) ' ' Mapping{electrode}] ,'FontSize',20)
