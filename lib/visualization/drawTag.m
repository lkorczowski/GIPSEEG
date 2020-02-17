function drawTag(Flash,Fs,Color,LineStyle)
if nargin<4
   LineStyle='-'; 
end
if nargin<3
    Color='r';
end
if nargin<2
   Fs=1; 
end
yL = get(gca,'YLim');
indF=find(Flash);
indF=indF/Fs;
for nf=1:length(indF)
line([indF(nf) indF(nf)],yL,'Color',Color,'LineStyle',LineStyle,'LineWidth',0.05);
end
end