function h=autoSubplot(NbColums,varargin)

nbplot=length(varargin);
NbRows=nbplot/NbColums;
if NbRows-floor(NbRows)
    disp('Number of argument not consistent with the number of Colums')
    return
end
for i=1:nbplot
h(i)=subplot(NbRows,NbColums,i);h=plot(varargin{i});
try
catch e
x = get(h,'XData'); % Get the plotted data
y = get(h,'YData');
imin = find(min(y) == y); % Find the index of the min and max
imax = find(max(y) == y);
text(x(imin),y(imin),[' Minimum = ',num2str(y(imin))],...
	'VerticalAlignment','middle',...
	'HorizontalAlignment','left',...
	'FontSize',max(14-2*nbplot,8))
text(x(imax),y(imax),['Maximum =  ',num2str(y(imax))],...
	'VerticalAlignment','bottom',...
	'HorizontalAlignment','right',...
	'FontSize',max(14-2*nbplot,8))
title(['MaxMinRatio=' num2str(y(imax)/y(imin))])
end
end

