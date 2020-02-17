function Scale=FindScaling(PLOTX,winElec,winTime)
if nargin<2
winElec=1:size(PLOTX,1);
end
if nargin<3
winTime=1:size(PLOTX,2);
end
Scale=max(max(abs(mean(PLOTX(winElec,winTime,:),3))))/2;