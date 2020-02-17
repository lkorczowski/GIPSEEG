function auc = area_under_curve(out,truelabels)

[tpr,fpr,~] = roc(truelabels,out);

auc =  ((1-sum((tpr(2:end)-tpr(1:end-1)).*(fpr(2:end)+fpr(1:end-1))))+1)/2;
