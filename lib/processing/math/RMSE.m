function [yrmse, RMSEelec]=RMSE(y_pred,y)
% Return for the 4D matrix y the root mean squared error according to the
% estimation y_pred (3D matrix)
% y_pred [nb_samples x nb_electrodes x nb_classes]
% y [nb_samples x nb_electrodes x nb_classes x nb_observations]
yrmse=[];
if nargin>1
for indC=1:size(y,3)
        yrmse(:,:,indC) = sqrt(mean((squeeze(y(:,:,indC,:))-repmat(y_pred(:,:,indC),[1,1,size(y,4)])),3).^2);
end
else
    y=y_pred;
    ypred=mean(y,4);
    for indC=1:size(y,3)
        yrmse(:,:,indC) = sqrt(mean((squeeze(y(:,:,indC,:))-repmat(y_pred(:,:,indC),[1,1,size(y,4)])),3).^2);
end

end