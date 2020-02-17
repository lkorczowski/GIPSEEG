function [yrmse]=RMSE2(y,y_pred)
% [yrmse, RMSEelec]=RMSE2(y_pred,y)
% Return for the root mean squared error between the prediction
%       y [nb_samples x nb_electrodes x nb_observations]
% and the prediction
%       y_pred [nb_samples x nb_electrodes] 
% by default y_pred=mean(y,3);
yrmse=[];
if nargin<2
    y_pred=mean(y,3);
end

    for indC=1:size(y,3)
        yrmse(indC) = sqrt(mean(mean((y(:,:,indC)-y_pred).^2)));
    end


end