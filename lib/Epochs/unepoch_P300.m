function Xout=unepoch_p300(X)
% Xout=unepoch_p300(X)
% this function generate reshape a 3D matrix into 2D matrix
% see also: epoch_p300, epoch_EEG_struct
l1=size(X,1);
l2=size(X,2);
l3=size(X,3);
Xout=zeros(l1*l3,l2);
for i=1:size(X,3)
    Xout((i-1)*l1+1:i*l1,:)=X(:,:,i);
end