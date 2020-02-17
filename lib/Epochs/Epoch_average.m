function P=Epoch_average(X,Y,NbUsers,P300_ref_orientation)


if nargin < 4
    P300_ref_orientation = 'commonP1';%or  'multiP1'
end
if nargin < 3
    NbUsers=1;
end

if NbUsers==1
    [COV P] = covariances_p300(X,Y,method_cov,arg_cov);
else
    labels = unique(Y);
    Nclass = length(labels);
        P = cell(Nclass,1);
    Nbe=size(X,1)/NbUsers;
    
    if size(Y(:,:,1))~=size(X(:,:,1))
        for i=1:Nclass;
        switch P300_ref_orientation
            case 'commonP1'
                P{i} = mean(X(:,:,Y==labels(i)),3);
                P{i}=transpose(mean(reshape(P{i}',size(P{i}',1),Nbe,NbUsers),3));
                %mean P1 on all users
            case 'multiP1'
                for Us=1:NbUsers
                    P{i}(:,:,Us) = mean(X((Nbe*(Us-1)+1):(Nbe*(Us)),:,Y==labels(i)),3);
                end
        end
        end
    else
        disp('Target file does not match epoch file')
    end
end