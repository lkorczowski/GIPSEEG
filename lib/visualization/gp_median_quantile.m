function [Y errY]=gp_median_quantile(X,p)
% [Y errY]=gp_median_quantile(X,p)
% return the median Y and the (p)-quantile lower and upper bound errY
% of data X

if nargin<2
    p=0.05;
end

if iscell(X)
    for indR=1:length(X)
        X{indR}=X{indR}(~any(isnan(X{indR}),2),:);
            for indC=1:size(X{indR},2)
            Y(indR,indC)=median(X{indR}(:,indC));
            errY(indR,indC,1)=Y(indR,indC)-quantile(X{indR}(:,indC),p);
            errY(indR,indC,2)=Y(indR,indC)-quantile(X{indR}(:,indC),1-p);
        end
        
    end
else
                    X=X(~any(isnan(X),2),:,:);

    for indR=1:size(X,3)

           for indC=1:size(X,2)
            Y(indR,indC)=median(X(:,indC,indR));
            errY(indR,indC,1)=Y(indR,indC)-quantile(X(:,indC,indR),p);
            errY(indR,indC,2)=Y(indR,indC)-quantile(X(:,indC,indR),1-p);
        end
        
    end
end


%nested
%     function [Y errY]=gp_median_quantile2(X,p)
%         for indC=1:size(X,2)
%             Y(1,indC)=median(X(:,indC));
%             errY(1,indC,1)=Y(1,indC)-quantile(X(:,indC),p);
%             errY(1,indC,2)=Y(1,indC)-quantile(X(:,indC),1-p);
%         end
%     end
end