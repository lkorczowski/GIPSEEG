function result = PercentShow(value, maxValue, prefix, suffix, goPreviousLine)

% return how many percents is value of maxValue
% 
% value          -- current value (x%, 0<=x<=100)  
% maxValue       -- max possible value (100%)
% prefix         -- text before percent (default == [])
% suffix         -- text after percent (default == [])
% goPreviousLine -- if equals FALSE writes percent value on the new line,
%                   otherwise writes it on the previous line (default == true)
%
% result == 0 if everyting is OK, and == 1, otherwise
%
% EXAMPLE
%
% PercentShow(0, 100, [],[], false);
% for ii = 1:100
%      % do something here
%      PercentShow(ii, 100);
% end

    if nargin < 3
        prefix = [];
    end
    if nargin < 4
        suffix = [];
    end
    if nargin < 5
        goPreviousLine = true;
    end

    if any([value<0 maxValue<0 value>maxValue])
        result = 1;
        return
    end
    
    str = [prefix sprintf('%3d%%',round(100*value/maxValue)) suffix];
    
    if goPreviousLine
        str = [sprintf(char(repmat(double('\b'),1,length(str)+1))) str];
    end
    disp(str);
    
    result = 0;
    
end