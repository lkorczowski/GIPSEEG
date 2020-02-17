function result = TextShow(text, goPreviousLine)

% show any text on the current or previous line
% 
% text           -- string which will be shown on the screen
% goPreviousLine -- if equals FALSE writes string on the new line,
%                   otherwise writes it on the previous line (default == true)
%
% result == 0 if everyting is OK, and == 1, otherwise
%
% EXAMPLE
%
% TextShow('test: ', false);
% for ii = 1:5
%      % do something here
%      TextShow([num2str(ii),', '], true);
% end

    if nargin < 2
        goPreviousLine = true;
    end
    try
        if goPreviousLine
            text = [sprintf(char(repmat(double('\b'),1,1))) text];
        end
        disp(text);
    catch ME
        errorMessageToFile('Err.log',ME.message);
        result = 1;
        return
    end
    
    result = 0;
    
end