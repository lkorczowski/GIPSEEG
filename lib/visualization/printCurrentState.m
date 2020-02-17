function printCurrentState(INFO,previousLine)
% *** History: 19-Mar-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
if nargin<2
    previousLine=true;
end

str = [sprintf(INFO)];
    if previousLine
        str = [sprintf(char(repmat(double('\b'),1,length(str)+1))) str];
    end
    disp(str);