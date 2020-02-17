function barPlot(VAL,pVAL,THRESH,POS)
% this function displays a graph bar with asterisks for statistical
% significance levels 
%
% TODOOOOOOOOOO: function that merges 'barPlot' and 'errorb'
%
% input:
% VAL       --> bar values. vector [1,n], or matrix [n,m] for n plots of
%               groups of m variables
% pVAL      --> (optionnal) associated significance values [1*n]
% THRESH    --> (optionnal) threshold for plotting an aterisk,
%                           can be a single value or a vector of up to three 
%                           increasing values (then up to 3 asterisks are drawn)
% POS       --> (optionnal) vertical position where to plot asterisks
%
%  *** History *** 
% First version: 20/07/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% default threshold
if(nargin < 3)||(isempty(THRESH))   THRESH = 0.01 ;  end
if (size(VAL,1)==1)
    isOneGroup = 1; % only one group to plot
else
    isOneGroup = 0;
end
    
% plot bars
if isOneGroup
    bar(VAL,'group','BarWidth',.6); 
    nbVal = length(VAL);
else
    bar(VAL,'group','BarWidth',1); 
    nbVal = size(VAL,1);
end

% plot asterisks
if(nargin > 1)
    if isOneGroup
        previous_maskThresh = false(1,1);
    else
        previous_maskThresh = false(1,length(VAL));
    end
    for thresh_ix = length(THRESH):-1:1
        if      (thresh_ix ==1)     
            asterisk_str = ' * '; 
        elseif  (thresh_ix ==2)     
            asterisk_str = '** ';
        elseif  (thresh_ix ==3)     
            asterisk_str = '***';
        else
            error('[barPlot] max 3 thresholds possible');  
        end
        maskThresh =  (pVAL<THRESH(thresh_ix)) & (~previous_maskThresh) ;
        previous_maskThresh = previous_maskThresh | maskThresh;
        if(nargin < 4)
            if isOneGroup
                vertPos = 1.5*max(VAL);
            else
                vertPos = 1.25*VAL(find(maskThresh));
            end
        else
            vertPos = repmat(POS,1,length(find(maskThresh)));
        end
        if find(maskThresh)
            if isOneGroup
                text(find(maskThresh)+.4,vertPos,asterisk_str, 'FontSize', 18,'FontWeight', 'Bold');
            else 
                text(find(maskThresh)-.05,vertPos,asterisk_str, 'FontSize', 18,'FontWeight', 'Bold');
            end
        end
    end
%     ylim([0 1.2*max(max(VAL))]),
end

   

end