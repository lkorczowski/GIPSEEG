function [] = dispBootstrap(VAL_OBS, VAL_BOOT, DISP_SIG_LEVELS)
%function [] = dispBootstrap(VAL_OBS, VAL_BOOT, DISP_SIG_LEVELS)

% ************************************************************************
% This function displays results of bootstrap analysis, in the form of a
% bar plot of bootstrap values (usually with pseudo-normal distribution), 
% along with observation value and its corresponding p-value.
%
% Input:
% ------
% VAL_OBS 	--> scalar, observation value
% VAL_BOOT	--> vector of values computed with bootstrap sampling
% DISP_SIG_LEVELS 	--> (optional) fill curves to display standard deviations  
%
%
% History:
% --------
% First version: 2013-06-20
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% *** assess inputs
if (nargin < 3) || (isempty(DISP_SIG_LEVELS))
    DISP_SIG_LEVELS = 1;
end
if (size(VAL_BOOT,2)==1)
    VAL_BOOT = VAL_BOOT';
end

% *** compute histogram values 
nbBoot  = size(VAL_BOOT,2);
xMin    = .8*abs(min([VAL_OBS VAL_BOOT]))*sign(min([VAL_OBS VAL_BOOT]));
xMax    = 1.2*abs(max([VAL_OBS VAL_BOOT]))*sign(max([VAL_OBS VAL_BOOT]));
bins    = [xMin:.005:xMax];
[histo_val, xout] = hist(VAL_BOOT,bins,'FaceColor','r','EdgeColor','w');
area    = sum(histo_val) * (xout(2)-xout(1));
p_val   = (1+sum(histo_val(xout>VAL_OBS)))/nbBoot;  % this is the p-value

% *** plot them
h = bar(bins,histo_val,'style','histc');
set(h,'FaceColor','g','EdgeColor','w'),
hold on,
line([VAL_OBS VAL_OBS],get(gca,'Ylim'),'Color','r','linewidth',2,'linestyle','-');
xlim([xMin xMax]),
legend('Surrogate','Experimental');
[fi,xi] = ksdensity(VAL_BOOT); % compute smooths density estimate
text(1.1*VAL_OBS,.5*max(get(gca,'ylim')),['p = ' num2str(p_val)],'color','r'),      % num2str(p_val,'%10.2e')

% *** fill area under curve for sig levels 0.05, 0.01 and 0.001
if (DISP_SIG_LEVELS ==1)
    boot_levels_ix  = nbBoot*[.05 .01 .001]; % Along Gaussian curve: +1sd = 84.13th prctile, +2sd = 97.72nd, and +3sd = 99.87th 
    sorted_val_boot = sort(VAL_BOOT,'descend'); 
    bootLevelVals   = sorted_val_boot(boot_levels_ix);
    fill([xi(xi>bootLevelVals(1)) xi(end) xi(find(xi>bootLevelVals(1),1))], ...
            [area*fi(xi>bootLevelVals(1)) 0 0] , [255 255 100]./255, 'FaceAlpha', 0.7); 
    fill([xi(xi>bootLevelVals(2)) xi(end) xi(find(xi>bootLevelVals(2),1))], ...
            [area*fi(xi>bootLevelVals(2)) 0 0] , [255 100 50]./255, 'FaceAlpha', 0.8); 
    fill([xi(xi>bootLevelVals(3)) xi(end) xi(find(xi>bootLevelVals(3),1))], ...
            [area*fi(xi>bootLevelVals(3)) 0 0] , [255 50 50]./255, 'FaceAlpha', 0.9); 
    legend('Surrogate','Experimental','p < 0.05','p < 0.01','p < 0.001'),
end

% plot smooth density estimate (must be done here for legend purpose)
plot(xi,area*fi,'color','k','linewidth',1);
        
