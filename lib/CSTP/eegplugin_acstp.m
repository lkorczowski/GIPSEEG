% eegplugin_acstp() - EEGLAB plugin for importing CTF data files.
%
% Usage:
%   >> eegplugin_acstp(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Louis Korczowski, Arnaud Delorme, Marco Congedo, GIPSA-Lab, 2015
% Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% See also: pop_acstp()

% Copyright (C) 2015 Arnaud Delorme
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_acstp(fig, trystrs, catchstrs)

    vers = 'acstp';
    if nargin < 3
        error('eegplugin_acstp requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if exist('pop_acstp', 'file')
        p = which('eegplugin_acstp.m');
        p = p(1:findstr(p,'eegplugin_acstp.m')-1);
        addpath(p);
        addpath(fullfile(p, 'function'));
    end;
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'tools');
    
    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_acstp(EEG);' catchstrs.new_and_hist ];
    
    % create menus
    % ------------
    uimenu( menu, 'label', 'Remove artifacts using ACSTP', 'callback', comcnt, 'separator', 'on');
