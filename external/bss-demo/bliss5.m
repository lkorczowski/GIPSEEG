function varargout = bliss
%    bliss5 Application M-file for bliss.fig, designed for matlab-5
%    fig = bliss5 launch bliss GUI.

% Adapted from the version for matlab-6 (Jitesh Shah) by Pham,
% August 11, 2002
% _requires bliss.m to function_

path(path,[pwd '/matlab-6']);
fig = openfig('bliss');
% Use system color scheme for figure:
    
set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
      
% Generate a structure of handles to pass to callbacks, and store it.
handles = guihandles(fig);
guidata(fig, handles);
   
if nargout > 0
  varargout{1} = fig;
end
 
global s1 s2 X;
s1 = [];
s2 = [];
X = [];
			    
%%%%%%%%%%%%%%% private version openfig of matlab-5 %%%%%%%%%%%%%%%%%

function h = openfig(filename)
%   OPENFIG   Loads HG object from FIG-file.
%   Helper function for OPEN.
%
%   See OPEN.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 1998/10/20 20:19:18 $
%   D. Foti  11/10/97

if isempty(find(filename == '.'))
  filename = [filename '.fig'];
end

fileVars = load(filename, '-mat');
varNames = fieldnames(fileVars);
index = strmatch('hgS', varNames);
if length(index) ~= 1
  error('invalid Figure file format');
end

varName = varNames{index};
versionStr = varName(find(varName == '_')+1:end);
versionNum = str2double(versionStr);

if versionNum > 50200
  fileVersion = sprintf('%d.%d.%d', str2double(versionStr(1:2)), ...
      str2double(versionStr(3:4)), str2double(versionStr(5:6)));
  warning(['Figure file created with a newer version (' fileVersion ') of MATLAB']);
end

hgS = getfield(fileVars, varNames{index});

numHandles = prod(size(hgS));
parent = zeros(numHandles, 1);
for i = 1:numHandles
  switch(hgS(i).type)
    case 'root'
      parent = 0;
    case 'figure'
      parent(i) = 0;
    case {'axes' 'uimenu' 'uicontextmenu' 'uicontrol' 'uitoolbar'}
      parent(i) = gcf;
    case {'uipushtool' 'uitoggletool'}
      parent(i) = gctb;
    otherwise
      parent(i) = gca;
  end
end

h = struct2handle(hgS, parent, 'convert');

% If we loaded only one top-level object, and it was a figure,
% remember the name of the file from which it was loaded (including
% full path and extension) so the Save option from the File menu
% can overwrite it without having to prompt for a file name.
if (numHandles == 1) & strcmp(get(h,'type'),'figure')
  set(h,'FileName',which(filename));
end

function tb = gctb
% find the first uitoolbar in the current figure, creating one if
% necessary

tb = findobj(gcf,'type','uitoolbar');
if ~isempty(tb)
  tb = tb(1);
else
  tb = uitoolbar;
end
