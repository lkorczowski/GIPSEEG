function varargout = bss(varargin)
%    bss Application M-file for bss.fig
%    fig = bliss launch bss GUI.
%    bliss('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 19-Aug-2002 15:03:04
% Adapted from the verion of Jitesh Shah) by Pham,  August 10, 2002

if nargin == 0  % LAUNCH GUI
  fig = openfig(mfilename, 'reuse');
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

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
  catch
    disp(lasterr);
  end
end

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.close.
close(handles.bliss);


% --------------------------------------------------------------------
function varargout = help_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.help.
%errordlg('help not implemented yet');
str={'This software implements three methods for blind source separation.'
'To start with the software, choose both sources from the list by'
'clicking on "Select Source". You are required to fill in the parameters'
'required for each of the sources.'
''
'In case you choose a wav file, a small window will pop-up with the'
'entire wave file plotted in it. You are to move your mouse and select'
'the starting point of the portion (1000 samples) of the source.'
''
'After selecting Source 1 and Source 2, you choose the compsition of'
'the mixture (manually or randomly) then click on "Accept".'
''
'Finally, you choose the separation algorithm to separate the sources.'
'The global matrix (the product of the computed separation matrix'
'and the mixing matrix) is then plotted, among others.'
''
'You can also apply a postnonlinear transform of the data before'
'applying the algorithm. In this case you should choose a fourth'
'method: "Mutual Information for postnonlinear mixture"'
''
'Software created by Dinh-Tuan Pham (Dinh-Tuan.Pham@imag.fr)'
'and Jitesh Shah (jiteshis@hotmail.com).'};
helpdlg(str,'Help for BSS-DEMO');

% --------------------------------------------------------------------
function varargout = popupmenu2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu2.
global s1 s2
if (get(handles.sourceno,'value')~=1)
  str = {'Sine wave'
         'Lin. combination of 2 sine waves'
	 'Aleluya (WAV)'
	 'Goodmorning (WAV)'
	 'CLIP (WAV)'
	 'Libertad (WAV)'
	 'ARMA process'};
  [selec,v] = listdlg('PromptString','Select the source:', ...
		      'SelectionMode','single','ListString',str,...
		      'Listsize',[200 150]);
  if(v~=0)
    switch selec
      case 1
        prompt={'Frequency (per 1000 samples) :'};
        def={'8'};
        dlgTitle='Input for sine wave';
  	answ=inputdlg(prompt,dlgTitle,1,def);
  	if(length(answ)~=0)
  	  f=str2double(answ)/1000;
  	  s=sin(2*pi*(f*(1:1000))+rand)*sqrt(2);% scale to have unit power 
  	end
      case 2
        prompt={'Frequency of 1st sine wave (per 1000 samples) :'
        	'Frequency of 2nd sine wave (per 1000 samples) :'
  	        'Amplitude of 2nd sine wave (relative to the 1st) :'};
        def={'9','10','1'};
        dlgTitle='Input for 2 sine waves';
        answ=inputdlg(prompt,dlgTitle,1,def);
        if(length(answ)~=0)
          f1=str2double(answ(1))/1000;
          f2=str2double(answ(2))/1000;
          a=str2double(answ(3));
          s=sin(2*pi*(f1*(1:1000)+rand))+a*sin(2*pi*(f2*(1:1000)+rand)) ...
	  	*sqrt(2/(1+a^2));          % scale to have unit power
         end
      case 3
        [in1,fs,bits]=wavread('aleluya.wav');
        figure('name','Click on the graph');
        plot(in1);
        title('Click on the graph from whence you want to load the signal');
        [a,b]=ginput(1);
        close('Click on the graph');
        a=min(fix(a),length(in1)-1000);
        s=in1(a+1:(a+1000))';
	s=s/sqrt(sum(s.^2/1000));	% scale to have unit power
      case 4
        [in1,fs,bits]=wavread('goodmorning.wav');
        figure('name','Click on the graph');
        plot(in1);
        title('Click on the graph from whence you want to load the signal');
        [a,b]=ginput(1);
        close('Click on the graph');
        a=min(fix(a),length(in1)-1000);
        s=in1(a+1:(a+1000))';
	s=s/sqrt(sum(s.^2/1000));	% scale to have unit power
      case 5
        [in1,fs,bits]=wavread('clip1.wav');
        figure('name','Click on the graph');
        plot(in1);
        title('Click on the graph from whence you want to load the signal');
        [a,b]=ginput(1);
        close('Click on the graph');
        a=min(fix(a),length(in1)-1000);
        s=in1(a+1:(a+1000))';
	s=s/sqrt(sum(s.^2/1000));      % scale to have unit power
      case 6
        [in1,fs,bits]=wavread('Libertad.wav');
        figure('name','Click on the graph');
        plot(in1);
        title('Click on the graph from whence you want to load the signal');
        [a,b]=ginput(1);
        close('Click on the graph');
        a=min(fix(a),length(in1)-1000);
        s=in1(a+1:(a+1000))';
	s=s/sqrt(sum(s.^2/1000));      % scale to have unit power
      case 7
%     x(i) = a(1)*x(i-1) + a(2)*x(i-2) + e(t) + b(1)e(t-1)
        prompt={'a1 (between -(1-a1) and 1-a1) :'
  	        'a2 (between -1 and 1) :'
  	        'b :'};
        dlgTitle='x(i)=a1x(i-1)+a2x(i-2)+e(t)+be(t-1)';
        defans={'1.4','-.9','0'};
        answ=inputdlg(prompt,dlgTitle,1,defans);
        if(length(answ)==0)
          errordlg('Incomplete data');
          return;
        end
        a=[str2double(answ(1)); str2double(answ(2))];
        b=str2double(answ(3));
        r1 = a(1)/(1-a(2));			% 1st corr of the AR model
        v = (1-r1^2)*(1-a(2)^2)/(1+b^2+2*r1*b);	% inn. variance to yield
        [s,e]=genearma(a,b,1000,v);		% unit power of the process
        if(e~=0)
          return;
        end
        s=s';
    end
    switch (get(handles.sourceno,'value'))
      case 2
        axes(handles.axes1);
        plot(s);
        s1 = s;
      case 3
        axes(handles.axes2);
        plot(s);
        s2 = s;
      otherwise
        errordlg('Invalid Selection');
    end
  end
end

% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit1.
%disp('edit1 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit2.
%disp('edit2 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit3.
%disp('edit3 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit4.
%disp('edit4 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = Random_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Random.
randnumbers = randn(2,1)./randn(2,1);
if prod(randnumbers) > 1; randnumbers = 1./randnumbers; end
set(handles.m12,'string',num2str(randnumbers(1)));
set(handles.m21,'string',num2str(randnumbers(2)));


% --------------------------------------------------------------------
function varargout = Accept_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Accept.
global s1 s2 X X0;
if (length(s1)==1000) & (length(s2)==1000)
  X = [1 str2num(get(handles.m12,'string'))
       str2num(get(handles.m21,'string')) 1]*[s1; s2];
  axes(handles.axes3);
  plot(X(1,:));
  axes(handles.axes4);
  plot(X(2,:));
  X0 = X;
  axes(handles.Axes2); cla
  axes(handles.Axes3); cla
else
  errordlg('Selection not completed');
end


% --------------------------------------------------------------------
function varargout = Postnonlin_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Postnonlin.
global X X0;
if size(X,2)~=1000
  errordlg('Selection not completed');
  return
end
obsno = get(handles.Postnonlin,'value')-1;
if obsno
  str = {'x -> tanh(ax) + bx'
	 'x -> x^3 + bx'};
  [selec,v] = listdlg('PromptString','Select la transformation:', ...
		      'SelectionMode','single','ListString',str,...
		      'Listsize',[160 50]);
  if v
    switch selec
      case 1
	answ=inputdlg({'Coefficient a :'
		       'Coefficient b :'},'x -> tanh(ax) + bx',1,{'4','0.1'});
	if (length(answ)~=0)
	  a = str2double(answ(1));
	  b = str2double(answ(2));
	  X(obsno,:) = tanh(a*X(obsno,:)) + b*X(obsno,:);
	end
      case 2
	answ=inputdlg('Coefficient a :','x -> x^3 + ax',1,{'0.1'});
	if (length(answ)~=0)
	  a = str2double(answ(1));
	  X(obsno,:) = X(obsno,:).^3 + a*X(obsno,:);
	end
      otherwise
	errordlg('Invalid Selection');
    end
  end
  switch obsno
    case 1
      axes(handles.Axes2);
      plot(sort(X0(1,:)),sort(X(1,:)));
      axis tight; grid on
      axes(handles.axes3);
      plot(X(1,:));
    case 2
      axes(handles.Axes3);
      plot(sort(X0(2,:)),sort(X(2,:)));
      axis tight; grid on
      axes(handles.axes4);
      plot(X(2,:));
  end
  set(handles.PostnonlinText,'String','Post nonlinear transform functions');
end

% --------------------------------------------------------------------
function varargout = algo_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton6.
global X X0;
algo = get(handles.algo,'value');
if algo == 1; return; end
if size(X,2)~=1000
  errordlg('Invalid selection');
  return;
else
  switch algo
    case 2
      prompt={'Maximum number of interations:'
	      'Number of frequency channels (< 500)'
	      'Stationary block length of the data:'};
      def={'15','10','50'};
      dlgTitle='Input for the parameters for the algorithm';
      answ=inputdlg(prompt,dlgTitle,1,def);
      if ~isempty(answ)
	maxiter=str2double(answ(1));
	nfreqc=str2double(answ(2));
	bloclenc=str2double(answ(3));
	if bloclenc<=(2*nfreqc)
	  bloclenc=2*nfreqc;
	  uiwait(msgbox('The block length must be at least double the number of frequency channels for a 2 sources system','Error','modal'));
	  end
	blklen = [0 2 4 5 10 20 40 50 100 200 500 1000];
	i=length(blklen);
	while (blklen(i-1) + blklen(i))/2 > bloclenc
	  i = i-1;
	end			% will finish with i >= 2 if bloclenc >=2
	if bloclenc ~= blklen(i)
	  if blklen(i) >=(2*nfreqc)
            bloclenc = blklen(i);
	  else
            bloclenc = blklen(i+1);
	  end
	  uiwait(msgbox({'The block length is adjusted to ',num2str(bloclenc)},'Warning','modal'));
	end
	espc = 2*1e-8;
	[sep, seps] = sepagaus(X',nfreqc,bloclenc,eye(2),espc,maxiter);

% Syntaxe  [sep, seps] = sepagaus(data, nfreq, bloclen, sep, eps, maxiter)
%
% Separation of source through the Gaussian mutual information criterion.
% * data is a n by K matrix, K being the number of sources
% * nfreq is the number of positive frequency channels to use
% * bloclen is the "stationarity length" of the data, it must be
%   greater than nfreq*K and preferably a divisor of n, otherwise the
%   last incomplete block of data is discard.
% * sep is the initial separation matrix.
% * seps, if present, is a (iter+1)*K by K matrix containing the sequence
%   of iterated separation matrices (iter = the number of iterations)
% * eps defines the stoping criterion, which is that the square norm
%   (with respect to a certain metric) of the relative "gradient" is
%   less it. This square norm also equals approximatively the decrease
%   of the criterion for this step
% * maxiter is maximum number of iteration of the algorithm
%
% nfreq defaults to 2, bloclen to fix(n/2), sep to the identity matrix, 
% eps to K*(K-1)*1e-8 and maxiter to 20

	Y=sep*X;
      end

    case 3
      prompt={'maximum number of iterations: '
	      'Kernel bandwidth/the optimal value for Gaussian density:'};
      def={'15','3'};
      dlgTitle='Input for the parameters for the algorithm';
      answ=inputdlg(prompt,dlgTitle,1,def);
      if ~isempty(answ)
        maxiter=str2double(answ(1));
        bdwidth=str2double(answ(2));
      end
      [sep, Y, seps] = icainf(X', eye(2), maxiter, bdwidth);

% Syntax:  [sep, source, seps] = icainf(data, sep, maxiter, bdwidth, nb)
%
% Separation of K linear instantaneous mixture of K sources based on
% the minimisation of the mutual information. The entropy is estimated
% is estimated via a kernel density estimator with the kernel being
% the density of the sum of 3 independent uniform random variables
% in [-1/2,1/2]. The kernel bandwidth is set to bdwidth times the
% optimal bandwith for Gausian density (= 2*(11*sqrt(pi)/(12*n))^.2 =
% 2.107683/n^.2) and the entropy is computed by discretizing at a step
% bandwidth/nb. 
%
% data must be a n by K matrix and sep a K by K matrix, sep, maxiter,
% bdwidth, nb are optional:
% sep at input is the the initial separation matrix (defaults to the
% identity matrix) and at output is the obtained separation matrix
% maxiter defaults to 15,
% bdwidth default 1
% nb defaults to 1
% seps, if present, is a iter*K by K matrix containing the sequence
% of iterated separation matrices, including the initial one.

      Y = Y';

    case 4
      prompt={'maximum number of iterations: '
	      'Stationary block length of the data:'};
      def={'40','100'};
      dlgTitle='Input for the parameters for the algorithm';
      answ=inputdlg(prompt,dlgTitle,1,def);
      if ~isempty(answ)
        maxiter=str2double(answ(1));
        bloclen=str2double(answ(2));
      end
      nbloc = round(1000/bloclen);
      [sep, Y, seps] = icansng(X, eye(2), nbloc, maxiter);

% Syntaxe       [sep, source, seps] = icansng(data, sep, nbloc, maxiter)
% Separation of non stationary non gaussian sources, batch method
% Use linear combinations of s^.5, s, s^3 as separating function
%
%  - data and source are m x n matrices, containing m mixtures and
%    m reconstructed sources.
%  - sep in input is the initial separating matrix
%    (default = identity matrix).
%  - sep in output is the estimated separating matrix.
%  - nbloc is the number of blocks in each of which the sources are
%    considered to be stationary (default = 3)
%  - maxiter = number of iteration
%  - eps: (squared) stoping threshold, default 1e-2*K(K-1)/n

    case 5
      prompt = {'Number of parameters for the nonlin. transf.'
		'Maximum number of iterations'
		'kernel binwidth/(standard value)'
		'learning step for the linear part'
		'learning step for the nonlinear part'};
      def={'12','40','1','0.6','0.3'};
      dlgTitle='Input for the parameters for the algorithm';
      answ=inputdlg(prompt,dlgTitle,1,def);
      if ~isempty(answ)
	np = str2double(answ(1));
        maxiter = str2double(answ(2));
        bin = str2double(answ(3));
	mu_B = str2double(answ(4));
	mu_z =  str2double(answ(5));
      end
      [Y, Z, sep, seps] = icainfpnl(X',np,maxiter,bin,mu_B,mu_z);
      Y = Y';
      axes(handles.Axes2);
      plot(sort(X0(1,:)),sort(Z(:,1)));
      axis tight; grid on
      axes(handles.Axes3);
      plot(sort(X0(2,:)),sort(Z(:,2)));
      axis tight; grid on
      set(handles.PostnonlinText,'String', ...
	  'Compenstators of post nonlinear transform');
      set(handles.text4,'string','X2 vs X1');
      axes(handles.axes9);
      plot(X(1,:),X(2,:), '.'); axis tight;
      set(handles.text3,'string','Compensated X1 versus X2 (ICA axes in red)');
      axes(handles.axes8);
      plot(Z(:,1),Z(:,2),'.'); axis tight;

    otherwise
      return
  end
end

if ~isempty(answ)

  if algo < 5
    set(handles.text3,'string','Y2 vs Y1');
    axes(handles.axes8);
    plot(Y(1,:),Y(2,:), '.'); axis tight;
    set(handles.PostnonlinText,'String','Post nonlinear transform functions');
    axes(handles.Axes2);
    plot(sort(X0(1,:)),sort(X(1,:)));
    axis tight; grid on;
    axes(handles.Axes3);
    plot(sort(X0(2,:)),sort(X(2,:)));
    axis tight; grid on
    set(handles.text4,'string','X2 vs X1 (with ICA axes in red)');
    axes(handles.axes9);
    plot(X(1,:),X(2,:), '.'); axis tight;
  end
  A = inv(sep);
  mini = min(Y'); maxi = max(Y');
  lh = line([[mini(1) 0; maxi(1) 0]*A(1,:)' ...
	     [0 mini(2); 0 maxi(2)]*A(1,:)'], ...
	    [[mini(1) 0; maxi(1) 0]*A(2,:)' ...
	     [0 mini(2); 0 maxi(2)]*A(2,:)']);
  set(lh(1),'Color', 'red');
  set(lh(2),'Color', 'red');

  axes(handles.axes5);
  plot(Y(1,:));
  %title('Estimated Source 1');
  %set(handles.text9,'string','Estimated Source 1 (Y1)');

  axes(handles.axes6);
  plot(Y(2,:));
  %title('Estimated Source 2');
  %set(handles.text5,'string','Estimated Source 2 (Y2)');

  axes(handles.axes7);
  iter = size(seps,1)/2;
  seps = reshape((seps*[1 str2double(get(handles.m12,'String'))
			str2double(get(handles.m21,'String')) 1])',4,iter);
  %plot(0:iter-1, seps);
  plot(0:iter-1, seps(1:2,:), 'b', 0:iter-1, seps(3:4,:), 'r');
  grid on;
  %set(handles.text10,'string','Global matrix: blue=row1, red=row2');
  xlabel('iteration step');

end
