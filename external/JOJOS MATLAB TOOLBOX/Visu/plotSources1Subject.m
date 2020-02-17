function plotSources1Subject(B,Cosp,NumIniFreq,NumFinalFreq,faxis,FreqShow,NumSources)
% PlotSources1Subject(B,Cosp,NumIniFreq,NumFinalFreq,faxis,FreqShow,NumSources)

% Description:
% ------------
% This function plot the energy of each source and each condition from 
% cospectra matrices computed previously (apply product B*Cosp*B')
% 
% Inputs:
% -------
% B             --> Demixing Matrix
% Cosp          --> Cosp we want to compute the power of the sources
% NumIniFreq    --> Index to start to compute the power spectra
% NumFinalFreq  --> Index to end to compute the power spectra
% faxis         --> array with the positions of the frequencies
% FreqShow      --> index with the postiton to start to plot
% NumSources    --> num of the sources we want to plot
%
% History:
% --------
% Last version: 2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
   

 error(nargchk(6,7,nargin));         % Assertion of inputs
 if(nargin < 7)
     NumSources = size(B,1);
 end

 %Calculate the diagonal elements of the product B*C*B'
 for f=NumIniFreq:NumFinalFreq
    DiagElems(:,:,f)=diag(B*Cosp(:,:,f)*B');
 end

 % Erase the dimension that is equal to 1
 % Result is M x f where M num components; f num freqs and i num trials 
  Energy=squeeze(DiagElems);
  Energy=Energy';

  % Plot the sources
  [x, y] = num_windows(NumSources);
  for m=1:NumSources
    subplot(x,y,m);
    plot (faxis(FreqShow:NumFinalFreq),Energy(FreqShow:NumFinalFreq,m),'blue');
    title (sprintf('Source %d',m));
    hold off
    axis tight;
  end