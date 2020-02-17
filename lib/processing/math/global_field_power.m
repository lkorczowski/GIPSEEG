function [gfp]=global_field_power(P,param)
%%
% *** History: 19-Mar-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%close all
for k=1:size(P,3);
u=P(:,:,k)';
%{
Ne=size(P,1);
Nt=size(P,2);
sumsqdif = zeros(Nt,1);
for i=1:Ne,
    
    %progress = sprintf('...electrode %6d of %6d\n',i,Ne);
    %if i>1, progress = [repmat('\b',1,length(progress)),progress]; end
    %fprintf(progress);
    
    % each iteration, extract an electrode waveform and then remove
    % it from further consideration, so
    % a) extract electrode i
    ui = u(:,i);
    ui = repmat(ui,1,length(i:Ne));
    % b) only consider this and the remaining electrodes
    uj = u(:,i:Ne);
    % differences of this waveform with all others (the
    % difference with itself is zero, it will not add to the sum)
    sqdif = (ui - uj).^2;
    % take the sum over electrodes
    sqdif = sum( sqdif, 2);
    % cumulative pairwise differences
    sumsqdif = sum([sumsqdif, sqdif], 2);
end
gfp_de(:,k)=sqrt(sumsqdif);

%}
gfp(:,k)=sqrt(mean(u.^2,2));
end
