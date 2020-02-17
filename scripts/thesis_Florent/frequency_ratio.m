function [scoreA scoreB]=frequency_ratio(Xtilde_fk)
%score=frequency_ratio(Xtilde_fk)
% compute the frequency ratio of Xtilde_fk
% INPUTS
% Xf is a 3D matrix : [N x F x K]
%                     N (electrodes)
%                     F (frequencies)
%                     K (epochs)
%
% OUTPUTS
% scoreA [F x 1] is the ratio of averaged frequency (between 0 and 1)
% scoreB [F x K] is the frequency ratio for each frequency
%                 and each epoch (between 0 and 1)

N=size(Xtilde_fk,1);
F=size(Xtilde_fk,2);
K=size(Xtilde_fk,3);

scoreA=NaN(F,1);
scoreB=NaN(F,K);
frac_den=NaN(K,1);
frac_num=NaN(F,K);
%compute denominator (total power of each epoch)
for k=1:size(Xtilde_fk,3)
     frac_den(k)=norm(Xtilde_fk(:,:,k),'fro');
end

%compute numerator and scoreB
for  k=1:size(Xtilde_fk,3)
     for f=1:size(Xtilde_fk,2)
         frac_num(f,k)=norm(Xtilde_fk(:,f,k),'fro');
         scoreB(f,k)=frac_num(f,k)/frac_den(k);
     end
end

     for f=1:size(Xtilde_fk,2)
scoreA(f)=sum(frac_num(f,k),2)/sum(frac_den(k));
     end