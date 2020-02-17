function [GradPSD R C Pval]=gradient_psd(PSD,freqs)

GradPSD=[];
for n=1:size(PSD,1)
[Rtmp,P]  =corrcoef(PSD(n,:),freqs);
Rtmp=Rtmp(2);
Pval(n)=P(2);
C(n)=(freqs-mean(freqs))'\(PSD(n,:)-mean(PSD(n,:)))';
end
