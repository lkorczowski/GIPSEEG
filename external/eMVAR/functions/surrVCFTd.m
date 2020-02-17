%% Generation of CFTd surrogates - useful to test the PDC - VERSION FOR STRICTLY CAUSAL MVAR MODEL
%%% ABSENCE OF DIRECT CAUSALITY FROM yj TO yi - Aij coefficients are forced to zero
% Note: the length of the series should be not be odd

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% Am: strictly causal MVAR coefficients
% Su: input covariance matrix of strictly causal model
% jj: input series
% ii: output series

%%% output:
% Ys, M*N matrix of surrogate time series (each time series is in a row)

function Ys=surrVCFTd(Y,Am,Su,ii,jj)

Y=Y'; %transpose series (the following code works with column vectors)
[N,M]=size(Y);
p = size(Am,2)/M; % p is the order of the MVAR model
Ys=NaN*ones(N,M);
for m=1:M%zero mean
    Y(:,m)=Y(:,m)-mean(Y(:,m));
end

% force to zero the coefficients Aij(r), for each r=1,..,p
for r=1:p %r>=1
    Am(ii,(r-1)*M+jj)=0;
end

% Fitted series: feed a strictly causal MVAR model, formed by the new coeffs, with gaussian noise with variance Su
U=InstModelfilter(N,Su,'StrictlyCausal'); % generate gaussian innovations with covariance Su
[Yc]=MVARfilter(Am,U); % use U as input for the MVAR model with some Am coeffs reduced to zero

Yc=Yc';

for cnt=1:M
    % phase from fittated series, channel cnt = Yc(:,cnt)
    xc=Yc(:,cnt);
    fxc=fft(xc);
    fasexc=angle(fxc);
    % modulus from original series, channel cnt =  Y(:,cnt)
    x=Y(:,cnt);
    fx=fft(x);
    mx=abs(fx);
    % Fourier Transform of surrogate series:
    fxs=mx.*(cos(fasexc)+j*sin(fasexc));
    % Surrogate series:
    xs=ifft(fxs);xs=real(xs);
% % %     xs=xs-mean(xs);
    Ys(:,cnt)=xs;
end
    
Ys=Ys'; % return to row vectors
    
    

