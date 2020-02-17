%% Generation of CFTf surrogates - useful to test the DC - VERSION FOR EXTENDED MVAR MODEL!
%%% ABSENCE OF FULL (direct&indirect) CAUSALITY FROM yj TO yi - Bim and Bmj coefficients are forced to zero
% Note: the length of the series should be not be odd

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% Bm: extended MVAR coefficients
% B0: instantaneous effects
% Sw: input covariance matrix of extended model
% jj: input series
% ii: output series
% flag: 'ExtendedGauss', 'ExtendedNonGauss' - to say whether identification has been performed with prior assumptions (use nongaussian residuals) or not (use gaussian residuals)

%%% output:
% Ys, M*N matrix of surrogate time series (each time series is in a row)

function Ys=surrVCFTf0(Y,Bm,B0,Sw,ii,jj)

Y=Y'; %transpose series (the following code works with column vectors)
[N,M]=size(Y);
p = size(Bm,2)/M; % p is the order of the MVAR model
Ys=NaN*ones(N,M);
for m=1:M %zero mean
    Y(:,m)=Y(:,m)-mean(Y(:,m));
end

% force to zero the coefficients Bim(r), for each r=0,1,..,p and each m~=i: blocks inflow into ii channel
for mm=1:M
    if mm ~= ii
        B0(ii,mm)=0; % r=0
        for r=1:p %r>=1
            Bm(ii,(r-1)*M+mm)=0;
        end
    end
end
% force to zero the coefficients Bmj(r), for each r=0,1,..,p and each m~=j: blocks outflow from jj channel
for mm=1:M
    if mm ~= jj
        B0(mm,jj)=0; % r=0
        for r=1:p %r>=1
            Bm(mm,(r-1)*M+jj)=0;
        end
    end
end

% Fitted series: feed an extended MVAR model, formed by the new coeffs, with noise with variance Sw
Am=diag_coeff_rev(Bm,B0,Sw); %strictly causal coeffs from extended coeffs
U=InstModelfilter(N,Sw,'ExtendedGauss',B0); % innovations W (gaussian or not depending on flag), U determined from B0
[Yc]=MVARfilter(Am,U);

Yc=Yc'; %column Yc

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
